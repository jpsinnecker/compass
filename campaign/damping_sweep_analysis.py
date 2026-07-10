#!/usr/bin/env python3
"""
damping_sweep_analysis.py
==========================

Post-processing pipeline for the damping/hysteresis study. Reads the
results produced by damping_sweep.py (manifest.csv + per-run CSVs) and
produces:

  1. For each hysteresis run: loop area, coercive field Bc, remanence Mr.
  2. Avalanche statistics: event detection via dM/dB within each field
     ramp segment, and power-law fit P(s) ~ s^{-tau} by maximum-likelihood
     estimation (Clauset-Shalizi-Newman method), with x_min selection by
     KS-distance minimisation.
  3. For each free-relaxation run (Stage 2): domain-size statistics
     (label_magnetic_domains, computed during simulation and saved in the
     meta JSON).
  4. Final plots:
       - area_vs_Q.png, Bc_vs_Q.png, Mr_vs_Q.png  (one curve per geometry)
       - tau_vs_Q.png                              (central test)
       - avalanche_histograms.png                  (P(s) per geometry, selected Q)
       - domain_size_vs_Q.png                      (mean domain size / count vs Q)
  5. A consolidated summary.csv (one row per hysteresis run, all metrics),
     for downstream inspection or tabulation.

Usage
-----
    python3 damping_sweep_analysis.py --sweep_dir /home/jps/sweep_results \
        --out_dir /home/jps/sweep_analysis

Dependencies: numpy, pandas, matplotlib, scipy.
Optional: powerlaw package (if absent, a built-in MLE equivalent to the
Clauset-Shalizi-Newman continuous method is used, with x_min by KS).
"""

import os
import sys
import json
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import powerlaw as _powerlaw_pkg
    _HAVE_POWERLAW = True
except ImportError:
    _HAVE_POWERLAW = False


def _bar(done, total, label="", extra="", width=36):
    """One-line overwriting progress bar."""
    frac   = done / total if total > 0 else 0.0
    filled = int(width * frac)
    bar    = "█" * filled + "░" * (width - filled)
    pct    = frac * 100.0
    suffix = f"  {extra}" if extra else ""
    line   = f"  {label:10s} [{bar}] {done:>4d}/{total}  {pct:5.1f}%{suffix}"
    end    = "\n" if done == total else "\r"
    print(line, end=end, flush=True)


# ────────────────────────────────────────────────────────────────────────
# Estilo de plot
# ────────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "font.size": 11,
})

GEOM_COLORS = {
    "square":     "#2E86AB",
    "triangular": "#E94560",
    "honeycomb":  "#3CB371",
}
GEOM_MARKERS = {
    "square":     "o",
    "triangular": "s",
    "honeycomb":  "^",
}


# ════════════════════════════════════════════════════════════════════════
# 1. Hysteresis loop metrics: area, coercivity, remanence
# ════════════════════════════════════════════════════════════════════════

def loop_area(B, M):
    """Area of the closed hysteresis loop via line integral
    (shoelace formula applied to the closed curve in the B-M plane).
    Equivalent to the integral of M dB over the full cycle."""
    B = np.asarray(B)
    M = np.asarray(M)
    # ensure the cycle is closed (last point ~ first point);
    # if not exactly closed (numerical noise on last step), close explicitly
    # so the line integral does not leak spurious area.
    if not np.isclose(B[0], B[-1], atol=1e-9) or not np.isclose(M[0], M[-1], atol=1e-6):
        B = np.append(B, B[0])
        M = np.append(M, M[0])
    # area of closed polygon in the (B, M) plane:
    area = 0.5 * np.abs(np.sum(B[:-1] * M[1:] - B[1:] * M[:-1]))
    return area


def find_zero_crossings(x, y):
    """Returns the values of x where y crosses zero, via linear interpolation
    between consecutive points of opposite sign."""
    sign = np.sign(y)
    sign[sign == 0] = 1  # avoid spurious crossings where y is exactly 0
    idx = np.where(np.diff(sign) != 0)[0]
    crossings = []
    for i in idx:
        x0, x1 = x[i], x[i + 1]
        y0, y1 = y[i], y[i + 1]
        if y1 == y0:
            continue
        xc = x0 + (0.0 - y0) * (x1 - x0) / (y1 - y0)
        crossings.append(xc)
    return np.array(crossings)


def coercivity_and_remanence(B, M):
    """
    Bc (coercive field)  : value of |B| where M crosses zero (median of
                            all crossings found, since the loop has more
                            than one crossing over the full cycle).
    Mr (remanence)        : value of |M| where B crosses zero.
    """
    B = np.asarray(B)
    M = np.asarray(M)

    # Bc: where M(t) crosses zero -> read the corresponding B value
    M_crossings_idx_B = find_zero_crossings(B, M)  # x=B, y=M -> returns B where M=0
    # find_zero_crossings(x,y) finds x where y=0. For Bc we want B where M=0
    # -> x must be B, y must be M. Already correct above.
    Bc_candidates = np.abs(M_crossings_idx_B)
    Bc = np.median(Bc_candidates) if len(Bc_candidates) > 0 else np.nan

    # Mr: where B(t) crosses zero -> read the corresponding M value. Using B
    # directly as interpolation axis doesn't work (B crosses zero multiple times);
    # instead interpolate M as a function of index/time at the points where B=0.
    B_crossings_idx = find_zero_crossings(np.arange(len(B)), B)
    if len(B_crossings_idx) > 0:
        M_at_Bzero = np.interp(B_crossings_idx, np.arange(len(M)), M)
        Mr = np.median(np.abs(M_at_Bzero))
    else:
        Mr = np.nan

    return Bc, Mr


# ════════════════════════════════════════════════════════════════════════
# 2. Avalanche detection via dM/dB within each field-ramp segment
# ════════════════════════════════════════════════════════════════════════

def split_hysteresis_segments(B):
    """
    The hysteresis cycle has 5 linear segments (0->+Bmax->0->-Bmax->0->+Bmax).
    Detects inflection points (where dB/dt changes sign) to split the series
    into monotone ramp segments -- avalanches are detected WITHIN each segment,
    never across an inflection (where dB/dt passes through zero and dM/dB would
    diverge numerically without being a real avalanche).
    """
    dB = np.diff(B)
    sign_dB = np.sign(dB)
    sign_dB[sign_dB == 0] = 0  # stationary B points (rare, numerical equilibrium)
    change_idx = np.where(np.diff(sign_dB) != 0)[0] + 1
    bounds = [0] + list(change_idx + 1) + [len(B)]
    bounds = sorted(set(b for b in bounds if 0 <= b <= len(B)))
    segments = [(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)
                if bounds[i + 1] - bounds[i] >= 5]  # skip degenerate segments (< 5 points)
    return segments


def detect_avalanches(t, B, M, mad_threshold=4.0, min_event_len=1):
    """
    For each monotone field ramp segment, computes dM/dB and marks as
    "avalanche" the steps where |dM/dB - local_median| exceeds
    mad_threshold * MAD (median absolute deviation, robust to outliers --
    appropriate here because the avalanche size distribution itself is
    heavy-tailed, which would inflate a standard deviation).

    Contiguous events (consecutive steps above threshold) are grouped
    into a single avalanche; the "size" of the avalanche is the sum of |dM|
    absorbed during the event (proxy for the magnetisation reversed
    in that avalanche).

    Returns: array of avalanche sizes (one value per detected event,
    aggregating all segments of the cycle).
    """
    segments = split_hysteresis_segments(B)
    sizes = []

    for (i0, i1) in segments:
        Bs = B[i0:i1]
        Ms = M[i0:i1]
        dB_seg = np.diff(Bs)
        dM_seg = np.diff(Ms)
        # avoid division by dB ~ 0 (should not occur inside a valid
        # monotone segment, but guards against numerical noise)
        valid = np.abs(dB_seg) > 1e-14
        if valid.sum() < 5:
            continue
        dMdB = np.full_like(dB_seg, np.nan)
        dMdB[valid] = dM_seg[valid] / dB_seg[valid]

        local = dMdB[valid]
        med = np.median(local)
        mad = np.median(np.abs(local - med)) + 1e-300  # guard against exact MAD=0
        excess = np.zeros_like(dMdB)
        excess[valid] = np.abs(dMdB[valid] - med) / mad

        is_event = excess > mad_threshold

        # agrupa eventos contiguos
        idx = np.where(is_event)[0]
        if len(idx) == 0:
            continue
        groups = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
        for g in groups:
            if len(g) < min_event_len:
                continue
            # avalanche size = sum of |dM| absorbed in the event
            size = np.sum(np.abs(dM_seg[g]))
            if size > 0:
                sizes.append(size)

    return np.array(sizes)


def detect_avalanches_from_S(t, S, min_event_len=1):
    """
    Detects avalanches using the 'S' column of field_log, which is the direct
    count of magnetic-moment flips (spin flips) recorded by compass.py at each
    integration step.  An avalanche is a maximal run of consecutive steps with
    S > 0; its size is the sum of S over the run.

    This is a cleaner alternative to the dM/dB proxy because:
      - S is an integer spin-flip count (no numerical differentiation noise).
      - It is not confused by numerical drift in B near inflection points.
      - It directly maps to the physical picture of discrete needle reversals.

    Returns (sizes, inter_arrival_times):
      sizes: array of avalanche sizes (sum of S per event)
      inter_arrival_times: array of waiting times between consecutive
          avalanche *starts* [seconds]
    """
    S = np.asarray(S)
    t = np.asarray(t)
    is_active = S > 0
    # find runs of True (S > 0)
    changes = np.diff(is_active.astype(int), prepend=0, append=0)
    starts = np.where(changes == 1)[0]
    ends   = np.where(changes == -1)[0]
    sizes = []
    event_times = []
    for s, e in zip(starts, ends):
        if (e - s) < min_event_len:
            continue
        sizes.append(int(S[s:e].sum()))
        event_times.append(float(t[s]))
    sizes = np.array(sizes, dtype=float)
    if len(event_times) > 1:
        iat = np.diff(event_times)
    else:
        iat = np.array([])
    return sizes, iat


def waiting_time_stats(inter_arrival_times):
    """
    Characterises the waiting-time distribution between avalanches.
    Fits two models:
      1. Exponential (Poisson process): P(tau) = lambda * exp(-lambda*tau)
         MLE: lambda_hat = 1/mean(tau)
      2. Power-law (scale-free dynamics):  P(tau) ~ tau^-alpha
         MLE via _mle_continuous_powerlaw.

    Returns dict with:
      n, mean_iat, cv (coefficient of variation = std/mean),
      lambda_exp (Poisson rate), log_likelihood_exp,
      alpha_pl, xmin_pl, log_likelihood_pl,
      lrt_statistic (2*(LL_pl - LL_exp)), lrt_favors_powerlaw (bool).

    High CV (>> 1) and a positive LRT statistic indicate bursty / scale-free
    dynamics (SOC-like); CV ~ 1 and negative LRT indicate Poisson.
    """
    tau = np.asarray(inter_arrival_times)
    tau = tau[tau > 0]
    n = len(tau)
    if n < 10:
        return dict(n=n, mean_iat=np.nan, cv=np.nan,
                    lambda_exp=np.nan, log_likelihood_exp=np.nan,
                    alpha_pl=np.nan, xmin_pl=np.nan, log_likelihood_pl=np.nan,
                    lrt_statistic=np.nan, lrt_favors_powerlaw=False)
    mean_tau = float(np.mean(tau))
    cv = float(np.std(tau) / mean_tau)
    # Exponential MLE
    lam = 1.0 / mean_tau
    ll_exp = float(n * np.log(lam) - lam * np.sum(tau))
    # Power-law MLE (use same fallback as fit_powerlaw_mle)
    pl = fit_powerlaw_mle(tau)
    alpha = pl["tau"]          # power-law alpha for waiting times
    xmin  = pl["xmin"]
    if not np.isnan(alpha) and xmin > 0:
        tail = tau[tau >= xmin]
        ll_pl = float(len(tail) * (np.log(alpha - 1) - np.log(xmin))
                      - alpha * np.sum(np.log(tail / xmin)))
    else:
        ll_pl = np.nan
    lrt = 2.0 * (ll_pl - ll_exp) if not np.isnan(ll_pl) else np.nan
    return dict(
        n=n, mean_iat=mean_tau, cv=cv,
        lambda_exp=lam, log_likelihood_exp=ll_exp,
        alpha_pl=alpha, xmin_pl=xmin, log_likelihood_pl=ll_pl,
        lrt_statistic=lrt, lrt_favors_powerlaw=bool(not np.isnan(lrt) and lrt > 0),
    )


def autocorrelation_M(t, M, max_lag_fraction=0.25):
    """
    Computes the temporal autocorrelation function C(lag) = <M(t)M(t+lag)> / var(M)
    for lags up to max_lag_fraction * total_duration.

    Returns (lags_t, C) and the estimated correlation time tau_corr (area
    under C until it first crosses zero, or the full range if it never crosses).

    Motivation: the 2025 echo-state paper (s41598-025-93189-w) shows that the
    memory capacity of ASI is a strong function of damping Q.  tau_corr vs Q
    is therefore a direct experimental prediction of that framework.
    """
    M = np.asarray(M, dtype=float)
    t = np.asarray(t, dtype=float)
    n = len(M)
    M_c = M - np.mean(M)
    var = float(np.var(M_c))
    if var < 1e-30 or n < 20:
        return np.array([]), np.array([]), np.nan

    # Use only the first max_lag_fraction fraction of the time series as lags
    dt_mean = (t[-1] - t[0]) / (n - 1)
    max_lag = max(1, int(max_lag_fraction * n))

    C = np.zeros(max_lag)
    for k in range(max_lag):
        C[k] = float(np.mean(M_c[:n - k] * M_c[k:])) / var

    lags_t = np.arange(max_lag) * dt_mean

    # Correlation time: integral up to first zero-crossing (or full range)
    zero_crossings = np.where(C < 0)[0]
    if len(zero_crossings) > 0:
        cutoff_idx = zero_crossings[0]
    else:
        cutoff_idx = max_lag
    tau_corr = float(np.trapz(C[:cutoff_idx], lags_t[:cutoff_idx]))
    return lags_t, C, tau_corr


def lrt_powerlaw_vs_lognormal(sizes):
    """
    Likelihood-ratio test: power-law vs log-normal, following Clauset, Shalizi
    & Newman 2009 §4.  Both distributions are fit to x >= xmin (the xmin chosen
    by the power-law MLE).

    Returns dict with:
      R (log-likelihood ratio; R>0 favors power-law, R<0 favors log-normal),
      p (two-sided p-value, ~significance of the sign of R; computed by a
         Vuong-style normalization using the standard error of R),
      favors_powerlaw (bool: R > 0 AND p < 0.1).
    """
    sizes = np.asarray(sizes)
    sizes = sizes[sizes > 0]
    if len(sizes) < 20:
        return dict(R=np.nan, p=np.nan, favors_powerlaw=False)

    pl = fit_powerlaw_mle(sizes)
    if np.isnan(pl["tau"]):
        return dict(R=np.nan, p=np.nan, favors_powerlaw=False)

    xmin = pl["xmin"]
    xs = sizes[sizes >= xmin]
    n = len(xs)
    if n < 10:
        return dict(R=np.nan, p=np.nan, favors_powerlaw=False)

    # Power-law log-likelihoods
    alpha = pl["tau"]
    ll_pl_i = np.log(alpha - 1.0) - np.log(xmin) - alpha * np.log(xs / xmin)

    # Log-normal MLE on xs
    log_xs = np.log(xs)
    mu_ln  = float(np.mean(log_xs))
    sig_ln = float(np.std(log_xs, ddof=1))
    if sig_ln < 1e-14:
        return dict(R=np.nan, p=np.nan, favors_powerlaw=False)
    from scipy.stats import norm as _norm
    ll_ln_i = _norm.logpdf(log_xs, loc=mu_ln, scale=sig_ln) - log_xs  # log-normal PDF

    log_ratio_i = ll_pl_i - ll_ln_i
    R = float(np.sum(log_ratio_i))
    se = float(np.std(log_ratio_i) * np.sqrt(n))
    if se < 1e-14:
        return dict(R=R, p=np.nan, favors_powerlaw=(R > 0))
    from scipy.special import erfc as _erfc
    p = float(_erfc(abs(R / se) / np.sqrt(2.0)))
    return dict(R=R, p=p, favors_powerlaw=bool(R > 0 and p < 0.1))


# ════════════════════════════════════════════════════════════════════════
# 3. Power-law fitting: MLE (Clauset-Shalizi-Newman method)
# ════════════════════════════════════════════════════════════════════════

def _mle_continuous_powerlaw(x, xmin):
    """
    Estimador de maxima verossimilhanca para o expoente tau de uma lei de
    potencia continua P(x) ~ x^-tau, x >= xmin (Clauset, Shalizi & Newman
    2009, eq. 3.1):

        tau_hat = 1 + n * [ sum( ln(x_i / xmin) ) ]^-1

    Valido para x continuo. Retorna (tau_hat, erro_padrao).
    """
    xs = x[x >= xmin]
    n = len(xs)
    if n < 10:
        return np.nan, np.nan, n
    tau_hat = 1.0 + n / np.sum(np.log(xs / xmin))
    sigma = (tau_hat - 1.0) / np.sqrt(n)
    return tau_hat, sigma, n


def _ks_distance_powerlaw(x, xmin, tau):
    """Distancia de Kolmogorov-Smirnov entre a CDF empirica de x (x>=xmin)
    e a CDF teorica de uma lei de potencia continua com expoente tau."""
    xs = np.sort(x[x >= xmin])
    n = len(xs)
    if n < 2:
        return np.inf
    cdf_emp = np.arange(1, n + 1) / n
    cdf_theo = 1.0 - (xs / xmin) ** (1.0 - tau)
    return np.max(np.abs(cdf_emp - cdf_theo))


def fit_powerlaw_mle(sizes, xmin_candidates=None):
    """
    Fits a power law to avalanche sizes via MLE, selecting
    x_min pela minimizacao da distancia KS (metodo padrao de Clauset,
    Shalizi & Newman 2009 -- "Power-law distributions in empirical data").

    Usa o pacote `powerlaw` se disponivel (implementacao de referencia);
    caso contrario, usa a implementacao propria acima, que segue a mesma
    formula MLE e o mesmo criterio de selecao de x_min por KS.

    Retorna dict com: tau, tau_err, xmin, ks, n_tail (pontos na cauda
    ajustada), n_total (pontos antes do corte).
    """
    sizes = np.asarray(sizes)
    sizes = sizes[sizes > 0]
    n_total = len(sizes)
    if n_total < 20:
        return dict(tau=np.nan, tau_err=np.nan, xmin=np.nan, ks=np.nan,
                     n_tail=0, n_total=n_total)

    if _HAVE_POWERLAW:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit = _powerlaw_pkg.Fit(sizes, discrete=False, verbose=False)
        tau = fit.power_law.alpha
        tau_err = fit.power_law.sigma
        xmin = fit.power_law.xmin
        # NOTA: em algumas versoes do pacote `powerlaw`, o metodo .KS() esta
        # quebrado (bug conhecido: chama uma funcao de modulo sem o prefixo
        # `self.`, lancando NameError). O atributo .D guarda a mesma
        # distancia KS, ja calculada internamente durante Fit() na escolha
        # de xmin -- usamos ele diretamente, com fallback para .KS() caso
        # uma versao futura remova .D.
        try:
            ks = fit.power_law.D
        except AttributeError:
            try:
                ks = fit.power_law.KS()
            except Exception:
                ks = np.nan
        n_tail = int(np.sum(sizes >= xmin))
        return dict(tau=tau, tau_err=tau_err, xmin=xmin, ks=ks,
                     n_tail=n_tail, n_total=n_total)

    # ── fallback proprio (mesma metodologia, sem dependencia externa) ──
    if xmin_candidates is None:
        # candidatos: quantis da propria amostra (evita testar todo ponto
        # unico quando n_total e grande, o que seria O(n^2) sem necessidade)
        xmin_candidates = np.unique(np.quantile(sizes, np.linspace(0.0, 0.9, 40)))
        xmin_candidates = xmin_candidates[xmin_candidates > 0]

    best = dict(tau=np.nan, tau_err=np.nan, xmin=np.nan, ks=np.inf, n_tail=0)
    for xmin in xmin_candidates:
        tau_hat, sigma, n_tail = _mle_continuous_powerlaw(sizes, xmin)
        if np.isnan(tau_hat) or n_tail < 15:
            continue
        ks = _ks_distance_powerlaw(sizes, xmin, tau_hat)
        if ks < best["ks"]:
            best = dict(tau=tau_hat, tau_err=sigma, xmin=xmin, ks=ks, n_tail=n_tail)

    best["n_total"] = n_total
    return best


# ════════════════════════════════════════════════════════════════════════
# 4. Main pipeline: reads manifest, processes each run, aggregates
# ════════════════════════════════════════════════════════════════════════

def process_hysteresis_run(row, mad_threshold):
    df = pd.read_csv(row["csv_path"])
    B, M = df["B"].values, df["M"].values
    t = df["t"].values
    # The 'S' column (spin-flip count) may not exist in older CSVs
    S = df["S"].values if "S" in df.columns else np.zeros_like(t)

    area = loop_area(B, M)
    Bc, Mr = coercivity_and_remanence(B, M)

    # Method 1: dM/dB proxy (MAD-based)
    avalanche_sizes_dmdB = detect_avalanches(t, B, M, mad_threshold=mad_threshold)

    # Method 2: direct spin-flip counter (S column) — cleaner, no differentiation
    avalanche_sizes_S, inter_arrival_times = detect_avalanches_from_S(t, S)

    # Waiting-time statistics (from S-column events)
    wt_stats = waiting_time_stats(inter_arrival_times)

    # Temporal autocorrelation of M(t)
    _, _, tau_corr = autocorrelation_M(t, M)

    with open(row["meta_path"]) as f:
        meta = json.load(f)

    return dict(
        tag=row["tag"], geometry=row["geometry"], damp_idx=row["damp_idx"],
        seed=row["seed"], damping=row["damping"], Q=row["Q"],
        area=area, Bc=Bc, Mr=Mr,
        n_avalanches=len(avalanche_sizes_dmdB),
        avalanche_sizes=avalanche_sizes_dmdB,
        n_avalanches_S=len(avalanche_sizes_S),
        avalanche_sizes_S=avalanche_sizes_S,
        inter_arrival_times=inter_arrival_times,
        wt_mean_iat=wt_stats["mean_iat"],
        wt_cv=wt_stats["cv"],
        wt_lrt=wt_stats["lrt_statistic"],
        wt_favors_powerlaw=wt_stats["lrt_favors_powerlaw"],
        tau_corr=tau_corr,
        B_max=meta.get("B_max", np.nan),
    )


def fit_powerlaw_per_geom_damp(summary_df, size_col="avalanche_sizes"):
    """
    Aggregates avalanche sizes from ALL seeds of each (geometry, damp_idx)
    (geometria, damp_idx) antes de ajustar a lei de potencia -- um unico
    ajuste MLE por ponto (geometria, Q), em vez de um ajuste fraco por
    individual run (which typically has few avalanches and produces unstable tau
    instavel / NaN por falta de amostra).

    Adds a LRT column comparing power-law vs log-normal on the aggregated
    sample (Clauset, Shalizi & Newman 2009 §4).

    Retorna um DataFrame com uma linha por (geometria, damp_idx).
    """
    rows = []
    for (geom, di), grp in summary_df.groupby(["geometry", "damp_idx"]):
        all_sizes = np.concatenate(grp[size_col].values) \
            if len(grp) else np.array([])
        fit = fit_powerlaw_mle(all_sizes)
        lrt = lrt_powerlaw_vs_lognormal(all_sizes)
        rows.append(dict(
            geometry=geom, damp_idx=di, Q_mean=grp["Q"].mean(),
            tau=fit["tau"], tau_err=fit["tau_err"], xmin=fit["xmin"],
            ks=fit["ks"], n_tail=fit["n_tail"], n_total=fit.get("n_total", 0),
            lrt_R=lrt["R"], lrt_p=lrt["p"],
            lrt_favors_powerlaw=lrt["favors_powerlaw"],
            n_seeds=len(grp),
        ))
    return pd.DataFrame(rows).sort_values(["geometry", "damp_idx"])


def process_relax_run(row):
    with open(row["meta_path"]) as f:
        meta = json.load(f)
    sizes = np.array(meta.get("domain_sizes", []))
    sizes = sizes[sizes > 0]
    if len(sizes) == 0:
        mean_size, max_size, n_dom = np.nan, np.nan, 0
    else:
        mean_size = float(np.mean(sizes))
        max_size = float(np.max(sizes))
        n_dom = len(sizes)
    return dict(
        tag=row["tag"], geometry=row["geometry"], damp_idx=row["damp_idx"],
        seed=row["seed"], damping=row["damping"], Q=row["Q"],
        n_domains=n_dom, mean_domain_size=mean_size, max_domain_size=max_size,
    )


def aggregate_by_geom_damp(df, value_cols):
    """Agrega (media, desvio padrao) sobre seeds, para cada (geometria, damp_idx),
    preservando Q (que e ~constante entre seeds, ja que so depende de
    geometria/damping/r_nn -- pequenas variacoes de r_nn entre geometrias
    are expected and the Q of each individual run is already correct per row;
    aqui usamos a media de Q tambem, por consistencia)."""
    agg_funcs = {c: ["mean", "std"] for c in value_cols}
    agg_funcs["Q"] = ["mean"]
    g = df.groupby(["geometry", "damp_idx"]).agg(agg_funcs)
    g.columns = ["_".join(c) for c in g.columns]
    g = g.reset_index().sort_values(["geometry", "damp_idx"])
    return g


# ════════════════════════════════════════════════════════════════════════
# 5. Plots
# ════════════════════════════════════════════════════════════════════════

def plot_metric_vs_Q(agg_df, value_col, ylabel, title, out_path, logy=False):
    fig, ax = plt.subplots(figsize=(7, 5))
    for geom in agg_df["geometry"].unique():
        sub = agg_df[agg_df["geometry"] == geom].sort_values("Q_mean")
        ax.errorbar(sub["Q_mean"], sub[f"{value_col}_mean"],
                     yerr=sub[f"{value_col}_std"],
                     marker=GEOM_MARKERS.get(geom, "o"),
                     color=GEOM_COLORS.get(geom, None),
                     label=geom, capsize=3, lw=1.5, ms=6)
    ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel(r"$Q = \omega_0 I / b$  (quality factor)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(title="Geometry")
    ax.axvline(1.0, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(1.0, ax.get_ylim()[1], "  Q=1\n  (crossover)", fontsize=8,
            color="gray", va="top")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_avalanche_histograms(summary_df, out_path, n_Q_examples=3):
    """For each geometry, plots P(s) in log-log for a few representative Q values
    (underdamped, intermediate, overdamped), with the MLE-fitted line overlaid."""
    geometries = sorted(summary_df["geometry"].unique())
    fig, axes = plt.subplots(1, len(geometries), figsize=(5.5 * len(geometries), 5),
                              sharey=True)
    if len(geometries) == 1:
        axes = [axes]

    for ax, geom in zip(axes, geometries):
        sub = summary_df[summary_df["geometry"] == geom]
        damp_idxs = sorted(sub["damp_idx"].unique())
        chosen = [damp_idxs[0], damp_idxs[len(damp_idxs) // 2], damp_idxs[-1]] \
            if len(damp_idxs) >= 3 else damp_idxs

        for di in chosen:
            rows = sub[sub["damp_idx"] == di]
            all_sizes = np.concatenate(rows["avalanche_sizes"].values) \
                if len(rows) else np.array([])
            all_sizes = all_sizes[all_sizes > 0]
            Q_label = rows["Q"].mean()
            if len(all_sizes) < 20:
                ax.plot([], [], " ", label=f"Q≈{Q_label:.2f}  (n={len(all_sizes)}, "
                                            f"too few for fit)")
                continue

            # log-binned histogram (for visualisation only -- the tau fit
            # used in the other plots comes from the MLE, not this histogram)
            bins = np.geomspace(all_sizes.min(), all_sizes.max(), 25)
            hist, edges = np.histogram(all_sizes, bins=bins, density=True)
            centers = np.sqrt(edges[:-1] * edges[1:])
            mask = hist > 0
            pts = ax.loglog(centers[mask], hist[mask], "o", ms=4,
                             label=f"Q≈{Q_label:.2f}  (n={len(all_sizes)})")
            color = pts[0].get_color()

            fit = fit_powerlaw_mle(all_sizes)
            if not np.isnan(fit["tau"]):
                xs_line = np.geomspace(fit["xmin"], all_sizes.max(), 50)
                ys_line = (fit["tau"] - 1.0) / fit["xmin"] * (xs_line / fit["xmin"]) ** (-fit["tau"])
                # anchor the MLE line at the histogram point closest to xmin
                # (rather than the last bin), where the fit actually starts --
                # avoids visual misalignment of the line when the tail has
                # few log-spaced bins.
                idx_anchor = np.argmin(np.abs(centers[mask] - fit["xmin"]))
                anchor_hist_y = hist[mask][idx_anchor]
                anchor_line_y = (fit["tau"] - 1.0) / fit["xmin"] * \
                    (centers[mask][idx_anchor] / fit["xmin"]) ** (-fit["tau"])
                scale = anchor_hist_y / (anchor_line_y + 1e-300)
                ax.loglog(xs_line, ys_line * scale, "--", lw=1.2,
                           color=color, alpha=0.7,
                           label=f"  MLE: τ={fit['tau']:.2f}±{fit['tau_err']:.2f}")

        ax.set_title(f"{geom}")
        ax.set_xlabel("avalanche size $s$")
        if ax is axes[0]:
            ax.set_ylabel(r"$P(s)$")
        ax.legend(fontsize=8)

    fig.suptitle("Avalanche size distribution (points: log-binned histogram; "
                  "dashed line: MLE tail fit)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_tau_vs_Q(tau_df, out_path):
    """Central plot of the study: tau (power-law exponent of the avalanche
    size distribution) as a function of Q, per geometry. The tau fit uses ALL
    avalanches aggregated across seeds for a given (geometry, damp_idx)
    -- see fit_powerlaw_per_geom_damp -- yielding one point (with MLE error bar)
    per geometry/Q, instead of a weak fit per individual run.

    Flat tau (tau ~ const) supports the 'incremental' result; systematic
    variation of tau with Q, or qualitative differences between geometries,
    is the signature of the 'surprising' result (universality-class breaking
    by anisotropic long-range dipolar coupling)."""
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for geom in tau_df["geometry"].unique():
        sub = tau_df[tau_df["geometry"] == geom].dropna(subset=["tau"]).sort_values("Q_mean")
        if len(sub) == 0:
            continue
        ax.errorbar(sub["Q_mean"], sub["tau"], yerr=sub["tau_err"],
                     marker=GEOM_MARKERS.get(geom, "o"),
                     color=GEOM_COLORS.get(geom, None),
                     label=f"{geom}  (mean n_tail={sub['n_tail'].mean():.0f})",
                     capsize=3, lw=1.5, ms=7)

    ax.set_xscale("log")
    ax.set_xlabel(r"$Q = \omega_0 I / b$  (quality factor)")
    ax.set_ylabel(r"power-law exponent $\tau$  [$P(s) \sim s^{-\tau}$]")
    ax.set_title("Central test: avalanche exponent vs. quality factor\n"
                  "(avalanches aggregated across seeds before MLE fit)")
    ax.axvline(1.0, color="gray", ls="--", lw=1, alpha=0.6)
    ax.legend(title="Geometry", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_domain_stats_vs_Q(relax_df, out_path):
    agg = aggregate_by_geom_damp(relax_df, ["mean_domain_size", "n_domains"])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for geom in agg["geometry"].unique():
        sub = agg[agg["geometry"] == geom].sort_values("Q_mean")
        axes[0].errorbar(sub["Q_mean"], sub["mean_domain_size_mean"],
                          yerr=sub["mean_domain_size_std"],
                          marker=GEOM_MARKERS.get(geom, "o"),
                          color=GEOM_COLORS.get(geom, None),
                          label=geom, capsize=3, lw=1.5, ms=6)
        axes[1].errorbar(sub["Q_mean"], sub["n_domains_mean"],
                          yerr=sub["n_domains_std"],
                          marker=GEOM_MARKERS.get(geom, "o"),
                          color=GEOM_COLORS.get(geom, None),
                          label=geom, capsize=3, lw=1.5, ms=6)

    for ax, ylabel, title in zip(
            axes,
            ["mean domain size (# needles)", "number of domains"],
            ["Mean domain size vs. Q", "Number of domains vs. Q"]):
        ax.set_xscale("log")
        ax.set_xlabel(r"$Q = \omega_0 I / b$")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.axvline(1.0, color="gray", ls="--", lw=1, alpha=0.6)
        ax.legend(title="Geometry")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_waiting_time_vs_Q(summary_df, out_path):
    """Plots mean inter-avalanche waiting time and the LRT result vs Q.
    Two panels: (left) mean IAT and CV; (right) LRT statistic.
    High LRT (>0) indicates bursty / scale-free timing (SOC-like)."""
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    agg_wt = aggregate_by_geom_damp(summary_df, ["wt_mean_iat", "wt_cv", "wt_lrt", "tau_corr"])

    for geom in agg_wt["geometry"].unique():
        sub = agg_wt[agg_wt["geometry"] == geom].sort_values("Q_mean")
        axes[0].errorbar(sub["Q_mean"], sub["wt_mean_iat_mean"],
                         yerr=sub["wt_mean_iat_std"],
                         marker=GEOM_MARKERS.get(geom, "o"),
                         color=GEOM_COLORS.get(geom, None),
                         label=geom, capsize=3, lw=1.5, ms=6)
        axes[1].plot(sub["Q_mean"], sub["wt_lrt_mean"],
                     marker=GEOM_MARKERS.get(geom, "o"),
                     color=GEOM_COLORS.get(geom, None),
                     label=geom, lw=1.5, ms=6)

    for ax, ylabel, title in zip(
            axes,
            ["mean inter-avalanche time (s)", "LRT: power-law vs log-normal (R)"],
            ["Mean waiting time vs. Q",
             "LRT result (R>0: power-law preferred)\n"
             "[Clauset, Shalizi & Newman 2009 §4]"]):
        ax.set_xscale("log")
        ax.set_xlabel(r"$Q = \omega_0 I / b$")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.axvline(1.0, color="gray", ls="--", lw=1, alpha=0.6)
        ax.legend(title="Geometry")
    axes[1].axhline(0.0, color="k", ls=":", lw=1)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_autocorrelation_vs_Q(summary_df, out_path):
    """Plots the M(t) correlation time tau_corr vs Q.
    Predicted by the echo-state/memory-capacity framework (s41598-025-93189-w)
    to peak near Q~1 (critical damping) and decrease in both limits."""
    fig, ax = plt.subplots(figsize=(7, 5))
    agg_ac = aggregate_by_geom_damp(summary_df, ["tau_corr"])

    for geom in agg_ac["geometry"].unique():
        sub = agg_ac[agg_ac["geometry"] == geom].sort_values("Q_mean")
        ax.errorbar(sub["Q_mean"], sub["tau_corr_mean"],
                    yerr=sub["tau_corr_std"],
                    marker=GEOM_MARKERS.get(geom, "o"),
                    color=GEOM_COLORS.get(geom, None),
                    label=geom, capsize=3, lw=1.5, ms=6)

    ax.set_xscale("log")
    ax.set_xlabel(r"$Q = \omega_0 I / b$")
    ax.set_ylabel(r"correlation time $\tau_{\mathrm{corr}}$ of $M(t)$ [s]")
    ax.set_title("Temporal autocorrelation of M(t) vs. Q\n"
                 "(prediction of the echo-state / ASI memory framework)")
    ax.axvline(1.0, color="gray", ls="--", lw=1, alpha=0.6)
    ax.legend(title="Geometry")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ════════════════════════════════════════════════════════════════════════
# 6. FORC analysis
# ════════════════════════════════════════════════════════════════════════

def load_forc_data(csv_path):
    """Reads a FORC CSV (columns: t, B, M, S, B_r, sweep_dir, curve_idx)
    and returns a DataFrame containing only the ascending-branch rows."""
    try:
        df = pd.read_csv(csv_path)
        if df.empty or 'B_r' not in df.columns:
            return pd.DataFrame()
        return df[df['sweep_dir'] == 'up'].copy()
    except Exception:
        return pd.DataFrame()


def compute_forc_distribution(forc_df, n_grid=100, smooth_sigma=2.0):
    """Computes the FORC distribution rho(B_r, B) = -1/2 * d^2M / (dB_r dB)
    from a DataFrame of ascending-branch FORC data.

    Parameters
    ----------
    forc_df      : DataFrame with columns B, M, B_r, curve_idx
    n_grid       : number of grid points along each axis
    smooth_sigma : Gaussian smoothing sigma (in grid cells) before differencing

    Returns
    -------
    rho          : 2D array  shape (n_Br, n_B)  — the FORC distribution
    Br_vals      : 1D array of reversal field values [T]
    B_vals       : 1D array of applied field values [T]
    """
    from scipy.ndimage import gaussian_filter
    from scipy.interpolate import griddata

    if forc_df.empty or len(forc_df['curve_idx'].unique()) < 3:
        return None, None, None

    B_all  = forc_df['B'].values.astype(float)
    Br_all = forc_df['B_r'].values.astype(float)
    M_all  = forc_df['M'].values.astype(float)

    B_min, B_max   = B_all.min(),  B_all.max()
    Br_min, Br_max = Br_all.min(), Br_all.max()

    B_vals  = np.linspace(B_min,  B_max,  n_grid)
    Br_vals = np.linspace(Br_min, Br_max, n_grid)
    BB, BrBr = np.meshgrid(B_vals, Br_vals)   # shape (n_Br, n_B)

    # Interpolate scattered (B, B_r, M) onto regular grid
    M_grid = griddata(
        np.column_stack([B_all, Br_all]), M_all,
        (BB, BrBr), method='linear', fill_value=np.nan)

    # Fill NaN borders with nearest-neighbour for smoother gradients
    mask_nan = np.isnan(M_grid)
    if mask_nan.any() and not mask_nan.all():
        M_nn = griddata(
            np.column_stack([B_all, Br_all]), M_all,
            (BB, BrBr), method='nearest')
        M_grid[mask_nan] = M_nn[mask_nan]

    # Gaussian smoothing (reduces noise before differentiation)
    if smooth_sigma > 0:
        M_grid = gaussian_filter(M_grid, sigma=smooth_sigma)

    # Mixed second derivative: rho = -1/2 * d/dB_r (dM/dB)
    dB  = B_vals[1]  - B_vals[0]   if n_grid > 1 else 1.0
    dBr = Br_vals[1] - Br_vals[0]  if n_grid > 1 else 1.0
    dMdB   = np.gradient(M_grid, dB,  axis=1)    # d/dB  along axis 1
    d2MdBrdB = np.gradient(dMdB,  dBr, axis=0)   # d/dBr along axis 0
    rho = -0.5 * d2MdBrdB

    # Zero out the unphysical region B < B_r
    for i, br in enumerate(Br_vals):
        rho[i, B_vals < br] = np.nan

    return rho, Br_vals, B_vals


def plot_forc_diagram(rho, Br_vals, B_vals, out_path, title=""):
    """Plots one FORC distribution diagram.

    Horizontal axis: B (applied field during ascending branch) in mT
    Vertical axis:   B_r (reversal field) in mT
    Colour scale:    symmetric diverging (blue-white-red) centred on zero.
    Dashed diagonal lines show constant coercivity B_c = (B - B_r)/2 and
    constant interaction field B_u = (B + B_r)/2.
    """
    if rho is None:
        return

    B_mT  = B_vals  * 1e3
    Br_mT = Br_vals * 1e3

    vmax = np.nanpercentile(np.abs(rho), 98)
    if vmax == 0:
        vmax = 1.0

    fig, ax = plt.subplots(figsize=(6, 5))
    pcm = ax.pcolormesh(B_mT, Br_mT, rho,
                         cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                         shading='auto')
    plt.colorbar(pcm, ax=ax, label=r"$\rho(B_r, B)$  [a.u.]")

    # Physical constraint line: B = B_r (boundary of FORC region)
    lim = min(B_mT.max(), Br_mT.max())
    lo  = max(B_mT.min(), Br_mT.min())
    ax.plot([lo, lim], [lo, lim], 'k--', lw=0.8, alpha=0.5, label=r"$B = B_r$")

    ax.set_xlabel(r"$B$  [mT]  (applied field)")
    ax.set_ylabel(r"$B_r$  [mT]  (reversal field)")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_forc_vs_Q(forc_rows, args_out_dir):
    """Panel plot: one FORC diagram per Q value (selected dampings),
    one row per geometry.  Output: forc_diagrams.png"""
    geometries = sorted(forc_rows['geometry'].unique())
    n_geom = len(geometries)
    if n_geom == 0:
        return

    # Pick up to 4 Q values across the full range for display
    all_Q = sorted(forc_rows['Q'].unique())
    if len(all_Q) > 4:
        idxs = np.round(np.linspace(0, len(all_Q) - 1, 4)).astype(int)
        Q_display = [all_Q[i] for i in idxs]
    else:
        Q_display = all_Q
    n_Q = len(Q_display)

    fig, axes = plt.subplots(n_geom, n_Q,
                              figsize=(4.5 * n_Q, 4.0 * n_geom),
                              squeeze=False)

    for row_i, geom in enumerate(geometries):
        sub_geom = forc_rows[forc_rows['geometry'] == geom]
        for col_i, Q_target in enumerate(Q_display):
            ax = axes[row_i, col_i]
            # pick the damp_idx whose Q is closest to Q_target
            q_arr = sub_geom['Q'].values
            best_idx = np.argmin(np.abs(q_arr - Q_target))
            best_damp_idx = sub_geom.iloc[best_idx]['damp_idx']
            runs = sub_geom[sub_geom['damp_idx'] == best_damp_idx]
            Q_actual = runs['Q'].mean()

            # aggregate FORC data across seeds
            all_dfs = []
            for _, run_row in runs.iterrows():
                df = load_forc_data(run_row['csv_path'])
                if not df.empty:
                    all_dfs.append(df)
            if not all_dfs:
                ax.text(0.5, 0.5, 'no data', ha='center', va='center',
                        transform=ax.transAxes, fontsize=9)
                ax.set_title(f"{geom}  Q={Q_actual:.2f}", fontsize=9)
                continue

            combined = pd.concat(all_dfs, ignore_index=True)
            rho, Br_vals, B_vals = compute_forc_distribution(combined)

            if rho is None:
                ax.text(0.5, 0.5, 'too few curves', ha='center', va='center',
                        transform=ax.transAxes, fontsize=9)
                ax.set_title(f"{geom}  Q={Q_actual:.2f}", fontsize=9)
                continue

            B_mT  = B_vals  * 1e3
            Br_mT = Br_vals * 1e3
            vmax = np.nanpercentile(np.abs(rho), 98)
            if vmax == 0:
                vmax = 1.0
            pcm = ax.pcolormesh(B_mT, Br_mT, rho,
                                 cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                                 shading='auto')
            plt.colorbar(pcm, ax=ax)
            lo  = max(B_mT.min(),  Br_mT.min())
            lim = min(B_mT.max(), Br_mT.max())
            ax.plot([lo, lim], [lo, lim], 'k--', lw=0.7, alpha=0.5)
            ax.set_title(f"{geom}   Q = {Q_actual:.2f}", fontsize=9)
            if col_i == 0:
                ax.set_ylabel(r"$B_r$  [mT]", fontsize=8)
            if row_i == n_geom - 1:
                ax.set_xlabel(r"$B$  [mT]", fontsize=8)

    fig.suptitle("FORC distributions  —  ascending branch  (ρ = −½ ∂²M/∂Bₐ∂B)",
                 fontsize=11)
    fig.tight_layout()
    out_path = os.path.join(args_out_dir, "forc_diagrams.png")
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_path}")


# ════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep_dir", type=str, required=True,
                     help="Directory written by --out_dir in damping_sweep.py")
    ap.add_argument("--out_dir", type=str, required=True)
    ap.add_argument("--mad_threshold", type=float, default=4.0,
                     help="MAD multiplier (threshold in units of MAD) for dM/dB avalanche detection")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    manifest_path = os.path.join(args.sweep_dir, "manifest.csv")
    if not os.path.exists(manifest_path):
        print(f"ERROR: manifest.csv not found in {args.sweep_dir}")
        sys.exit(1)
    manifest = pd.read_csv(manifest_path)

    hyst_rows  = manifest[manifest["stage"] == "hysteresis"]
    relax_rows = manifest[manifest["stage"] == "relax"]
    forc_rows  = manifest[manifest["stage"] == "forc"]

    print(f"Manifest: {len(hyst_rows)} hysteresis runs, "
          f"{len(relax_rows)} relaxation runs, "
          f"{len(forc_rows)} FORC runs.")

    # ── Stage 1: process each hysteresis run ────────────────────────────
    print("\nProcessing hysteresis runs (area, Bc, Mr, avalanches, MLE)...")
    hyst_results = []
    n_hyst = len(hyst_rows)
    _bar(0, n_hyst, label="HYSTERESIS")
    for i, (_, row) in enumerate(hyst_rows.iterrows()):
        hyst_results.append(process_hysteresis_run(row, args.mad_threshold))
        _bar(i + 1, n_hyst, label="HYSTERESIS", extra=row['tag'])

    summary_df = pd.DataFrame(hyst_results)
    # CSV plano (sem a coluna de arrays de avalanche, que nao serializa bem em CSV)
    summary_csv = summary_df.drop(columns=["avalanche_sizes", "avalanche_sizes_S",
                                           "inter_arrival_times"])
    summary_csv_path = os.path.join(args.out_dir, "summary_hysteresis.csv")
    summary_csv.to_csv(summary_csv_path, index=False)
    print(f"Saved: {summary_csv_path}")

    # ── ajuste de lei de potencia (dM/dB method): agrega avalanches entre seeds
    print("\nFitting power law [dM/dB method] (avalanches aggregated across seeds)...")
    tau_df = fit_powerlaw_per_geom_damp(summary_df, size_col="avalanche_sizes")
    tau_csv_path = os.path.join(args.out_dir, "powerlaw_fits.csv")
    tau_df.to_csv(tau_csv_path, index=False)
    print(f"Saved: {tau_csv_path}")
    print(tau_df.to_string(index=False))

    # ── ajuste de lei de potencia (S-column method)
    print("\nFitting power law [S-column method] (avalanches aggregated across seeds)...")
    tau_df_S = fit_powerlaw_per_geom_damp(summary_df, size_col="avalanche_sizes_S")
    tau_S_csv_path = os.path.join(args.out_dir, "powerlaw_fits_S.csv")
    tau_df_S.to_csv(tau_S_csv_path, index=False)
    print(f"Saved: {tau_S_csv_path}")

    # ── Stage 2: process each relaxation run ────────────────────────────
    print("\nProcessing relaxation runs (domain statistics)...")
    relax_results = []
    n_relax = len(relax_rows)
    _bar(0, n_relax, label="RELAX")
    for i, (_, row) in enumerate(relax_rows.iterrows()):
        relax_results.append(process_relax_run(row))
        _bar(i + 1, n_relax, label="RELAX", extra=row['tag'])
    relax_df = pd.DataFrame(relax_results)
    relax_csv_path = os.path.join(args.out_dir, "summary_relax.csv")
    relax_df.to_csv(relax_csv_path, index=False)
    print(f"Saved: {relax_csv_path}")

    # ── aggregate (mean+std over seeds) of loop metrics ──────────────────
    print("\nGenerating plots...")
    agg_loop = aggregate_by_geom_damp(summary_df, ["area", "Bc", "Mr"])
    agg_loop.to_csv(os.path.join(args.out_dir, "aggregated_loop_metrics.csv"), index=False)

    plot_metric_vs_Q(agg_loop, "area", "hysteresis loop area  [T]",
                      "Loop area vs. quality factor",
                      os.path.join(args.out_dir, "area_vs_Q.png"))
    plot_metric_vs_Q(agg_loop, "Bc", r"coercive field $B_c$  [T]",
                      "Coercive field vs. quality factor",
                      os.path.join(args.out_dir, "Bc_vs_Q.png"))
    plot_metric_vs_Q(agg_loop, "Mr", r"remanence $M_r$",
                      "Remanence vs. quality factor",
                      os.path.join(args.out_dir, "Mr_vs_Q.png"))

    # ── central plot: tau(Q) ─────────────────────────────────────────
    plot_tau_vs_Q(tau_df, os.path.join(args.out_dir, "tau_vs_Q.png"))
    plot_tau_vs_Q(tau_df_S, os.path.join(args.out_dir, "tau_vs_Q_Scol.png"))

    # ── avalanche histograms (visual, selected geometries/Q) ────────────
    plot_avalanche_histograms(summary_df, os.path.join(args.out_dir, "avalanche_histograms.png"))

    # ── domain statistics vs Q ──────────────────────────────────────
    plot_domain_stats_vs_Q(relax_df, os.path.join(args.out_dir, "domain_size_vs_Q.png"))

    # ── inter-avalanche waiting times (S-column) ────────────────────
    plot_waiting_time_vs_Q(summary_df, os.path.join(args.out_dir, "waiting_time_vs_Q.png"))

    # ── temporal autocorrelation of M(t) ───────────────────────────
    plot_autocorrelation_vs_Q(summary_df, os.path.join(args.out_dir, "autocorr_vs_Q.png"))

    # ── Stage 3: FORC analysis ────────────────────────────────────
    if len(forc_rows) > 0:
        print(f"\nProcessing FORC runs ({len(forc_rows)} curves)...")
        _bar(0, len(forc_rows), label="FORC")
        for i, (_, row) in enumerate(forc_rows.iterrows(), 1):
            _bar(i, len(forc_rows), label="FORC",
                 extra=f"{row['tag']}  Q={row['Q']:.3f}")
        print()  # finalise bar
        plot_forc_vs_Q(forc_rows, args.out_dir)
    else:
        print("\nNo FORC runs in manifest — skipping FORC analysis.")

    print(f"\nDone. Results saved to: {args.out_dir}")
    print("Files written:")
    for fname in ["summary_hysteresis.csv", "summary_relax.csv",
                  "aggregated_loop_metrics.csv",
                  "powerlaw_fits.csv", "powerlaw_fits_S.csv",
                  "area_vs_Q.png", "Bc_vs_Q.png",
                  "Mr_vs_Q.png", "tau_vs_Q.png", "tau_vs_Q_Scol.png",
                  "avalanche_histograms.png",
                  "domain_size_vs_Q.png",
                  "waiting_time_vs_Q.png",
                  "autocorr_vs_Q.png",
                  "forc_diagrams.png"]:
        if os.path.exists(os.path.join(args.out_dir, fname)):
            print(f"  - {fname}")



if __name__ == "__main__":
    main()
