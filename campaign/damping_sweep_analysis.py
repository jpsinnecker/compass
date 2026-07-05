#!/usr/bin/env python3
"""
damping_sweep_analysis.py
==========================

Parte 2 do estudo de damping/histerese: le os resultados produzidos por
damping_sweep.py (manifest.csv + CSVs por rodada) e produz:

  1. Para cada rodada de histerese: area do laco, campo coercivo Bc,
     remanencia Mr.
  2. Estatistica de avalanches: deteccao de eventos em dM/dB dentro de
     cada segmento de rampa, e ajuste de lei de potencia P(s) ~ s^-tau
     via maxima verossimilhanca (metodo Clauset-Shalizi-Newman), com
     selecao de x_min por minimizacao da distancia KS.
  3. Para cada rodada de relaxacao livre (Etapa 2): estatistica de
     tamanho de dominios (label_magnetic_domains, ja calculada no
     momento da simulacao e salva no meta JSON).
  4. Graficos finais:
       - area_vs_Q.png, Bc_vs_Q.png, Mr_vs_Q.png   (uma curva por geometria)
       - tau_vs_Q.png                               (o teste central)
       - avalanche_histograms.png                   (P(s) por geometria, alguns Q)
       - domain_size_vs_Q.png                       (tamanho medio/expoente de dominio vs Q)
  5. Um summary.csv consolidado (uma linha por rodada de histerese, com
     todas as metricas), para inspecao/tabulacao posterior.

Uso
---
    python3 damping_sweep_analysis.py --sweep_dir /home/jps/sweep_results \
        --out_dir /home/jps/sweep_analysis

Dependencias: numpy, pandas, matplotlib, scipy.
Opcional: powerlaw (se ausente, usa um MLE proprio equivalente ao metodo
Clauset-Shalizi-Newman para o caso continuo, com busca de x_min por KS).
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
# 1. Metricas de laco de histerese: area, coercividade, remanencia
# ════════════════════════════════════════════════════════════════════════

def loop_area(B, M):
    """Area do laco de histerese fechado, via integral de linha
    (formula do shoelace aplicada a curva fechada no plano B-M).
    Equivalente a integral de M dB ao longo do ciclo completo."""
    B = np.asarray(B)
    M = np.asarray(M)
    # garante que o ciclo esta fechado (ultimo ponto ~ primeiro ponto);
    # se nao estiver exatamente fechado (ruido numerico no ultimo passo),
    # fecha explicitamente para a integral de linha nao vazar area espuria.
    if not np.isclose(B[0], B[-1], atol=1e-9) or not np.isclose(M[0], M[-1], atol=1e-6):
        B = np.append(B, B[0])
        M = np.append(M, M[0])
    # area de poligono fechado no plano (B, M):
    area = 0.5 * np.abs(np.sum(B[:-1] * M[1:] - B[1:] * M[:-1]))
    return area


def find_zero_crossings(x, y):
    """Retorna os valores de x onde y cruza zero, via interpolacao linear
    entre pontos consecutivos de sinais opostos."""
    sign = np.sign(y)
    sign[sign == 0] = 1  # evita falsos cruzamentos em y exatamente 0
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
    Bc (campo coercivo)  : valor de |B| onde M cruza zero (media dos
                            cruzamentos encontrados, ja que o laco tem
                            mais de um cruzamento ao longo do ciclo).
    Mr (remanencia)       : valor de |M| onde B cruza zero.
    """
    B = np.asarray(B)
    M = np.asarray(M)

    # Bc: onde M(t) cruza zero -> le o B correspondente
    M_crossings_idx_B = find_zero_crossings(B, M)  # x=B, y=M -> retorna B onde M=0... espera, precisamos inverter
    # find_zero_crossings(x,y) acha x onde y=0. Para Bc queremos o valor de
    # B onde M=0 -> x deve ser B, y deve ser M. Jah esta correto acima.
    Bc_candidates = np.abs(M_crossings_idx_B)
    Bc = np.median(Bc_candidates) if len(Bc_candidates) > 0 else np.nan

    # Mr: onde B(t) cruza zero -> le o M correspondente. Aqui x=B nao serve
    # como eixo de interpolacao (B passa por zero varias vezes); melhor
    # interpolar M em funcao do indice/tempo nos pontos onde B=0.
    B_crossings_idx = find_zero_crossings(np.arange(len(B)), B)
    if len(B_crossings_idx) > 0:
        M_at_Bzero = np.interp(B_crossings_idx, np.arange(len(M)), M)
        Mr = np.median(np.abs(M_at_Bzero))
    else:
        Mr = np.nan

    return Bc, Mr


# ════════════════════════════════════════════════════════════════════════
# 2. Deteccao de avalanches via dM/dB dentro de cada segmento de rampa
# ════════════════════════════════════════════════════════════════════════

def split_hysteresis_segments(B):
    """
    O ciclo de histerese tem 5 segmentos lineares (0->+Bmax->0->-Bmax->0->+Bmax).
    Detecta os pontos de inflexao (onde dB/dt muda de sinal) para dividir
    a serie em segmentos de rampa monotona -- avalanches sao detectadas
    DENTRO de cada segmento, nunca atravessando uma inflexao (onde dB/dt
    passa por zero e dM/dB diverge numericamente sem ser avalanche real).
    """
    dB = np.diff(B)
    sign_dB = np.sign(dB)
    sign_dB[sign_dB == 0] = 0  # pontos de B estacionario (raro, equilibrio numerico)
    change_idx = np.where(np.diff(sign_dB) != 0)[0] + 1
    bounds = [0] + list(change_idx + 1) + [len(B)]
    bounds = sorted(set(b for b in bounds if 0 <= b <= len(B)))
    segments = [(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1)
                if bounds[i + 1] - bounds[i] >= 5]  # ignora segmentos degenerados (<5 pontos)
    return segments


def detect_avalanches(t, B, M, mad_threshold=4.0, min_event_len=1):
    """
    Para cada segmento de rampa monotona em B, calcula dM/dB e marca como
    "avalanche" os pontos onde |dM/dB - mediana_local| excede
    mad_threshold * MAD (median absolute deviation, robusto a outliers --
    apropriado aqui porque a propria distribuicao de avalanches e de
    cauda pesada, o que inflaria um desvio-padrao comum).

    Eventos contiguos (passos consecutivos acima do limiar) sao agrupados
    em uma unica avalanche; o "tamanho" da avalanche e a soma de |dM|
    absorvida durante o evento (proxy para a magnetizacao revertida
    naquela avalanche).

    Retorna: array de tamanhos de avalanche (um valor por evento detectado,
    agregando todos os segmentos do ciclo).
    """
    segments = split_hysteresis_segments(B)
    sizes = []

    for (i0, i1) in segments:
        Bs = B[i0:i1]
        Ms = M[i0:i1]
        dB_seg = np.diff(Bs)
        dM_seg = np.diff(Ms)
        # evita divisao por dB ~ 0 (nao deveria ocorrer dentro de um
        # segmento monotono valido, mas protege contra ruido numerico)
        valid = np.abs(dB_seg) > 1e-14
        if valid.sum() < 5:
            continue
        dMdB = np.full_like(dB_seg, np.nan)
        dMdB[valid] = dM_seg[valid] / dB_seg[valid]

        local = dMdB[valid]
        med = np.median(local)
        mad = np.median(np.abs(local - med)) + 1e-300  # evita MAD=0 exato
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
            # tamanho da avalanche = soma de |dM| absorvida no evento
            size = np.sum(np.abs(dM_seg[g]))
            if size > 0:
                sizes.append(size)

    return np.array(sizes)


# ════════════════════════════════════════════════════════════════════════
# 3. Ajuste de lei de potencia: MLE (metodo Clauset-Shalizi-Newman)
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
    Ajusta lei de potencia aos tamanhos de avalanche via MLE, escolhendo
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
# 4. Pipeline principal: le manifest, processa cada rodada, agrega
# ════════════════════════════════════════════════════════════════════════

def process_hysteresis_run(row, mad_threshold):
    df = pd.read_csv(row["csv_path"])
    B, M = df["B"].values, df["M"].values
    t = df["t"].values

    area = loop_area(B, M)
    Bc, Mr = coercivity_and_remanence(B, M)
    avalanche_sizes = detect_avalanches(t, B, M, mad_threshold=mad_threshold)

    with open(row["meta_path"]) as f:
        meta = json.load(f)

    return dict(
        tag=row["tag"], geometry=row["geometry"], damp_idx=row["damp_idx"],
        seed=row["seed"], damping=row["damping"], Q=row["Q"],
        area=area, Bc=Bc, Mr=Mr,
        n_avalanches=len(avalanche_sizes),
        avalanche_sizes=avalanche_sizes,
        B_max=meta.get("B_max", np.nan),
    )


def fit_powerlaw_per_geom_damp(summary_df):
    """
    Agrega os tamanhos de avalanche de TODAS as seeds de cada combinacao
    (geometria, damp_idx) antes de ajustar a lei de potencia -- um unico
    ajuste MLE por ponto (geometria, Q), em vez de um ajuste fraco por
    rodada individual (que tipicamente tem poucas avalanches e gera tau
    instavel / NaN por falta de amostra).

    Retorna um DataFrame com uma linha por (geometria, damp_idx): Q medio,
    tau, tau_err, xmin, ks, n_tail, n_total (agregado).
    """
    rows = []
    for (geom, di), grp in summary_df.groupby(["geometry", "damp_idx"]):
        all_sizes = np.concatenate(grp["avalanche_sizes"].values) \
            if len(grp) else np.array([])
        fit = fit_powerlaw_mle(all_sizes)
        rows.append(dict(
            geometry=geom, damp_idx=di, Q_mean=grp["Q"].mean(),
            tau=fit["tau"], tau_err=fit["tau_err"], xmin=fit["xmin"],
            ks=fit["ks"], n_tail=fit["n_tail"], n_total=fit.get("n_total", 0),
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
    sao esperadas e o Q de cada rodada individual ja esta correto por linha;
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
    ax.set_xlabel(r"$Q = \omega_0 I / b$  (fator de qualidade)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(title="Geometria")
    ax.axvline(1.0, color="gray", ls="--", lw=1, alpha=0.6)
    ax.text(1.0, ax.get_ylim()[1], "  Q=1\n  (crossover)", fontsize=8,
            color="gray", va="top")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_avalanche_histograms(summary_df, out_path, n_Q_examples=3):
    """Para cada geometria, plota P(s) em log-log para alguns valores de Q
    representativos (subamortecido, intermediario, sobreamortecido), com
    a reta MLE ajustada sobreposta."""
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
                                            f"insuficiente p/ ajuste)")
                continue

            # histograma log-binado (so para visualizacao -- o ajuste tau
            # usado nos demais graficos vem do MLE, nao deste histograma)
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
                # ancora a linha MLE no ponto do histograma mais proximo de
                # xmin (em vez do ultimo bin), que e onde o ajuste de fato
                # comeca a valer -- evita desalinhamento visual da reta
                # quando a cauda tem poucos bins log-espacados.
                idx_anchor = np.argmin(np.abs(centers[mask] - fit["xmin"]))
                anchor_hist_y = hist[mask][idx_anchor]
                anchor_line_y = (fit["tau"] - 1.0) / fit["xmin"] * \
                    (centers[mask][idx_anchor] / fit["xmin"]) ** (-fit["tau"])
                scale = anchor_hist_y / (anchor_line_y + 1e-300)
                ax.loglog(xs_line, ys_line * scale, "--", lw=1.2,
                           color=color, alpha=0.7,
                           label=f"  MLE: τ={fit['tau']:.2f}±{fit['tau_err']:.2f}")

        ax.set_title(f"{geom}")
        ax.set_xlabel("tamanho da avalanche $s$")
        if ax is axes[0]:
            ax.set_ylabel(r"$P(s)$")
        ax.legend(fontsize=8)

    fig.suptitle("Distribuição de tamanhos de avalanche (pontos: histograma log-binado; "
                  "linha tracejada: ajuste MLE da cauda)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_tau_vs_Q(tau_df, out_path):
    """O grafico central do estudo: tau (expoente da lei de potencia de
    avalanches) em funcao de Q, por geometria. O ajuste de tau usa TODAS
    as avalanches agregadas entre seeds de uma mesma (geometria, damp_idx)
    -- ver fit_powerlaw_per_geom_damp -- o que da um unico ponto (com
    barra de erro do MLE) por geometria/Q, em vez de um ajuste fraco por
    rodada individual.

    Achatamento (tau ~ const) favorece o resultado 'incremental'; variacao
    sistematica de tau com Q, ou diferenca qualitativa entre geometrias, e
    a assinatura do resultado 'surpreendente' (quebra de universalidade
    por acoplamento dipolar anisotropico de longo alcance)."""
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    for geom in tau_df["geometry"].unique():
        sub = tau_df[tau_df["geometry"] == geom].dropna(subset=["tau"]).sort_values("Q_mean")
        if len(sub) == 0:
            continue
        ax.errorbar(sub["Q_mean"], sub["tau"], yerr=sub["tau_err"],
                     marker=GEOM_MARKERS.get(geom, "o"),
                     color=GEOM_COLORS.get(geom, None),
                     label=f"{geom}  (n_tail médio={sub['n_tail'].mean():.0f})",
                     capsize=3, lw=1.5, ms=7)

    ax.set_xscale("log")
    ax.set_xlabel(r"$Q = \omega_0 I / b$  (fator de qualidade)")
    ax.set_ylabel(r"expoente da lei de potência $\tau$  [$P(s) \sim s^{-\tau}$]")
    ax.set_title("Teste central: expoente de avalanche vs. fator de qualidade\n"
                  "(avalanches agregadas entre seeds antes do ajuste MLE)")
    ax.axvline(1.0, color="gray", ls="--", lw=1, alpha=0.6)
    ax.legend(title="Geometria", fontsize=9)
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
            ["tamanho médio de domínio (nº agulhas)", "número de domínios"],
            ["Tamanho médio de domínio vs. Q", "Número de domínios vs. Q"]):
        ax.set_xscale("log")
        ax.set_xlabel(r"$Q = \omega_0 I / b$")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.axvline(1.0, color="gray", ls="--", lw=1, alpha=0.6)
        ax.legend(title="Geometria")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep_dir", type=str, required=True,
                     help="Diretorio passado como --out_dir ao damping_sweep.py")
    ap.add_argument("--out_dir", type=str, required=True)
    ap.add_argument("--mad_threshold", type=float, default=4.0,
                     help="Limiar (em multiplos de MAD) para deteccao de avalanche em dM/dB")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    manifest_path = os.path.join(args.sweep_dir, "manifest.csv")
    if not os.path.exists(manifest_path):
        print(f"ERRO: manifest.csv nao encontrado em {args.sweep_dir}")
        sys.exit(1)
    manifest = pd.read_csv(manifest_path)

    hyst_rows = manifest[manifest["stage"] == "hysteresis"]
    relax_rows = manifest[manifest["stage"] == "relax"]

    print(f"Manifest: {len(hyst_rows)} rodadas de histerese, "
          f"{len(relax_rows)} rodadas de relaxacao.")

    # ── Etapa 1: processa cada rodada de histerese ──────────────────────
    print("\nProcessando rodadas de histerese (area, Bc, Mr, avalanches, MLE)...")
    hyst_results = []
    for i, (_, row) in enumerate(hyst_rows.iterrows()):
        print(f"  [{i+1}/{len(hyst_rows)}] {row['tag']}", end="\r")
        hyst_results.append(process_hysteresis_run(row, args.mad_threshold))
    print()

    summary_df = pd.DataFrame(hyst_results)
    # CSV plano (sem a coluna de arrays de avalanche, que nao serializa bem em CSV)
    summary_csv = summary_df.drop(columns=["avalanche_sizes"])
    summary_csv_path = os.path.join(args.out_dir, "summary_hysteresis.csv")
    summary_csv.to_csv(summary_csv_path, index=False)
    print(f"Salvo: {summary_csv_path}")

    # ── ajuste de lei de potencia: agrega avalanches entre seeds ANTES do
    #    MLE, para cada (geometria, damp_idx) -- ver docstring da funcao.
    print("\nAjustando lei de potencia (avalanches agregadas entre seeds)...")
    tau_df = fit_powerlaw_per_geom_damp(summary_df)
    tau_csv_path = os.path.join(args.out_dir, "powerlaw_fits.csv")
    tau_df.to_csv(tau_csv_path, index=False)
    print(f"Salvo: {tau_csv_path}")
    print(tau_df.to_string(index=False))

    # ── Etapa 2: processa cada rodada de relaxacao ───────────────────────
    print("\nProcessando rodadas de relaxacao (estatistica de dominios)...")
    relax_results = [process_relax_run(row) for _, row in relax_rows.iterrows()]
    relax_df = pd.DataFrame(relax_results)
    relax_csv_path = os.path.join(args.out_dir, "summary_relax.csv")
    relax_df.to_csv(relax_csv_path, index=False)
    print(f"Salvo: {relax_csv_path}")

    # ── agregados (media+std sobre seeds) das metricas de laco ──────────
    print("\nGerando graficos...")
    agg_loop = aggregate_by_geom_damp(summary_df, ["area", "Bc", "Mr"])
    agg_loop.to_csv(os.path.join(args.out_dir, "aggregated_loop_metrics.csv"), index=False)

    plot_metric_vs_Q(agg_loop, "area", "área do laço de histerese  [T]",
                      "Área do laço vs. fator de qualidade",
                      os.path.join(args.out_dir, "area_vs_Q.png"))
    plot_metric_vs_Q(agg_loop, "Bc", r"campo coercivo $B_c$  [T]",
                      "Campo coercivo vs. fator de qualidade",
                      os.path.join(args.out_dir, "Bc_vs_Q.png"))
    plot_metric_vs_Q(agg_loop, "Mr", r"remanência $M_r$",
                      "Remanência vs. fator de qualidade",
                      os.path.join(args.out_dir, "Mr_vs_Q.png"))

    # ── o grafico central: tau(Q) ─────────────────────────────────────
    plot_tau_vs_Q(tau_df, os.path.join(args.out_dir, "tau_vs_Q.png"))

    # ── histogramas de avalanche (visual, algumas geometrias/Q) ─────────
    plot_avalanche_histograms(summary_df, os.path.join(args.out_dir, "avalanche_histograms.png"))

    # ── dominios vs Q ────────────────────────────────────────────────
    plot_domain_stats_vs_Q(relax_df, os.path.join(args.out_dir, "domain_size_vs_Q.png"))

    print(f"\nConcluido. Resultados em: {args.out_dir}")
    print("Arquivos gerados:")
    for fname in ["summary_hysteresis.csv", "summary_relax.csv",
                  "aggregated_loop_metrics.csv", "powerlaw_fits.csv",
                  "area_vs_Q.png", "Bc_vs_Q.png",
                  "Mr_vs_Q.png", "tau_vs_Q.png", "avalanche_histograms.png",
                  "domain_size_vs_Q.png"]:
        print(f"  - {fname}")


if __name__ == "__main__":
    main()
