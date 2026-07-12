#!/usr/bin/env python3
"""
avalanche_processor.py - offline spatio-temporal avalanche clustering and
MLE model comparison for compassV2.1 event logs.
Version: 1.0.0, July 2026

Pipeline
--------
1. Scan a results directory produced by compassV2_1.py for run triplets:
       meta/{tag}.json, data/{tag}_events.csv, states/{tag}_initial.npz
2. For each run and each activity channel SEPARATELY ("field", "angle"),
   cluster committed events into avalanches with the causal link rule:

       event j joins the cluster of event i  <=>
           site(j) is a spatial neighbour of site(i), r_ij <= r_link, AND
           0 <= t_j - t_i <= t_link.

   Links are transitive (union-find with path compression), so a cascade
   propagating needle-to-needle is one avalanche of any temporal length,
   while simultaneous but spatially disconnected activity stays separate.
   This is a PAIRWISE rule, not a fixed time bin: fixed bins conflate
   independent concurrent cascades and were deliberately rejected.

3. Aggregate avalanche sizes ACROSS SEEDS for each physical condition
   (default grouping: geometry, field_mode, Q rounded to 3 significant
   figures) before any statistical fit. Single-run fits are not produced.

4. For each (condition, channel, t_link) tail: Clauset-style discrete
   power-law MLE with KS-optimal xmin, compared against tail-conditioned
   log-normal and exponential alternatives via Vuong likelihood-ratio
   tests. Tails with fewer than --min_tail events are labelled
   "insufficient" and no exponent is reported for them.

5. t_link is a SWEPT parameter (default 1, 2, 4 T0). The causal
   neighbour response time is ~T0 only at high Q; overdamped arrays
   respond more slowly, so any single t_link imposes a Q-correlated bias.
   Only conclusions stable across the t_link sweep should be reported.

Outputs (in --out_dir)
----------------------
  avalanches_{condition}_{channel}_tlink{X}.npz : sizes, durations,
        participants, waiting times, per-run provenance.
  summary_avalanches.csv : one row per (condition, channel, t_link) with
        counts, s_max, xmin, alpha, KS distance, tail size, log-normal
        (mu, sigma), exponential lambda, Vuong LR statistics and p-values,
        and a conservative verdict string.
  optional CCDF plots with fitted models (--make_plot).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    from scipy.optimize import minimize
    from scipy.special import erfc as _erfc
    _SCIPY = True
except Exception:  # pragma: no cover
    _SCIPY = False

    def _erfc(x):
        return math.erfc(x)


# =============================================================================
# Union-find
# =============================================================================


class UnionFind:
    __slots__ = ("parent", "rank")

    def __init__(self, n: int):
        self.parent = np.arange(n, dtype=np.int64)
        self.rank = np.zeros(n, dtype=np.int8)

    def find(self, a: int) -> int:
        p = self.parent
        while p[a] != a:
            p[a] = p[p[a]]
            a = p[a]
        return a

    def union(self, a: int, b: int):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


# =============================================================================
# Run discovery and loading
# =============================================================================


def discover_runs(runs_dir: Path) -> List[Dict]:
    """Find all runs that have metadata, an event log, and an initial state."""
    runs = []
    for meta_path in sorted(runs_dir.glob("**/meta/*.json")):
        base = meta_path.parent.parent
        tag = meta_path.stem
        ev = base / "data" / f"{tag}_events.csv"
        st = base / "states" / f"{tag}_initial.npz"
        if not ev.exists():
            print(f"[skip] {tag}: no event log (run without --event_log?)")
            continue
        if not st.exists():
            print(f"[skip] {tag}: no initial-state NPZ (positions required)")
            continue
        meta = json.loads(meta_path.read_text())
        runs.append({"tag": tag, "meta": meta, "events": ev, "state": st})
    return runs


def condition_key(meta: Dict, group_by: Sequence[str]) -> Tuple:
    cfg = meta.get("config", {})
    der = meta.get("derived", {})
    key = []
    for g in group_by:
        if g == "Q":
            q = der.get("Q", float("nan"))
            key.append(float(f"{q:.3g}") if math.isfinite(q) else "inf")
        elif g in cfg:
            key.append(cfg[g])
        elif g in der:
            key.append(der[g])
        else:
            key.append("NA")
    return tuple(key)


def neighbour_lists(xs: np.ndarray, ys: np.ndarray, r_link: float) -> List[np.ndarray]:
    """Adjacency lists via cell binning; O(K) memory, near O(K) time."""
    K = xs.size
    cell = r_link
    ix = np.floor(xs / cell).astype(np.int64)
    iy = np.floor(ys / cell).astype(np.int64)
    grid: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for i in range(K):
        grid[(int(ix[i]), int(iy[i]))].append(i)
    r2 = r_link * r_link
    neigh: List[np.ndarray] = []
    for i in range(K):
        cand: List[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                cand.extend(grid.get((int(ix[i]) + dx, int(iy[i]) + dy), ()))
        cand = [j for j in cand if j != i
                and (xs[i] - xs[j]) ** 2 + (ys[i] - ys[j]) ** 2 <= r2]
        neigh.append(np.asarray(cand, dtype=np.int64))
    return neigh


# =============================================================================
# Clustering
# =============================================================================


def cluster_events(
    t: np.ndarray,
    nid: np.ndarray,
    neigh: List[np.ndarray],
    t_link: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Cluster time-sorted events with the pairwise causal link rule.

    Maintains, per needle, the index of its most recent event; each incoming
    event links to the most recent event on its own site and on neighbouring
    sites if that event lies within t_link in the past. Transitivity through
    union-find recovers arbitrarily long causal chains, so linking only to
    the MOST RECENT event per site loses nothing: any older in-window event
    on that site is already unioned with the newer one through the site's
    own event chain. Complexity O(E * z * alpha(E)).

    Returns (sizes, durations, participants, waiting_times).
    """
    E = t.size
    if E == 0:
        z = np.zeros(0)
        return z, z, z, z
    order = np.argsort(t, kind="stable")
    t = t[order]
    nid = nid[order]
    uf = UnionFind(E)
    last_ev = {}  # needle_id -> latest event index
    for e in range(E):
        i = int(nid[e])
        te = t[e]
        prev = last_ev.get(i)
        if prev is not None and te - t[prev] <= t_link:
            uf.union(e, prev)
        for j in neigh[i]:
            pj = last_ev.get(int(j))
            if pj is not None and te - t[pj] <= t_link:
                uf.union(e, pj)
        last_ev[i] = e

    roots = np.fromiter((uf.find(e) for e in range(E)), dtype=np.int64, count=E)
    clusters: Dict[int, List[int]] = defaultdict(list)
    for e, r in enumerate(roots):
        clusters[int(r)].append(e)

    sizes, durs, parts, t_start = [], [], [], []
    for members in clusters.values():
        m = np.asarray(members)
        sizes.append(m.size)
        durs.append(float(t[m].max() - t[m].min()))
        parts.append(int(np.unique(nid[m]).size))
        t_start.append(float(t[m].min()))
    t_start = np.sort(np.asarray(t_start))
    waits = np.diff(t_start) if t_start.size > 1 else np.zeros(0)
    return (np.asarray(sizes, dtype=np.int64), np.asarray(durs),
            np.asarray(parts, dtype=np.int64), waits)


# =============================================================================
# Discrete power-law MLE (Clauset-Shalizi-Newman), alternatives, Vuong LR
# =============================================================================


def hurwitz_zeta(a: float, q: float, N: int = 24) -> float:
    """Euler-Maclaurin Hurwitz zeta, adequate for a in (1, 8], q >= 1."""
    s = 0.0
    for k in range(N):
        s += (q + k) ** (-a)
    qn = q + N
    s += qn ** (1.0 - a) / (a - 1.0)
    s += 0.5 * qn ** (-a)
    s += a / 12.0 * qn ** (-a - 1.0)
    s -= a * (a + 1.0) * (a + 2.0) / 720.0 * qn ** (-a - 3.0)
    return s


def fit_powerlaw_discrete(s: np.ndarray, xmin: int) -> Tuple[float, float]:
    """Return (alpha_hat, loglik_per_event_sum) for tail s >= xmin."""
    tail = s[s >= xmin].astype(float)
    n = tail.size
    sum_log = float(np.sum(np.log(tail)))

    def nll(alpha: float) -> float:
        return alpha * sum_log + n * math.log(hurwitz_zeta(alpha, xmin))

    lo, hi = 1.005, 8.0
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c, d = b - gr * (b - a), a + gr * (b - a)
    fc, fd = nll(c), nll(d)
    for _ in range(80):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc = nll(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd = nll(d)
    alpha = 0.5 * (a + b)
    return alpha, -nll(alpha)


def pl_ccdf(svals: np.ndarray, alpha: float, xmin: int) -> np.ndarray:
    z = hurwitz_zeta(alpha, xmin)
    return np.array([hurwitz_zeta(alpha, float(sv)) / z for sv in svals])


def ks_powerlaw(s: np.ndarray, alpha: float, xmin: int) -> float:
    tail = np.sort(s[s >= xmin])
    n = tail.size
    uniq = np.unique(tail)
    model_ccdf = pl_ccdf(uniq, alpha, xmin)
    emp_ccdf = np.array([np.sum(tail >= u) for u in uniq], dtype=float) / n
    return float(np.max(np.abs(emp_ccdf - model_ccdf)))


def scan_xmin(s: np.ndarray, min_tail: int) -> Tuple[int, float, float, int]:
    """KS-optimal xmin. Returns (xmin, alpha, KS, n_tail)."""
    uniq = np.unique(s)
    best = None
    for xm in uniq:
        n_tail = int(np.sum(s >= xm))
        if n_tail < min_tail:
            break
        alpha, _ = fit_powerlaw_discrete(s, int(xm))
        ks = ks_powerlaw(s, alpha, int(xm))
        if best is None or ks < best[2]:
            best = (int(xm), alpha, ks, n_tail)
    if best is None:
        xm = int(uniq[0])
        alpha, _ = fit_powerlaw_discrete(s, xm)
        return xm, alpha, ks_powerlaw(s, alpha, xm), int(np.sum(s >= xm))
    return best


def fit_lognormal_tail(s: np.ndarray, xmin: int) -> Tuple[float, float, np.ndarray]:
    """Tail-conditioned continuous log-normal MLE. Returns (mu, sigma, pointwise ll)."""
    tail = s[s >= xmin].astype(float)
    ln = np.log(tail)
    mu0, sg0 = float(np.mean(ln)), max(float(np.std(ln)), 1e-3)
    lx = math.log(xmin)

    def norm_tail(mu, sg):
        # P(S >= xmin) = 0.5 * erfc((ln xmin - mu)/(sg sqrt2))
        return max(0.5 * _erfc((lx - mu) / (sg * math.sqrt(2.0))), 1e-300)

    def pointwise(mu, sg):
        return (-np.log(tail * sg * math.sqrt(2.0 * math.pi))
                - (ln - mu) ** 2 / (2.0 * sg * sg)
                - math.log(norm_tail(mu, sg)))

    def nll(p):
        mu, sg = p[0], abs(p[1]) + 1e-6
        return -float(np.sum(pointwise(mu, sg)))

    if _SCIPY:
        res = minimize(nll, np.array([mu0, sg0]), method="Nelder-Mead",
                       options={"xatol": 1e-6, "fatol": 1e-8, "maxiter": 2000})
        mu, sg = res.x[0], abs(res.x[1]) + 1e-6
    else:  # coarse grid + refinement fallback
        mu, sg, best = mu0, sg0, nll([mu0, sg0])
        for _ in range(3):
            for m_try in np.linspace(mu - 2, mu + 2, 21):
                for s_try in np.linspace(max(sg * 0.3, 1e-3), sg * 3, 21):
                    v = nll([m_try, s_try])
                    if v < best:
                        best, mu, sg = v, m_try, s_try
    return mu, sg, pointwise(mu, sg)


def fit_exponential_tail(s: np.ndarray, xmin: int) -> Tuple[float, np.ndarray]:
    tail = s[s >= xmin].astype(float)
    lam = 1.0 / max(float(np.mean(tail - xmin)), 1e-12)
    ll = np.log(lam) - lam * (tail - xmin)
    return lam, ll


def vuong(ll1: np.ndarray, ll2: np.ndarray) -> Tuple[float, float]:
    """Vuong LR test. Positive R favours model 1. Returns (R, p_two_sided)."""
    d = ll1 - ll2
    n = d.size
    R = float(np.sum(d))
    sd = float(np.std(d))
    if sd < 1e-12 or n < 2:
        return R, 1.0
    p = float(_erfc(abs(R) / (math.sqrt(2.0 * n) * sd)))
    return R, p


def pl_pointwise_ll(s: np.ndarray, alpha: float, xmin: int) -> np.ndarray:
    tail = s[s >= xmin].astype(float)
    return -alpha * np.log(tail) - math.log(hurwitz_zeta(alpha, xmin))


# =============================================================================
# Main analysis
# =============================================================================


def analyse_condition(
    sizes: np.ndarray, min_tail: int
) -> Dict[str, object]:
    out: Dict[str, object] = {
        "n_avalanches": int(sizes.size),
        "s_max": int(sizes.max()) if sizes.size else 0,
        "s_mean": float(sizes.mean()) if sizes.size else 0.0,
    }
    if sizes.size < min_tail or sizes.max() <= sizes.min():
        out["verdict"] = "insufficient"
        return out
    xmin, alpha, ks, n_tail = scan_xmin(sizes, min_tail)
    out.update({"xmin": xmin, "alpha": round(alpha, 4),
                "KS": round(ks, 4), "n_tail": n_tail})
    ll_pl = pl_pointwise_ll(sizes, alpha, xmin)
    mu, sg, ll_ln = fit_lognormal_tail(sizes, xmin)
    lam, ll_ex = fit_exponential_tail(sizes, xmin)
    R_ln, p_ln = vuong(ll_pl, ll_ln)
    R_ex, p_ex = vuong(ll_pl, ll_ex)
    out.update({"lognorm_mu": round(mu, 4), "lognorm_sigma": round(sg, 4),
                "exp_lambda": round(lam, 5),
                "LR_pl_vs_ln": round(R_ln, 3), "p_pl_vs_ln": round(p_ln, 4),
                "LR_pl_vs_exp": round(R_ex, 3), "p_pl_vs_exp": round(p_ex, 4)})
    # Conservative verdict logic: a power-law claim requires the PL to be
    # significantly favoured over BOTH alternatives AND >= ~1.5 decades of
    # tail. Otherwise report the honest state of the evidence.
    decades = math.log10(sizes.max() / max(xmin, 1))
    if n_tail < max(min_tail, 100):
        verdict = "insufficient"
    elif decades < 1.5:
        verdict = "tail_too_narrow_for_exponent_claim"
    elif R_ln > 0 and p_ln < 0.1 and R_ex > 0 and p_ex < 0.1:
        verdict = "powerlaw_favoured"
    elif R_ln < 0 and p_ln < 0.1:
        verdict = "lognormal_favoured"
    else:
        verdict = "indistinguishable"
    out["tail_decades"] = round(decades, 2)
    out["verdict"] = verdict
    return out


def make_ccdf_plot(sizes: np.ndarray, res: Dict, path: Path, title: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    uniq = np.unique(sizes)
    ccdf = np.array([np.sum(sizes >= u) for u in uniq], dtype=float) / sizes.size
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    ax.loglog(uniq, ccdf, "o", ms=4, mfc="none", label="data (CCDF)")
    if "alpha" in res:
        xm, al = res["xmin"], res["alpha"]
        sv = uniq[uniq >= xm].astype(float)
        frac = np.sum(sizes >= xm) / sizes.size
        ax.loglog(sv, frac * pl_ccdf(sv, al, xm), "-",
                  label=f"PL alpha={al:.2f}, xmin={xm}")
        mu, sg = res["lognorm_mu"], res["lognorm_sigma"]
        z0 = 0.5 * _erfc((math.log(xm) - mu) / (sg * math.sqrt(2)))
        ln_ccdf = np.array([0.5 * _erfc((math.log(v) - mu) / (sg * math.sqrt(2)))
                            for v in sv]) / max(z0, 1e-300)
        ax.loglog(sv, frac * ln_ccdf, "--", label="log-normal")
    ax.set_xlabel("avalanche size s")
    ax.set_ylabel("P(S >= s)")
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main(argv: Optional[Sequence[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Spatio-temporal avalanche clustering "
                                            "and MLE analysis for compassV2.1 event logs.")
    p.add_argument("--runs_dir", required=True, help="directory containing run outputs")
    p.add_argument("--out_dir", default="avalanche_analysis")
    p.add_argument("--channels", default="field,angle")
    p.add_argument("--t_link_T0", default="1,2,4",
                   help="comma list of causal link windows in units of T0 (sensitivity sweep)")
    p.add_argument("--r_link_rnn", type=float, default=1.05,
                   help="spatial link radius in units of r_nn")
    p.add_argument("--group_by", default="geometry,field_mode,Q",
                   help="metadata keys defining a physical condition for seed aggregation")
    p.add_argument("--min_tail", type=int, default=100,
                   help="minimum tail events for any exponent fit")
    p.add_argument("--make_plot", action="store_true")
    args = p.parse_args(argv)

    runs_dir = Path(args.runs_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    channels = [c.strip() for c in args.channels.split(",")]
    t_links = [float(x) for x in args.t_link_T0.split(",")]
    group_by = [g.strip() for g in args.group_by.split(",")]

    runs = discover_runs(runs_dir)
    if not runs:
        print("No complete runs found.")
        return 1
    print(f"Found {len(runs)} runs.")

    # Group runs by physical condition.
    groups: Dict[Tuple, List[Dict]] = defaultdict(list)
    for r in runs:
        groups[condition_key(r["meta"], group_by)].append(r)

    summary_rows: List[Dict] = []
    for ckey, members in sorted(groups.items(), key=lambda kv: str(kv[0])):
        cname = "_".join(str(k) for k in ckey)
        seeds = [m["meta"]["config"]["seed"] for m in members]
        print(f"\nCondition {ckey}: {len(members)} run(s), seeds {seeds}")

        # Per-run precomputation: positions, neighbours, T0, events.
        prepared = []
        for m in members:
            st = np.load(m["state"], allow_pickle=True)
            xs = np.asarray(st["xs"], float)
            ys = np.asarray(st["ys"], float)
            r_nn = float(st["r_nn"])
            T0 = float(m["meta"]["derived"]["T0_s"])
            neigh = neighbour_lists(xs, ys, args.r_link_rnn * r_nn)
            ev_t, ev_id, ev_ch = [], [], []
            with open(m["events"], newline="") as fh:
                for row in csv.DictReader(fh):
                    ev_t.append(float(row["t_s"]))
                    ev_id.append(int(row["needle_id"]))
                    ev_ch.append(row["channel"])
            prepared.append({
                "tag": m["tag"], "T0": T0, "neigh": neigh,
                "t": np.asarray(ev_t), "nid": np.asarray(ev_id, dtype=np.int64),
                "ch": np.asarray(ev_ch),
            })

        for channel in channels:
            for tl in t_links:
                all_sizes, all_durs, all_parts, all_waits = [], [], [], []
                prov = []
                for pr in prepared:
                    mask = pr["ch"] == channel
                    s, d, pa, w = cluster_events(
                        pr["t"][mask], pr["nid"][mask], pr["neigh"], tl * pr["T0"])
                    all_sizes.append(s)
                    all_durs.append(d)
                    all_parts.append(pa)
                    all_waits.append(w)
                    prov.extend([pr["tag"]] * s.size)
                sizes = np.concatenate(all_sizes) if all_sizes else np.zeros(0, dtype=np.int64)
                res = analyse_condition(sizes, args.min_tail)
                res_row = {"condition": cname, "channel": channel,
                           "t_link_T0": tl, "n_runs": len(members), **res}
                summary_rows.append(res_row)
                print(f"  [{channel}, t_link={tl}T0] n={res['n_avalanches']} "
                      f"s_max={res['s_max']} verdict={res['verdict']}"
                      + (f" alpha={res['alpha']} (xmin={res['xmin']}, "
                         f"n_tail={res['n_tail']})" if "alpha" in res else ""))

                np.savez_compressed(
                    out_dir / f"avalanches_{cname}_{channel}_tlink{tl:g}.npz",
                    sizes=sizes,
                    durations=np.concatenate(all_durs) if all_durs else np.zeros(0),
                    participants=np.concatenate(all_parts) if all_parts else np.zeros(0),
                    waiting_times=np.concatenate(all_waits) if all_waits else np.zeros(0),
                    run_tags=np.asarray(prov),
                    condition=cname, channel=channel, t_link_T0=tl,
                )
                if args.make_plot and sizes.size:
                    make_ccdf_plot(sizes, res,
                                   out_dir / f"ccdf_{cname}_{channel}_tlink{tl:g}.png",
                                   f"{cname} | {channel} | t_link={tl} T0")

    if summary_rows:
        cols: List[str] = []
        for r in summary_rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
        with open(out_dir / "summary_avalanches.csv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(summary_rows)
        print(f"\nWrote {out_dir/'summary_avalanches.csv'} ({len(summary_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
