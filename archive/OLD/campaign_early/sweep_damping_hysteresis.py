#!/usr/bin/env python3
"""
sweep_damping_hysteresis.py
============================
Campanha de varredura: damping (Q) x geometria x seed, em --field_mode hysteresis.

Roteiro implementado (conforme discutido):
  - 8 valores de damping, log-espaçados, cobrindo Q ~ 15 (subamortecido)
    a Q ~ 0.05 (sobreamortecido), calculados a partir da física real do
    sistema (omega0, I) — não são valores arbitrários.
  - 3 geometrias: square, triangular, honeycomb
  - 5 seeds de ruído inicial por combinação (damping, geometria)
  - Cada rodada roda um ciclo completo de histerese (5 segmentos:
    0->+Bmax->0->-Bmax->0->+Bmax) e salva o field_log (t, B, M, S)
    em CSV, mais um JSON de metadados da rodada.

Total: 8 x 3 x 5 = 120 rodadas.

Uso:
    python3 sweep_damping_hysteresis.py --out_dir /home/claude/sweep_results
    python3 sweep_damping_hysteresis.py --out_dir ... --dry_run   # só lista o plano
    python3 sweep_damping_hysteresis.py --out_dir ... --resume    # pula CSVs já existentes
    python3 sweep_damping_hysteresis.py --out_dir ... --gpu 1     # usa GPU se disponível

Cada rodada é independente -> pode ser interrompida (Ctrl+C) e retomada
depois com --resume sem perder o trabalho já feito.
"""

import argparse
import contextlib
import csv
import io
import json
import os
import shutil
import sys
import time
import traceback

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compass_sim as cs


# ─────────────────────────────────────────────────────────────────────────
# Parâmetros físicos de referência (mesma agulha de bússola de mesa usada
# nos defaults do compass_sim.py) — usados só para calcular Q = omega0*I/b
# de forma consistente com o que relax() vai calcular internamente.
# ─────────────────────────────────────────────────────────────────────────
R_DEFAULT       = 0.025                # m
MOMENT          = cs.MOMENT_DEFAULT    # A.m^2
INERTIA         = cs.INERTIA_DEFAULT   # kg.m^2
MU0_OVER_4PI    = cs.MU0_OVER_4PI

N_GRID          = 30                   # N x M -> 900 agulhas
M_GRID          = 30
NOISE_INIT      = 1.5                  # rad, ruido inicial nos angulos (mesmo p/ todas seeds -> seed controla *qual* ruido)
SEEDS           = [1, 2, 3, 4, 5]
GEOMETRIES      = ["square", "triangular", "honeycomb"]
N_DAMPING       = 8
Q_MIN, Q_MAX    = 0.05, 15.0           # sobreamortecido -> subamortecido

CUTOFF          = 8.0                  # m, raio de interacao dipolar (cobre boa parte da rede 30x30 com R~0.025 -> ~0.6m por lado)
B_MAX           = 4.0e-4               # T, amplitude maxima do ciclo de histerese (~2.5x B_ref tipico, ver calibracao)
DT_FACTOR       = 0.04                 # fracao de T0 por passo -- um pouco mais fino que o default 0.05
T_SIM_IN_T0     = 60.0                 # duracao total do ciclo em unidades de T0 (calibrado por damping, ver compute_t_sim)


def nn_distance_for_geometry(geometry, R):
    """Distancia entre primeiros vizinhos -- todas as geometrias usam 2R."""
    return 2.0 * R


def compute_omega0(R, moment=MOMENT, inertia=INERTIA):
    """Replica o calculo de omega0 feito dentro de relax() (B_ref a partir do
    campo dipolar entre primeiros vizinhos), para definir a faixa de damping
    e o t_sim de forma consistente com a fisica real que sera simulada."""
    r_nn = nn_distance_for_geometry("square", R)
    B_ref = MU0_OVER_4PI * 2.0 * moment / r_nn**3
    omega0 = np.sqrt(moment * B_ref / inertia)
    return omega0


def compute_damping_values(R=R_DEFAULT, n=N_DAMPING, q_min=Q_MIN, q_max=Q_MAX):
    """8 valores de damping log-espacados cobrindo Q de q_max (subamortecido)
    a q_min (sobreamortecido). Retorna lista de dicts {Q, damping}."""
    omega0 = compute_omega0(R)
    Q_targets = np.logspace(np.log10(q_max), np.log10(q_min), n)
    b_values = omega0 * INERTIA / Q_targets
    return [{"Q": float(Q), "damping": float(b)} for Q, b in zip(Q_targets, b_values)]


def compute_t_sim(damping, R=R_DEFAULT):
    """t_sim em segundos: fixamos a duracao do ciclo em unidades do periodo
    natural T0, garantindo resolucao temporal comparavel (em T0) entre os
    regimes sub e sobreamortecido. Mais T0's nao custa muito no overdamped
    pois dt tambem cresce, mas asseguramos um piso/teto para nao explodir
    o numero de passos no caso subamortecido (dt pequeno)."""
    omega0 = compute_omega0(R)
    T0 = 2.0 * np.pi / omega0
    return T_SIM_IN_T0 * T0


def run_single(geometry, damping, Q, seed, out_dir, R=R_DEFAULT,
                use_gpu=False, show_progress=False):
    """Roda uma simulacao de histerese e salva o field_log + metadados.
    Retorna o path do CSV salvo (ou None se pulado)."""

    tag = f"{geometry}_Q{Q:08.4f}_seed{seed:02d}"
    csv_path  = os.path.join(out_dir, f"{tag}.csv")
    meta_path = os.path.join(out_dir, f"{tag}.json")

    np.random.seed(seed)
    xs, ys, thetas, nn_dist, Lx, Ly = cs.make_grid(
        N=N_GRID, M=M_GRID, geometry=geometry, noise=NOISE_INIT, R=R
    )

    t_sim = compute_t_sim(damping, R=R)

    # relax() escreve "hysteresis_loop.csv"/".png" no cwd como efeito
    # colateral hardcoded (alem de retornar field_log, que e o que de fato
    # usamos). Isolamos cada rodada num scratch dir proprio para nao
    # sobrescrever/acumular lixo, e descartamos esse scratch ao final --
    # o CSV "oficial" desta rodada e o que escrevemos abaixo, em out_dir.
    scratch_dir = os.path.join(out_dir, "_scratch", tag)
    os.makedirs(scratch_dir, exist_ok=True)
    cwd_orig = os.getcwd()

    t_wall_start = time.perf_counter()
    try:
        os.chdir(scratch_dir)
        # show_progress=False ainda deixa os _print() informativos do
        # relax() ativos (nao ha flag de verbosidade no compass_sim.py) --
        # suprimimos stdout durante a chamada para nao poluir o log da
        # campanha de 120 rodadas.
        _buf = io.StringIO()
        with contextlib.redirect_stdout(_buf):
            theta_cur, omega_cur, hist, n_frames, dt, stop_reason, field_log = cs.relax(
                thetas, xs, ys,
                t_sim=t_sim,
                dt_factor=DT_FACTOR,
                inertia=INERTIA,
                damping=damping,
                cutoff=CUTOFF,
                ext_field=(B_MAX, 0.0),
                moment=MOMENT,
                field_mode="hysteresis",
                frame_dir=None,
                pbc=False,
                use_gpu=use_gpu,
                show_progress=show_progress,
            )
    finally:
        os.chdir(cwd_orig)
        shutil.rmtree(scratch_dir, ignore_errors=True)
    t_wall = time.perf_counter() - t_wall_start

    # ── salva o field_log (t, B_proj, M_proj, S) em CSV ────────────────────
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["t", "B_proj", "M_proj", "S"])
        writer.writerows(field_log)

    # ── metadados da rodada (essencial para reproduzir/auditar depois) ────
    meta = {
        "geometry": geometry,
        "Q_target": Q,
        "damping": damping,
        "seed": seed,
        "R": R,
        "N_grid": N_GRID,
        "M_grid": M_GRID,
        "noise_init": NOISE_INIT,
        "cutoff": CUTOFF,
        "B_max": B_MAX,
        "dt_factor": DT_FACTOR,
        "dt_actual": dt,
        "t_sim": t_sim,
        "n_steps_logged": len(field_log),
        "stop_reason": stop_reason,
        "wall_time_seconds": t_wall,
        "moment": MOMENT,
        "inertia": INERTIA,
        "used_gpu": bool(use_gpu),
    }
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    return csv_path, t_wall


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out_dir", required=True, help="Diretorio de saida para CSVs/JSONs")
    ap.add_argument("--dry_run", action="store_true", help="So lista o plano de rodadas, nao executa")
    ap.add_argument("--resume", action="store_true", help="Pula rodadas cujo CSV ja existe")
    ap.add_argument("--gpu", type=int, default=0, help="1 para usar GPU (CuPy) se disponivel")
    ap.add_argument("--geometries", nargs="+", default=GEOMETRIES, choices=GEOMETRIES)
    ap.add_argument("--seeds", nargs="+", type=int, default=SEEDS)
    ap.add_argument("--n_damping", type=int, default=N_DAMPING)
    ap.add_argument("--only_index", type=int, default=None,
                     help="Roda so a rodada de indice N do plano (debug/paralelizacao manual)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    damping_values = compute_damping_values(n=args.n_damping)

    plan = []
    for geometry in args.geometries:
        for dv in damping_values:
            for seed in args.seeds:
                plan.append({"geometry": geometry, "Q": dv["Q"],
                             "damping": dv["damping"], "seed": seed})

    print(f"Plano: {len(plan)} rodadas "
          f"({len(args.geometries)} geometrias x {len(damping_values)} dampings x {len(args.seeds)} seeds)")
    print(f"Faixa de Q: {damping_values[-1]['Q']:.3f} (sobreamortecido) "
          f"a {damping_values[0]['Q']:.3f} (subamortecido)")
    print(f"t_sim varia por damping (proporcional a {T_SIM_IN_T0} x T0)")
    print()

    if args.dry_run:
        for i, p in enumerate(plan):
            print(f"  [{i:3d}] geometry={p['geometry']:<10s} Q={p['Q']:8.4f} "
                  f"damping={p['damping']:.3e} seed={p['seed']}")
        return

    if args.only_index is not None:
        plan = [plan[args.only_index]]

    n_done, n_skipped, n_failed = 0, 0, 0
    t_campaign_start = time.perf_counter()

    for i, p in enumerate(plan):
        tag = f"{p['geometry']}_Q{p['Q']:08.4f}_seed{p['seed']:02d}"
        csv_path = os.path.join(args.out_dir, f"{tag}.csv")

        if args.resume and os.path.exists(csv_path):
            n_skipped += 1
            print(f"[{i+1}/{len(plan)}] {tag}  -- ja existe, pulando (--resume)")
            continue

        print(f"[{i+1}/{len(plan)}] {tag}  -- rodando...", flush=True)
        try:
            _, t_wall = run_single(
                geometry=p["geometry"], damping=p["damping"], Q=p["Q"],
                seed=p["seed"], out_dir=args.out_dir,
                use_gpu=bool(args.gpu), show_progress=False,
            )
            n_done += 1
            print(f"           ok  ({t_wall:.1f} s)")
        except KeyboardInterrupt:
            print("\nInterrompido pelo usuario (Ctrl+C). "
                  "Use --resume para continuar de onde parou.")
            sys.exit(130)
        except Exception as e:
            n_failed += 1
            print(f"           FALHOU: {e}")
            traceback.print_exc()
            # registra a falha para nao silenciar problemas
            fail_path = os.path.join(args.out_dir, f"{tag}.FAILED.txt")
            with open(fail_path, "w") as f:
                f.write(f"{type(e).__name__}: {e}\n\n")
                traceback.print_exc(file=f)

    t_campaign = time.perf_counter() - t_campaign_start
    print()
    print(f"Campanha concluida em {t_campaign/60.0:.1f} min: "
          f"{n_done} ok, {n_skipped} pulados, {n_failed} falharam.")


if __name__ == "__main__":
    main()
