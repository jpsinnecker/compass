#!/usr/bin/env python3
"""
damping_sweep — V02
====================

Campanha de simulacao para o estudo central:
    "Como o coeficiente de amortecimento (damping) controla a forma do
    laco de histerese e a estatistica de avalanches em redes dipolares
    classicas, e isso depende da geometria (quadrada / triangular /
    honeycomb)?"

Novidades da V02 (em relacao a V01):
    --n_workers N : executa N rodadas em paralelo (processos
                    independentes; ganho ~linear no numero de nucleos).
                    Cada worker roda em um diretorio de trabalho proprio
                    para evitar colisao dos arquivos incidentais
                    (hysteresis_loop.csv) gravados pelo simulador no
                    diretorio corrente.
    manifest incremental : o manifest.csv e atualizado apos CADA rodada
                    concluida (nao apenas ao final) — uma interrupcao no
                    meio da campanha nao perde mais o registro.
    --resume 1    : retoma uma campanha interrompida — rodadas ja
                    registradas no manifest.csv sao puladas.

Este script NAO modifica compass_sim.py. Ele importa o modulo e chama
relax() programaticamente, varrendo:

    damping  x  seed  x  geometria

para o modo de campo 'hysteresis' (Etapa 1) e 'static' B=0 (Etapa 2,
relaxacao livre para estatistica de dominios).

Saida
-----
Para cada rodada, salva um CSV com o field_log completo:
    out_dir/hysteresis/<geom>_damp<idx>_seed<seed>.csv
    out_dir/relax/<geom>_damp<idx>_seed<seed>.csv      (colunas: t,B,M,S)

E um JSON de metadados por rodada (parametros fisicos derivados, Q, etc.):
    out_dir/meta/<geom>_damp<idx>_seed<seed>.json

Um manifest.csv no topo do out_dir resume todas as rodadas (uma linha
por rodada), atualizado incrementalmente.

Uso
---
    python3 damping_sweepV02.py --out_dir /home/jps/sweep_results \
        --n_seeds 5 --grid_n 30 --n_workers 8

    --quick_test 1   roda uma versao reduzida (poucos dampings/seeds,
                      t_sim curto) so para validar que o pipeline
                      funciona ponta-a-ponta antes da campanha completa.
"""

import os
import sys
import json
import time
import argparse
import tempfile
import itertools
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compass_sim as cs


# ────────────────────────────────────────────────────────────────────────
# Parametros fisicos fixos do estudo (mesma agulha "de mesa" em todas as
# rodadas; so geometria/R/damping/seed variam)
# ────────────────────────────────────────────────────────────────────────
R_DEFAULT          = 0.025      # m
NEEDLE_FRAC        = 0.80       # mesmo default do CLI
DT_FACTOR          = 0.05       # mesmo default do CLI
CUTOFF_DEFAULT     = 3.5        # mesmo default do CLI (em unidades de R? -> ver assinatura relax)
NOISE_DEFAULT      = 1.5        # rad, ruido inicial nos angulos (default make_grid)


def derive_needle_geometry(R):
    needle_len   = NEEDLE_FRAC * 2.0 * R
    needle_width = needle_len * 0.22
    thickness    = cs.NEEDLE_THICKNESS_DEFAULT
    return needle_len, needle_width, thickness


def compute_physical_params(R):
    """Replica exatamente o que main() faz para derivar inertia/moment
    a partir da geometria, para que o estudo seja fisicamente consistente
    entre as 3 geometrias (mesma agulha real, mesma densidade de aco)."""
    needle_len, needle_width, thickness = derive_needle_geometry(R)
    inertia = cs.compute_inertia_from_geometry(needle_len, needle_width, thickness)
    moment  = cs.compute_moment_from_geometry(needle_len, needle_width, thickness)
    return inertia, moment, needle_len, needle_width, thickness


def compute_Q(inertia, moment, damping, r_nn):
    """Replica o calculo interno de relax() para Q = omega0*I/b,
    usando B_ref dipolar entre vizinhos mais proximos como campo de escala."""
    B_ref  = cs.MU0_OVER_4PI * 2.0 * moment / r_nn**3
    omega0 = np.sqrt(moment * B_ref / inertia)
    Q      = omega0 * inertia / damping if damping > 0 else np.inf
    return Q, omega0, B_ref


def make_damping_grid(inertia, moment, r_nn, n_points=8, Q_min=0.05, Q_max=15.0):
    """Constroi a grade de damping log-espacada em Q, depois converte para
    valores de `damping` reais via a relacao Q = omega0*I/b."""
    B_ref  = cs.MU0_OVER_4PI * 2.0 * moment / r_nn**3
    omega0 = np.sqrt(moment * B_ref / inertia)
    Q_vals = np.geomspace(Q_max, Q_min, n_points)   # do mais subamortecido ao mais sobreamortecido
    damping_vals = omega0 * inertia / Q_vals
    return damping_vals, Q_vals, omega0


def run_one_hysteresis(geometry, R, damping, seed, t_sim, B_max, n_grid,
                        pbc=False, dt_factor=DT_FACTOR, use_gpu=False,
                        show_progress=False):
    np.random.seed(seed)
    inertia, moment, needle_len, needle_width, thickness = compute_physical_params(R)

    xs, ys, thetas, r_nn, Lx, Ly = cs.make_grid(
        N=n_grid, M=n_grid, geometry=geometry, noise=NOISE_DEFAULT, R=R)

    Q, omega0, B_ref = compute_Q(inertia, moment, damping, r_nn)

    phi_ext_deg = 0.0  # direcao do campo fixa; eixo de simetria escolhido por geometria nao importa p/ M_proj

    t0 = time.perf_counter()
    theta_f, omega_f, hist, n_frames, dt, stop_reason, field_log = cs.relax(
        thetas, xs, ys,
        t_sim=t_sim, dt_factor=dt_factor,
        inertia=inertia, damping=damping,
        cutoff=CUTOFF_DEFAULT, ext_field=(B_max, 0.0),
        moment=moment, field_mode='hysteresis',
        pbc=pbc, Lx=Lx if pbc else None, Ly=Ly if pbc else None,
        use_gpu=use_gpu, show_progress=show_progress,
        make_images=False,
    )
    wall_s = time.perf_counter() - t0

    meta = dict(
        geometry=geometry, R=R, n_grid=n_grid, seed=seed,
        damping=damping, Q=Q, omega0=omega0, B_ref=B_ref, B_max=B_max,
        inertia=inertia, moment=moment, r_nn=r_nn,
        needle_len=needle_len, needle_width=needle_width, thickness=thickness,
        t_sim=t_sim, dt=dt, dt_factor=dt_factor, n_steps=len(field_log),
        pbc=pbc, stop_reason=stop_reason, wall_seconds=wall_s,
        field_mode='hysteresis',
    )
    return field_log, meta


def run_one_relax(geometry, R, damping, seed, t_sim, n_grid,
                   pbc=False, dt_factor=DT_FACTOR, use_gpu=False,
                   show_progress=False):
    np.random.seed(seed)
    inertia, moment, needle_len, needle_width, thickness = compute_physical_params(R)

    xs, ys, thetas, r_nn, Lx, Ly = cs.make_grid(
        N=n_grid, M=n_grid, geometry=geometry, noise=NOISE_DEFAULT, R=R)

    Q, omega0, B_ref = compute_Q(inertia, moment, damping, r_nn)

    t0 = time.perf_counter()
    theta_f, omega_f, hist, n_frames, dt, stop_reason, field_log = cs.relax(
        thetas, xs, ys,
        t_sim=t_sim, dt_factor=dt_factor,
        inertia=inertia, damping=damping,
        cutoff=CUTOFF_DEFAULT, ext_field=(0.0, 0.0),
        moment=moment, field_mode='static',
        pbc=pbc, Lx=Lx if pbc else None, Ly=Ly if pbc else None,
        use_gpu=use_gpu, show_progress=show_progress,
        make_images=False,
    )
    wall_s = time.perf_counter() - t0

    # estatistica de dominios no estado final
    domain_labels, n_domains = cs.label_magnetic_domains(theta_f, tol_deg=15.0)
    sizes = np.bincount(domain_labels.ravel())
    sizes = sizes[sizes > 0]

    meta = dict(
        geometry=geometry, R=R, n_grid=n_grid, seed=seed,
        damping=damping, Q=Q, omega0=omega0, B_ref=B_ref,
        inertia=inertia, moment=moment, r_nn=r_nn,
        needle_len=needle_len, needle_width=needle_width, thickness=thickness,
        t_sim=t_sim, dt=dt, dt_factor=dt_factor, n_steps=len(field_log),
        pbc=pbc, stop_reason=stop_reason, wall_seconds=wall_s,
        field_mode='static_relax',
        n_domains=int(n_domains),
        domain_sizes=sizes.tolist(),
    )
    return field_log, meta


def save_field_log_csv(field_log, path):
    arr = np.array(field_log)  # colunas: t, B, M, S
    header = "t,B,M,S"
    np.savetxt(path, arr, delimiter=",", header=header, comments="")


# ────────────────────────────────────────────────────────────────────────
# V02: execucao por "jobs" — serial ou paralela
# ────────────────────────────────────────────────────────────────────────

def _worker_init(scratch_root):
    """Inicializador de cada processo worker: muda para um diretorio de
    trabalho proprio, para que os arquivos incidentais que o simulador
    grava no diretorio corrente (hysteresis_loop.csv) nao colidam entre
    workers concorrentes."""
    d = tempfile.mkdtemp(prefix="worker_", dir=scratch_root)
    os.chdir(d)


def run_job(job):
    """Executa UMA rodada (histerese ou relaxacao) e retorna a linha do
    manifest. Funcao de nivel de modulo (picklavel) para uso com
    ProcessPoolExecutor. Grava o CSV e o meta JSON da rodada."""
    tag = job["tag"]
    if job["stage"] == "hysteresis":
        field_log, meta = run_one_hysteresis(
            job["geometry"], job["R"], job["damping"], job["seed"],
            job["t_sim"], job["B_max"], job["n_grid"],
            pbc=job["pbc"], use_gpu=job["use_gpu"], show_progress=False)
        csv_path = os.path.join(job["hyst_dir"], tag + ".csv")
        meta_path = os.path.join(job["meta_dir"], tag + "_hyst.json")
    else:  # relax
        field_log, meta = run_one_relax(
            job["geometry"], job["R"], job["damping"], job["seed"],
            job["t_sim"], job["n_grid"],
            pbc=job["pbc"], use_gpu=job["use_gpu"], show_progress=False)
        csv_path = os.path.join(job["relax_dir"], tag + ".csv")
        meta_path = os.path.join(job["meta_dir"], tag + "_relax.json")

    save_field_log_csv(field_log, csv_path)
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    row = dict(
        tag=tag, stage=job["stage"], geometry=job["geometry"],
        damp_idx=job["damp_idx"], seed=job["seed"], damping=job["damping"],
        Q=meta["Q"], n_steps=meta["n_steps"],
        wall_seconds=meta["wall_seconds"],
        csv_path=csv_path, meta_path=meta_path,
    )
    return row


MANIFEST_FIELDS = ["tag", "stage", "geometry", "damp_idx", "seed",
                   "damping", "Q", "n_steps", "wall_seconds",
                   "csv_path", "meta_path"]


def _append_manifest_row(manifest_path, row):
    """Escrita INCREMENTAL do manifest: acrescenta a linha imediatamente
    apos a rodada concluir (com flush), para que uma interrupcao no meio
    da campanha nao perca o registro das rodadas ja feitas."""
    import csv as csv_module
    new_file = not os.path.exists(manifest_path)
    with open(manifest_path, "a", newline="") as f:
        writer = csv_module.DictWriter(f, fieldnames=MANIFEST_FIELDS)
        if new_file:
            writer.writeheader()
        writer.writerow(row)
        f.flush()


def _load_done_set(manifest_path):
    """Para --resume: le o manifest existente e retorna o conjunto de
    (stage, tag) ja concluidos."""
    import csv as csv_module
    done = set()
    if os.path.exists(manifest_path):
        with open(manifest_path, newline="") as f:
            for row in csv_module.DictReader(f):
                done.add((row["stage"], row["tag"]))
    return done


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out_dir', type=str, required=True)
    ap.add_argument('--grid_n', type=int, default=30,
                     help='N=M do grid (default 30 -> 900 agulhas)')
    ap.add_argument('--n_seeds', type=int, default=5)
    ap.add_argument('--n_dampings', type=int, default=8)
    ap.add_argument('--Q_min', type=float, default=0.05)
    ap.add_argument('--Q_max', type=float, default=15.0)
    ap.add_argument('--geometries', type=str, default='square,triangular,honeycomb')
    ap.add_argument('--R', type=float, default=R_DEFAULT)
    ap.add_argument('--B_max_factor', type=float, default=8.0,
                     help='B_max = factor * B_ref (campo dipolar entre vizinhos), '
                          'para garantir saturacao no pico do ciclo')
    ap.add_argument('--t_sim_periods', type=float, default=40.0,
                     help='t_sim em unidades de T0 = 2*pi/omega0 (escala natural). '
                          'Ciclo de histerese tem 5 segmentos; este e o tempo TOTAL do ciclo.')
    ap.add_argument('--pbc', type=int, default=0)
    ap.add_argument('--use_gpu', type=int, default=0)
    ap.add_argument('--skip_relax_stage', type=int, default=0,
                     help='Se 1, roda so a Etapa 1 (histerese), pula a Etapa 2 (relaxacao/dominios)')
    ap.add_argument('--n_workers', type=int, default=1,
                     help='Numero de rodadas executadas em paralelo (processos). '
                          '1 = serial (padrao). Ganho ~linear no numero de nucleos. '
                          'Evite combinar >1 com --use_gpu 1 (contencao da GPU).')
    ap.add_argument('--resume', type=int, default=0,
                     help='Se 1, retoma campanha interrompida: rodadas ja presentes '
                          'no manifest.csv do out_dir sao puladas.')
    ap.add_argument('--quick_test', type=int, default=0,
                     help='Se 1, sobrescreve n_seeds=1, n_dampings=3, t_sim_periods=8, '
                          'grid_n=12 -- so para validar o pipeline ponta-a-ponta')
    args = ap.parse_args()

    if args.quick_test:
        args.n_seeds = 1
        args.n_dampings = 3
        args.t_sim_periods = 8.0
        args.grid_n = 12
        print(">>> MODO QUICK_TEST: parametros reduzidos para validacao do pipeline <<<")

    geometries = args.geometries.split(',')

    out_dir = args.out_dir
    hyst_dir   = os.path.join(out_dir, 'hysteresis')
    relax_dir  = os.path.join(out_dir, 'relax')
    meta_dir   = os.path.join(out_dir, 'meta')
    for d in (hyst_dir, relax_dir, meta_dir):
        os.makedirs(d, exist_ok=True)

    manifest_path = os.path.join(out_dir, 'manifest.csv')

    # ── grade de damping: calculada uma vez com a geometria de referencia
    #    (square), e reaplicada nas outras geometrias. Isso significa que
    #    o MESMO conjunto de valores fisicos de damping [N.m.s/rad] e usado
    #    em todas as geometrias -- o que muda entre geometrias e r_nn (e
    #    portanto Q efetivo varia ligeiramente entre geometrias, o que sera
    #    registrado no meta de cada rodada e usado na analise, nao escondido).
    inertia_ref, moment_ref, _, _, _ = compute_physical_params(args.R)
    xs_ref, ys_ref, _, r_nn_ref, _, _ = cs.make_grid(
        N=args.grid_n, M=args.grid_n, geometry='square', noise=0.0, R=args.R)
    damping_vals, Q_nominal, omega0_ref = make_damping_grid(
        inertia_ref, moment_ref, r_nn_ref,
        n_points=args.n_dampings, Q_min=args.Q_min, Q_max=args.Q_max)

    B_ref_square = cs.MU0_OVER_4PI * 2.0 * moment_ref / r_nn_ref**3
    B_max = args.B_max_factor * B_ref_square
    T0_ref = 2.0 * np.pi / omega0_ref
    t_sim = args.t_sim_periods * T0_ref

    print(f"Grade de damping (Q nominal, geometria square como referencia):")
    for d, q in zip(damping_vals, Q_nominal):
        print(f"    damping={d:.4e}  N.m.s/rad   ->  Q_nominal={q:.3f}")
    print(f"B_max = {B_max:.5f} T  ({args.B_max_factor}x B_ref)")
    print(f"t_sim = {t_sim:.4f} s  ({args.t_sim_periods} periodos T0={T0_ref:.5f}s)")
    print()

    total_runs = len(geometries) * len(damping_vals) * args.n_seeds
    t_campaign_start = time.perf_counter()

    # ── V02: monta a lista completa de jobs (Etapas 1 e 2) ──────────────
    common = dict(R=args.R, t_sim=t_sim, B_max=B_max, n_grid=args.grid_n,
                  pbc=bool(args.pbc), use_gpu=bool(args.use_gpu),
                  hyst_dir=hyst_dir, relax_dir=relax_dir, meta_dir=meta_dir)
    jobs = []
    stages = ['hysteresis'] + ([] if args.skip_relax_stage else ['relax'])
    for stage in stages:
        for geometry in geometries:
            for damp_idx, damping in enumerate(damping_vals):
                for seed in range(args.n_seeds):
                    tag = f"{geometry}_damp{damp_idx:02d}_seed{seed:02d}"
                    jobs.append(dict(common, stage=stage, tag=tag,
                                     geometry=geometry, damp_idx=damp_idx,
                                     seed=seed, damping=damping))

    # ── V02: resume — pula jobs ja registrados no manifest ──────────────
    if args.resume:
        done = _load_done_set(manifest_path)
        n_before = len(jobs)
        jobs = [j for j in jobs if (j['stage'], j['tag']) not in done]
        print(f"Resume: {n_before - len(jobs)} rodadas ja concluidas "
              f"(puladas); {len(jobs)} restantes.")
    elif os.path.exists(manifest_path):
        # campanha nova sobre out_dir antigo: recomeca o manifest do zero
        os.remove(manifest_path)

    n_total = len(jobs)
    if n_total == 0:
        print("Nada a fazer (todas as rodadas ja concluidas).")
        return

    # ── V02: execucao — serial ou paralela — com manifest incremental ───
    scratch_root = os.path.join(out_dir, '_scratch')
    os.makedirs(scratch_root, exist_ok=True)

    if args.n_workers <= 1:
        # serial: roda no diretorio de trabalho atual (comportamento V01)
        for i, job in enumerate(jobs, 1):
            print(f"[{i}/{n_total}] {job['stage'].upper():10s} {job['tag']}  "
                  f"(damping={job['damping']:.3e})  ...", end=' ', flush=True)
            row = run_job(job)
            _append_manifest_row(manifest_path, row)
            print(f"ok  ({row['wall_seconds']:.1f}s, Q={row['Q']:.3f})")
    else:
        print(f"Executando com {args.n_workers} workers em paralelo...")
        n_done = 0
        with ProcessPoolExecutor(max_workers=args.n_workers,
                                 initializer=_worker_init,
                                 initargs=(scratch_root,)) as ex:
            futures = {ex.submit(run_job, job): job for job in jobs}
            for fut in as_completed(futures):
                job = futures[fut]
                row = fut.result()   # propaga excecoes do worker
                _append_manifest_row(manifest_path, row)
                n_done += 1
                print(f"[{n_done}/{n_total}] {job['stage'].upper():10s} "
                      f"{job['tag']}  ok  ({row['wall_seconds']:.1f}s, "
                      f"Q={row['Q']:.3f})", flush=True)

    total_wall = time.perf_counter() - t_campaign_start
    print()
    print(f"Campanha concluida em {total_wall/60.0:.1f} min.")
    print(f"Manifest: {manifest_path}")


if __name__ == '__main__':
    main()
