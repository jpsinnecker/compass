"""
============================================================
compass_sim.py — Simulação de rede de agulhas de bússola
============================================================

Modela uma grade 2D de dipolos magnéticos clássicos (agulhas de
bússola) que interagem pelo campo que cada uma gera sobre as vizinhas.
A simulação busca o estado de equilíbrio por meio de uma dinâmica
amortecida inspirada na equação de Landau-Lifshitz-Gilbert (LLG).

Geometrias disponíveis
----------------------
  square      : grade retangular com espaçamentos independentes dx e dy.
  triangular  : grade triangular com offset de dx/2 em linhas ímpares.
  honeycomb   : rede tipo colmeia (grafeno), construída por linhas de
                átomos com espaçamento vertical alternado a/2 e a,
                onde a = dx é o comprimento do lado do hexágono.

Unidades do Sistema Internacional (SI)
---------------------------------------
  Todas as grandezas físicas estão em SI:

    Grandeza            Parâmetro       Unidade
    ─────────────────── ─────────────── ───────
    Momento magnético   --moment        A·m²
    Espaçamento horiz.  --dx            m
    Espaçamento vert.   --dy            m
    Campo externo       --B_ext         T  (Tesla)
    Direção do campo    --phi_ext       graus (0 = +x, 90 = +y)

  Valores típicos para agulhas de bússola de mesa (~5 cm de comprimento):
    moment ≈ 0.1  A·m²    (agulha pequena)  a  1.0  A·m² (agulha grande)
    dx, dy ≈ 0.03–0.10 m  (3 a 10 cm entre centros)
    B_ext  ≈ 50e-6 T      (campo terrestre ≈ 50 µT)
           ≈ 1e-3  T      (campo de imã de geladeira ≈ 1 mT a 5 cm)

  O campo dipolar gerado por uma agulha de moment=0.1 A·m² sobre sua
  vizinha a dx=0.05 m (ao longo do eixo) é:
    B = (μ₀/4π) · 2m/r³ = 1e-7 · 2·0.1 / 0.05³ ≈ 1.6 mT

Campo externo
-------------
  Especificado por intensidade B_ext (Tesla) e ângulo phi_ext (graus,
  medido a partir do eixo +x no sentido anti-horário).
  Uma seta dourada é desenhada em cada sítio da rede nas figuras.

Uso rápido
----------
  python compass_sim.py                                    # padrões SI
  python compass_sim.py --geometry honeycomb --N 12 --M 8
  python compass_sim.py --B_ext 50e-6 --phi_ext 0         # campo terrestre
  python compass_sim.py --B_ext 1e-3  --phi_ext 45        # 1 mT a 45°
  python compass_sim.py --video relaxacao.mp4             # gera vídeo MP4

Dependências: numpy, matplotlib, ffmpeg (para vídeo)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import argparse
import time

# ══════════════════════════════════════════════════════════════════════════════
# 1. CONSTANTES FÍSICAS (Sistema Internacional — SI)
# ══════════════════════════════════════════════════════════════════════════════

# μ₀/(4π) = 1×10⁻⁷  T·m/A  (permeabilidade do vácuo / 4π)
MU0_OVER_4PI = 1.0e-7   # T·m/A

# ── Parâmetros físicos padrão de uma agulha de bússola de mesa (~5 cm) ──────
#
# Momento magnético
MOMENT_DEFAULT = 0.1      # A·m²   (agulha de bússola de mesa típica)
#
# Massa e geometria da agulha (barra fina de aço, comprimento L, seção circular)
#   massa   m ≈ 0.5 g  = 5×10⁻⁴ kg
#   comprimento L ≈ 5 cm = 0.05 m
#   Momento de inércia de barra fina em torno do centro:
#       I = (1/12) · m · L²  =  (1/12) · 5×10⁻⁴ · (0.05)²  ≈ 1.04×10⁻⁷  kg·m²
INERTIA_DEFAULT = 1.0e-7  # kg·m²
#
# Amortecimento viscoso (resistência do ar ao giro da agulha).
# O Q efetivo depende do campo dominante B_eff = max(B_dipolar, B_externo):
#   Q = ω₀·I/b   onde  ω₀ = sqrt(m·B_eff/I)
#
# Com B_ext = 0.1 T:   ω₀ ≈ 316 rad/s  →  b = 8e-6 dá Q ≈ 4  (suave)
# Com B_ext = 0  :     ω₀ ≈  13 rad/s  →  b = 5e-8 dá Q ≈ 25 (bússola)
#
# O padrão 5e-8 é adequado para simular sem campo externo ou campo fraco.
# Para campos fortes (> 1 mT), use --damping 1e-6 a 1e-5 para movimento suave.
DAMPING_DEFAULT = 5.0e-8  # N·m·s/rad


# ══════════════════════════════════════════════════════════════════════════════
# 2. CAMPO DIPOLAR MAGNÉTICO (SI)
# ══════════════════════════════════════════════════════════════════════════════

def dipole_field_2d(rx, ry, theta_src, moment):
    """
    Calcula o campo magnético gerado por um dipolo no plano XY (unidades SI).

    Fórmula exata do campo de dipolo magnético em 3D, avaliada no plano z=0:

        B = (μ₀/4π) · [ 3(m̂·r̂)r̂ − m̂ ] / r³

    Em componentes cartesianas, com m = moment·(cos θ, sin θ, 0):

        Bx = (μ₀/4π) · [ 3(m·r)·rx / r⁵  −  mx / r³ ]
        By = (μ₀/4π) · [ 3(m·r)·ry / r⁵  −  my / r³ ]

    Unidades SI:
        moment  [A·m²]
        rx, ry  [m]
        Bx, By  [T]

    Parâmetros
    ----------
    rx, ry    : componentes do vetor r = (ponto receptor) − (dipolo fonte)  [m]
    theta_src : ângulo do dipolo fonte em relação ao eixo +x  [rad]
    moment    : módulo do momento de dipolo magnético da agulha fonte  [A·m²]

    Retorna
    -------
    (Bx, By)  : componentes do campo magnético no ponto (rx, ry)  [T]
    """
    r2 = rx**2 + ry**2
    if r2 < 1e-24:
        # evita divisão por zero quando o ponto coincide com a fonte
        return 0.0, 0.0

    r  = np.sqrt(r2)          # distância  [m]
    r5 = r2 * r2 * r          # r⁵  [m⁵]

    mx = moment * np.cos(theta_src)   # componente x do momento  [A·m²]
    my = moment * np.sin(theta_src)   # componente y do momento  [A·m²]

    mdotr = mx * rx + my * ry         # m · r  [A·m³]

    Bx = MU0_OVER_4PI * (3.0 * mdotr * rx / r5  -  mx / (r2 * r))
    By = MU0_OVER_4PI * (3.0 * mdotr * ry / r5  -  my / (r2 * r))
    return Bx, By   # [T]


# ══════════════════════════════════════════════════════════════════════════════
# 3. CAMPO TOTAL SOBRE UMA AGULHA (SI)
# ══════════════════════════════════════════════════════════════════════════════

def total_field_on(i, j, thetas, xs, ys, cutoff, moment):
    """
    Soma o campo dipolar de todas as agulhas vizinhas sobre a agulha (i, j).

    Para eficiência computacional, apenas agulhas dentro do raio `cutoff`
    são consideradas. Isso é fisicamente justificado porque o campo dipolar
    decai como 1/r³ — contribuições de agulhas distantes são desprezíveis.

    Parâmetros
    ----------
    i, j    : índices da agulha receptora na grade (N×M)
    thetas  : array 2D dos ângulos atuais de todas as agulhas  [rad]
    xs, ys  : arrays 2D das posições fixas das agulhas  [m]
    cutoff  : raio máximo de interação  [m]
    moment  : momento magnético de cada agulha  [A·m²]

    Retorna
    -------
    (Bx_tot, By_tot) : campo total na posição da agulha (i, j)  [T]
    """
    Bx_tot, By_tot = 0.0, 0.0
    N, M = thetas.shape
    xi, yi = xs[i, j], ys[i, j]

    for ni in range(N):
        for nj in range(M):
            if ni == i and nj == j:
                continue  # a agulha não interage com ela mesma

            rx = xi - xs[ni, nj]   # [m]
            ry = yi - ys[ni, nj]   # [m]
            dist = np.sqrt(rx**2 + ry**2)

            if dist > cutoff:
                continue  # fora do raio de corte: ignora

            bx, by = dipole_field_2d(rx, ry, thetas[ni, nj], moment)
            Bx_tot += bx
            By_tot += by

    return Bx_tot, By_tot   # [T]


# ══════════════════════════════════════════════════════════════════════════════
# 4. DINÂMICA INERCIAL (2ª Lei de Newton para rotação — sem atrito no pino)
# ══════════════════════════════════════════════════════════════════════════════

def relax(thetas, xs, ys, t_sim=2.0, dt_factor=0.05,
          inertia=INERTIA_DEFAULT, damping=DAMPING_DEFAULT,
          cutoff=3.5, ext_field=(0.0, 0.0), moment=MOMENT_DEFAULT,
          callback=None,
          frame_dir=None, frame_every=10,
          needle_len=0.042, needle_width=0.010, r_halo=None,
          B_ext=0.0, phi_ext_deg=0.0):
    """
    Integra a dinâmica inercial real de cada agulha sem atrito no pino.

    Modelo físico
    -------------
    Cada agulha é um corpo rígido com momento de inércia I girando em torno
    do pino (eixo z). A equação de movimento é a 2ª Lei de Newton rotacional:

        I · d²θ/dt²  =  τ_z(θ, t)  −  b · dθ/dt

    onde:
        I        [kg·m²]  : momento de inércia (barra fina: I = mL²/12)
        τ_z      [N·m]    : torque magnético planar = mx·By − my·Bx
        b        [N·m·s]  : amortecimento viscoso do ar
        dθ/dt    [rad/s]  : velocidade angular

    Atrito no pino = ZERO.

    Tempo de simulação
    ------------------
    `t_sim` define o tempo físico total [s]. O número de passos é calculado
    internamente: n_steps = ceil(t_sim / dt). A simulação pode terminar
    antes se convergir (ω_max < 0.1% de ω₀ por 50 passos consecutivos).

    Cronômetro nos frames
    ---------------------
    Cada PNG salvo exibe no canto superior esquerdo:
      - Tempo atual  t = X.XXXX s
      - Barra de progresso proporcional a t / t_sim
      - Tempo total  / X.XXX s

    Parâmetros
    ----------
    thetas      : array N×M de ângulos iniciais  [rad]
    xs, ys      : posições fixas das agulhas  [m]
    t_sim       : tempo total de simulação  [s]
    dt_factor   : fração do período natural T₀ usada como passo de tempo
    inertia     : momento de inércia de cada agulha  [kg·m²]
    damping     : coeficiente de amortecimento viscoso b  [N·m·s/rad]
    cutoff      : raio máximo de interação dipolar  [m]
    ext_field   : tupla (Bext_x, Bext_y)  [T]
    moment      : momento magnético de cada agulha  [A·m²]
    callback    : função opcional callback(step, thetas, omegas)
    frame_dir   : diretório onde salvar PNGs de frames (ou None)
    frame_every : intervalo em passos entre frames salvos
    needle_len  : tamanho das agulhas para renderização  [m]
    needle_width: largura das agulhas  [m]
    B_ext       : intensidade do campo externo para rótulo nos frames  [T]
    phi_ext_deg : direção do campo externo  [graus]

    Retorna
    -------
    theta_cur : array N×M — ângulos no estado final  [rad]
    omega_cur : array N×M — velocidades angulares  [rad/s]
    hist      : lista de tuplas (thetas, omegas) a cada 20 passos
    n_frames  : número de frames PNG salvos
    """
    import os

    N, M = thetas.shape
    Bext_x, Bext_y = ext_field
    n_frames = 0

    # ── distância entre primeiros vizinhos ─────────────────────────────────
    r_nn = np.inf
    for i in range(N):
        for j in range(M):
            for ni in range(N):
                for nj in range(M):
                    if ni == i and nj == j:
                        continue
                    d = np.sqrt((xs[i,j]-xs[ni,nj])**2 + (ys[i,j]-ys[ni,nj])**2)
                    if d < r_nn:
                        r_nn = d

    # ── campo efetivo, frequência natural, dt e n_steps ───────────────────
    B_ref     = MU0_OVER_4PI * 2.0 * moment / r_nn**3
    B_ext_mag = np.sqrt(Bext_x**2 + Bext_y**2)
    B_eff     = max(B_ref, B_ext_mag)
    omega0    = np.sqrt(moment * B_eff / inertia)
    T0        = 2.0 * np.pi / omega0
    dt        = dt_factor * T0
    # número de passos calculado a partir do tempo total pedido
    n_steps   = max(1, int(np.ceil(t_sim / dt)))
    Q         = omega0 * inertia / damping if damping > 0 else np.inf

    print(f"  Dinâmica inercial:")
    print(f"    r_nn   = {r_nn*100:.2f} cm")
    print(f"    B_ref  = {B_ref*1e3:.4f} mT  (dipolar entre vizinhos)")
    if B_ext_mag > 0:
        print(f"    B_ext  = {B_ext_mag*1e3:.4f} mT  (campo externo)")
    print(f"    B_eff  = {B_eff*1e3:.4f} mT  (campo dominante → define dt)")
    print(f"    ω₀     = {omega0:.2f} rad/s   T₀ = {T0:.5f} s")
    print(f"    dt     = {dt:.6f} s  ({dt_factor:.0%} de T₀)")
    print(f"    t_sim  = {t_sim:.3f} s  →  {n_steps} passos")
    print(f"    Q      = {Q:.1f}  "
          + ("(sub-amortecido — oscila)" if Q > 2
             else "(criticamente amortecido)" if Q > 0.5
             else "(super-amortecido — sem oscilações)"))
    print()

    # ── condições iniciais ─────────────────────────────────────────────────
    theta_cur = thetas.copy()
    omega_cur = np.zeros((N, M))
    hist      = [(theta_cur.copy(), omega_cur.copy())]

    if frame_dir is not None:
        os.makedirs(frame_dir, exist_ok=True)

    # ── funções auxiliares de formatação ──────────────────────────────────
    def _fmt_B(B):
        if B == 0:    return ""
        if B >= 0.1:  return f"B={B:.3f} T"
        if B >= 1e-4: return f"B={B*1e3:.3f} mT"
        return            f"B={B*1e6:.1f} µT"

    def _draw_clock(ax, wall_elapsed, needle_len, stop_label=None):
        """
        Desenha um cronômetro no canto superior esquerdo do eixo mostrando
        o tempo REAL de integração (tempo de relógio desde o início do loop).

        Elementos:
          - Texto grande: tempo decorrido em mm:ss.ss ou ss.ss s
          - Barra de progresso baseada em passos concluídos / n_steps
          - Rótulo de parada quando a simulação termina antes do tempo
          - Caixa semitransparente de fundo
        """
        xlim  = ax.get_xlim()
        ylim  = ax.get_ylim()
        xspan = xlim[1] - xlim[0]
        yspan = ylim[1] - ylim[0]

        # posição do painel — canto superior esquerdo
        px    = xlim[0] + 0.02 * xspan
        py    = ylim[1] - 0.03 * yspan

        bar_w = 0.30 * xspan
        bar_h = 0.018 * yspan
        bar_y = py - 0.085 * yspan

        extra_h = 0.030 * yspan if stop_label else 0.0

        # caixa de fundo semitransparente
        pad = needle_len * 0.25
        ax.add_patch(plt.Rectangle(
            (px - pad * 0.3, bar_y - pad * 0.6 - extra_h),
            bar_w + pad, 0.125 * yspan + pad + extra_h,
            facecolor='#080818', edgecolor='#3A3A6A',
            linewidth=0.8, alpha=0.80, zorder=19,
            transform=ax.transData))

        # ── barra de progresso baseada em passos (não em tempo de relógio) ─
        # usa n_frames como proxy do progresso — avança a cada frame salvo
        frac = min(n_frames / max(n_steps // frame_every, 1), 1.0)

        # trilha cinza
        ax.add_patch(plt.Rectangle(
            (px, bar_y), bar_w, bar_h,
            facecolor='#252540', edgecolor='none',
            zorder=20, transform=ax.transData))

        # preenchimento: verde → amarelo → vermelho
        if frac > 0:
            r_col = min(2.0 * frac, 1.0)
            g_col = min(2.0 * (1.0 - frac), 1.0)
            ax.add_patch(plt.Rectangle(
                (px, bar_y), bar_w * frac, bar_h,
                facecolor=(r_col, g_col, 0.15), edgecolor='none',
                zorder=21, transform=ax.transData))

        # borda
        ax.add_patch(plt.Rectangle(
            (px, bar_y), bar_w, bar_h,
            facecolor='none', edgecolor='#5A5A8A',
            linewidth=0.7, zorder=22, transform=ax.transData))

        # marcador branco
        if 0 < frac < 1.0:
            mx = px + bar_w * frac
            ax.plot([mx, mx], [bar_y, bar_y + bar_h],
                    color='white', lw=1.2, zorder=24,
                    transform=ax.transData)

        # ── texto do cronômetro ────────────────────────────────────────────
        # formata: mm:ss.ss se >= 60 s, senão ss.ssss s
        if wall_elapsed >= 60.0:
            mins = int(wall_elapsed // 60)
            secs = wall_elapsed - mins * 60
            time_str = f"{mins:02d}:{secs:05.2f}"
        else:
            time_str = f"{wall_elapsed:07.4f} s"

        ax.text(px + pad * 0.2, py,
                time_str,
                color='#E8E8FF', fontsize=11, fontweight='bold',
                fontfamily='monospace', va='top', ha='left',
                zorder=25, transform=ax.transData)

        # rótulo "tempo real" abaixo do número
        ax.text(px + pad * 0.2, py - 0.032 * yspan,
                "tempo de integração",
                color='#5555AA', fontsize=6,
                fontfamily='monospace', va='top', ha='left',
                zorder=25, transform=ax.transData)

        # rótulo de parada
        if stop_label:
            ax.text(px + pad * 0.2,
                    bar_y - pad * 0.5 - extra_h * 0.3,
                    stop_label,
                    color='#FFD700', fontsize=7, fontweight='bold',
                    fontfamily='monospace', va='top', ha='left',
                    zorder=25, transform=ax.transData)

    def _save_frame(step, th, om, stop_label=None):
        nonlocal n_frames
        wall_elapsed = time.time() - _wall_start   # tempo real desde o início [s]
        S      = np.abs(np.mean(np.exp(1j * th)))
        om_max = np.max(np.abs(om))
        b_str  = _fmt_B(B_ext)
        title  = (f"S = {S:.4f}   ω_max = {om_max:.2f} rad/s"
                  + (f"   {b_str} @ {phi_ext_deg:.0f}°" if b_str else ""))

        fig, ax = plot_state(th, xs, ys, title=title,
                             needle_len=needle_len, needle_width=needle_width,
                             r_halo=r_halo,
                             B_ext=B_ext, phi_ext_deg=phi_ext_deg)
        plt.tight_layout()

        # cronômetro com tempo real de integração
        _draw_clock(ax, wall_elapsed, needle_len, stop_label=stop_label)

        fpath = os.path.join(frame_dir, f"frame_{n_frames:05d}.png")
        plt.savefig(fpath, dpi=100, bbox_inches='tight', facecolor='#1A1A2E')
        plt.close(fig)
        n_frames += 1

    _wall_start = time.time()   # marca o início real da integração

    if frame_dir is not None:
        _save_frame(0, theta_cur, omega_cur)

    # ── thread de escuta de teclado (Ctrl+I = interrupção interativa) ──────
    # Ctrl+I no terminal envia o caractere ASCII \t (TAB, código 9).
    # Uma thread separada lê stdin em modo raw e seta uma flag compartilhada
    # quando detecta \t, sem interromper o processo principal.
    # Funciona em Linux e macOS. No Windows usa msvcrt.
    import threading, sys, os as _os

    _stop_flag = threading.Event()   # setado pela thread ou por S=1

    def _keyboard_listener():
        """
        Lê stdin caractere a caractere em modo raw.
        Ao detectar Ctrl+I (\\t, ASCII 9), seta _stop_flag.
        A thread termina automaticamente quando o processo principal encerra.
        """
        try:
            import tty, termios
            fd = sys.stdin.fileno()
            old = termios.tcgetattr(fd)
            try:
                tty.setraw(fd)
                while not _stop_flag.is_set():
                    ch = _os.read(fd, 1)
                    if ch == b'\t':          # Ctrl+I = TAB = ASCII 9
                        _stop_flag.set()
                        break
            finally:
                termios.tcsetattr(fd, termios.TCSADRAIN, old)
        except Exception:
            # stdin não é um terminal (redirecionamento, pipe, Windows sem
            # msvcrt) — desabilita silenciosamente a escuta de teclado
            pass

    # inicia thread como daemon para não bloquear o encerramento do programa
    _kb_thread = threading.Thread(target=_keyboard_listener, daemon=True)
    _kb_thread.start()

    # imprime instrução de uso
    print("  [Ctrl+I para interromper e salvar o vídeo agora]")

    # ── função auxiliar: calcula torques SI para um estado θ ───────────────
    def _torques(th):
        """Retorna array N×M de torques [N·m] para o estado angular th."""
        tau = np.zeros((N, M))
        for i in range(N):
            for j in range(M):
                Bx, By = total_field_on(i, j, th, xs, ys, cutoff, moment)
                Bx += Bext_x
                By += Bext_y
                mx = moment * np.cos(th[i, j])
                my = moment * np.sin(th[i, j])
                tau[i, j] = mx * By - my * Bx
        return tau

    # ── loop de integração Velocity-Verlet ────────────────────────────────
    tau_cur          = _torques(theta_cur)
    _converged_count = 0   # passos consecutivos com ω_max < tolerância
    _S1_count        = 0   # passos consecutivos com S ≥ 0.9999
    _stop_reason     = "tempo total atingido"

    for step in range(1, n_steps + 1):

        # ── verifica interrupção por teclado (Ctrl+I) ──────────────────────
        if _stop_flag.is_set():
            t_now = step * dt
            print(f"\n  Interrompido pelo usuário (Ctrl+I) em t={t_now:.4f}s"
                  f"  (passo {step}/{n_steps})")
            _stop_reason = "interrompido pelo usuário (Ctrl+I)"
            if frame_dir is not None:
                _save_frame(step, theta_cur, omega_cur,
                            stop_label="■ interrompido (Ctrl+I)")
            break

        # passo 1: atualiza θ
        accel_cur  = (tau_cur - damping * omega_cur) / inertia
        theta_new  = theta_cur + omega_cur * dt + 0.5 * accel_cur * dt**2
        theta_new  = (theta_new + np.pi) % (2.0 * np.pi) - np.pi

        # passo 2: calcula τ no novo estado
        tau_new = _torques(theta_new)

        # passo 3: atualiza ω (Velocity-Verlet implícito com amortecimento)
        b_half    = damping * dt / (2.0 * inertia)
        omega_new = (omega_cur * (1.0 - b_half)
                     + dt * (tau_cur + tau_new) / (2.0 * inertia)) \
                    / (1.0 + b_half)

        theta_cur = theta_new
        omega_cur = omega_new
        tau_cur   = tau_new

        if callback:
            callback(step, theta_cur.copy(), omega_cur.copy())

        if step % 20 == 0:
            hist.append((theta_cur.copy(), omega_cur.copy()))

        # ── calcula S e ω_max neste passo ─────────────────────────────────
        S_now     = np.abs(np.mean(np.exp(1j * theta_cur)))
        omega_max = np.max(np.abs(omega_cur))
        t_now     = step * dt

        # salva frame periódico
        if frame_dir is not None and step % frame_every == 0:
            _save_frame(step, theta_cur, omega_cur)
            if step % (frame_every * 10) == 0:
                print(f"  frame {n_frames-1:4d}  t={t_now:.4f}s  "
                      f"S={S_now:.4f}  ω_max={omega_max:.2f} rad/s")

        # ── condição 1: S = 1.00 (alinhamento total) ──────────────────────
        # A condição é satisfeita quando S ≥ 0.9999 por 30 passos
        # consecutivos, garantindo que não é flutuação momentânea.
        if S_now >= 0.9999:
            _S1_count += 1
        else:
            _S1_count = 0

        if _S1_count >= 30:
            print(f"\n  S = 1.00 atingido em t={t_now:.4f}s"
                  f"  (passo {step}/{n_steps})")
            _stop_reason = f"S = 1.00 atingido em t = {t_now:.4f} s"
            _stop_flag.set()
            if frame_dir is not None:
                _save_frame(step, theta_cur, omega_cur,
                            stop_label="★ S = 1.00  alinhamento total")
            break

        # ── condição 2: rede em repouso (ω_max → 0) ───────────────────────
        if omega_max < omega0 * 1e-3:
            _converged_count += 1
        else:
            _converged_count = 0

        if _converged_count >= 50:
            print(f"\n  Rede em repouso em t={t_now:.4f}s"
                  f"  (passo {step}/{n_steps})  S={S_now:.4f}")
            _stop_reason = f"rede em repouso em t = {t_now:.4f} s"
            _stop_flag.set()
            if frame_dir is not None:
                _save_frame(step, theta_cur, omega_cur,
                            stop_label="● rede em repouso")
            break

    _stop_flag.set()   # garante que a thread de teclado encerra
    return theta_cur, omega_cur, hist, n_frames, dt, _stop_reason

def next_available_path(path):
    """
    Se `path` não existe, retorna `path` sem alteração.
    Se `path` já existe, insere um sufixo numérico crescente antes da extensão:
        rapido.mp4  →  rapido0001.mp4  →  rapido0002.mp4  …

    Parâmetros
    ----------
    path : str — caminho desejado para o arquivo

    Retorna
    -------
    str — primeiro caminho disponível (que não existe ainda no disco)
    """
    import os
    if not os.path.exists(path):
        return path
    base, ext = os.path.splitext(path)
    n = 1
    while True:
        candidate = f"{base}{n:04d}{ext}"
        if not os.path.exists(candidate):
            return candidate
        n += 1


def render_video(frame_dir, output_path, fps=24, crf=20):
    """
    Agrega os PNGs numerados de `frame_dir` em um vídeo MP4 usando ffmpeg.

    Tenta três estratégias em ordem, parando na primeira que funcionar:

    1. H.264 (libx264) com -f mp4 explícito — melhor qualidade, máxima
       compatibilidade com players modernos.
    2. MPEG-4 (mpeg4) com -f mp4 — fallback para builds antigas do ffmpeg
       (conda macOS, por exemplo) que não incluem libx264.
    3. Codec nativo do ffmpeg (sem -c:v) — último recurso.

    O argumento -f mp4 força o formato de saída explicitamente, contornando
    o bug em versões antigas do ffmpeg que não inferem o formato pela extensão
    do arquivo de saída ("Unable to find a suitable output format").

    Parâmetros
    ----------
    frame_dir   : diretório com frame_00000.png, frame_00001.png, …
    output_path : caminho do arquivo MP4 de saída
    fps         : quadros por segundo (padrão 24)
    crf         : qualidade H.264 (0=lossless … 51=pior; padrão 20)

    Retorna
    -------
    True  se o ffmpeg gerou o arquivo com sucesso
    False se todas as tentativas falharam
    """
    import subprocess, shutil, os

    if not shutil.which('ffmpeg'):
        print("AVISO: ffmpeg não encontrado. Instale com:")
        print("  Ubuntu/Debian : sudo apt install ffmpeg")
        print("  macOS (brew)  : brew install ffmpeg")
        print("  conda         : conda install -c conda-forge ffmpeg")
        print("  Windows       : https://ffmpeg.org/download.html")
        return False

    input_pattern = os.path.join(frame_dir, "frame_%05d.png")
    vf = 'pad=ceil(iw/2)*2:ceil(ih/2)*2'

    strategies = [
        ("H.264 (libx264)",
         ['-c:v', 'libx264', '-preset', 'fast', '-crf', str(crf),
          '-pix_fmt', 'yuv420p']),
        ("MPEG-4 (mpeg4) — fallback para ffmpeg antigo",
         ['-c:v', 'mpeg4', '-q:v', '5',
          '-pix_fmt', 'yuv420p']),
        ("codec padrão do ffmpeg",
         ['-pix_fmt', 'yuv420p']),
    ]

    print(f"\nMontando vídeo MP4: {output_path}")
    print(f"  Fonte : {frame_dir}/frame_*.png   FPS={fps}   CRF={crf}")

    for desc, codec_args in strategies:
        cmd = [
            'ffmpeg', '-y',
            '-framerate', str(fps),
            '-i', input_pattern,
        ] + codec_args + [
            '-vf', vf,
            '-f', 'mp4',
            output_path,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            size_mb = os.path.getsize(output_path) / 1024**2
            print(f"  Codec : {desc}")
            print(f"  Salvo : {output_path}  ({size_mb:.1f} MB)")
            return True
        else:
            last_err = result.stderr.strip().splitlines()
            short_err = last_err[-1] if last_err else "(sem mensagem)"
            print(f"  [{desc}] falhou: {short_err}")

    print(f"\nERRO: ffmpeg não conseguiu gerar o vídeo. Erro completo:")
    print(result.stderr[-600:])
    return False
    """
    Agrega os PNGs numerados de `frame_dir` em um vídeo MP4 usando ffmpeg.

    Tenta três estratégias em ordem, parando na primeira que funcionar:

    1. H.264 (libx264) com -f mp4 explícito — melhor qualidade, máxima
       compatibilidade com players modernos.
    2. MPEG-4 (mpeg4) com -f mp4 — fallback para builds antigas do ffmpeg
       (conda macOS, por exemplo) que não incluem libx264.
    3. Codec nativo do ffmpeg (sem -c:v) — último recurso.

    O argumento -f mp4 força o formato de saída explicitamente, contornando
    o bug em versões antigas do ffmpeg que não inferem o formato pela extensão
    do arquivo de saída ("Unable to find a suitable output format").

    Parâmetros
    ----------
    frame_dir   : diretório com frame_00000.png, frame_00001.png, …
    output_path : caminho do arquivo MP4 de saída
    fps         : quadros por segundo (padrão 24)
    crf         : qualidade H.264 (0=lossless … 51=pior; padrão 20)

    Retorna
    -------
    True  se o ffmpeg gerou o arquivo com sucesso
    False se todas as tentativas falharam
    """
    import subprocess, shutil, os

    if not shutil.which('ffmpeg'):
        print("AVISO: ffmpeg não encontrado. Instale com:")
        print("  Ubuntu/Debian : sudo apt install ffmpeg")
        print("  macOS (brew)  : brew install ffmpeg")
        print("  conda         : conda install -c conda-forge ffmpeg")
        print("  Windows       : https://ffmpeg.org/download.html")
        return False

    input_pattern = os.path.join(frame_dir, "frame_%05d.png")

    # filtro de vídeo: garante dimensões pares (obrigatório para H.264/MPEG-4)
    vf = 'pad=ceil(iw/2)*2:ceil(ih/2)*2'

    # lista de estratégias: (descrição, argumentos extras de codec)
    strategies = [
        ("H.264 (libx264)",
         ['-c:v', 'libx264', '-preset', 'fast', '-crf', str(crf),
          '-pix_fmt', 'yuv420p']),
        ("MPEG-4 (mpeg4) — fallback para ffmpeg antigo",
         ['-c:v', 'mpeg4', '-q:v', '5',
          '-pix_fmt', 'yuv420p']),
        ("codec padrão do ffmpeg",
         ['-pix_fmt', 'yuv420p']),
    ]

    print(f"\nMontando vídeo MP4: {output_path}")
    print(f"  Fonte : {frame_dir}/frame_*.png   FPS={fps}   CRF={crf}")

    for desc, codec_args in strategies:
        cmd = [
            'ffmpeg', '-y',
            '-framerate', str(fps),
            '-i', input_pattern,
        ] + codec_args + [
            '-vf', vf,
            '-f', 'mp4',          # força formato de saída — corrige bug em versões antigas
            output_path,
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            size_mb = os.path.getsize(output_path) / 1024**2
            print(f"  Codec : {desc}")
            print(f"  Salvo : {output_path}  ({size_mb:.1f} MB)")
            return True
        else:
            # extrai apenas a última linha do erro para não poluir o terminal
            last_err = result.stderr.strip().splitlines()
            short_err = last_err[-1] if last_err else "(sem mensagem)"
            print(f"  [{desc}] falhou: {short_err}")

    # todas as estratégias falharam — mostra erro completo da última tentativa
    print(f"\nERRO: ffmpeg não conseguiu gerar o vídeo. Erro completo:")
    print(result.stderr[-600:])
    return False


# ══════════════════════════════════════════════════════════════════════════════
# 5. GERAÇÃO DA GRADE DE AGULHAS
# ══════════════════════════════════════════════════════════════════════════════

def make_grid(N=8, M=8, geometry='square', noise=1.5, R=0.025):
    """
    Cria as posições e os ângulos iniciais de uma grade N×M de agulhas.

    O parâmetro unificador é R [m] — o raio do círculo que envolve cada agulha.
    As posições são calculadas de modo que círculos adjacentes se toquem
    tangencialmente (distância entre centros = 2R), garantindo uma rede
    visualmente correta e sem sobreposição de círculos em qualquer geometria.

    Parâmetros
    ----------
    N, M      : número de linhas e colunas de agulhas
    geometry  : tipo de rede (ver abaixo)
    noise     : amplitude do ruído aleatório nos ângulos iniciais  [rad]
                  0   → todas apontam para +x
                  π   → orientação totalmente aleatória
    R         : raio do círculo que envolve cada agulha  [m]
                Padrão: 0.025 m = 2.5 cm
                A distância entre centros de agulhas adjacentes = 2R.

    Geometrias e espaçamentos resultantes
    --------------------------------------
    'square'
        Grade retangular. Vizinhos em 4 direções, todos a 2R.
            dx_real = 2R    (horizontal)
            dy_real = 2R    (vertical)
        Cada sítio tem 4 vizinhos mais próximos.

    'triangular'
        Grade triangular equilateral. Vizinhos em 6 direções, todos a 2R.
            dx_real = 2R          (espaçamento horizontal entre colunas)
            dy_real = R·√3        (espaçamento vertical entre linhas)
            offset  = R           (deslocamento de linhas ímpares)
        Cada sítio tem 6 vizinhos mais próximos.

    'honeycomb'
        Rede colmeia (grafeno). Lado do hexágono = 2R.
        Construída por linhas horizontais com padrão A-B alternado:

            linha par  (tipo A): offset_x = 0,        Δy_próx = R
            linha ímpar(tipo B): offset_x = R·√3,     Δy_próx = 2R

        Espaçamentos:
            entre colunas na mesma linha : 2R·√3
            Δy(A→B)                      : R
            Δy(B→A)                      : 2R   (período = 3R no total)

        Verificação: distância entre vizinhos A e B:
            d = √((R√3)² + R²) = √(3R²+R²) = 2R  ✓

        Cada sítio tem exatamente 3 vizinhos mais próximos (coordenação 3).
        Total de agulhas = N × M.

    Retorna
    -------
    xs, ys  : arrays 2D (N×M) com as coordenadas das agulhas  [m]
    thetas  : array 2D (N×M) com os ângulos iniciais  [rad]
    nn_dist : distância entre primeiros vizinhos  [m]  (= 2R em todos os casos)
    """
    s3 = np.sqrt(3.0)

    # ── grade quadrada ─────────────────────────────────────────────────────
    if geometry == 'square':
        d = 2.0 * R          # distância entre vizinhos
        xs = np.zeros((N, M))
        ys = np.zeros((N, M))
        for i in range(N):
            for j in range(M):
                xs[i, j] = j * d
                ys[i, j] = i * d
        thetas = noise * np.random.randn(N, M)
        return xs, ys, thetas, d

    # ── grade triangular equilateral ───────────────────────────────────────
    elif geometry == 'triangular':
        d      = 2.0 * R          # distância entre vizinhos
        dx_col = d                # espaçamento entre colunas
        dy_row = R * s3           # espaçamento entre linhas = R·√3 = d·(√3/2)
        offset = R                # offset de linhas ímpares = R = d/2
        xs = np.zeros((N, M))
        ys = np.zeros((N, M))
        for i in range(N):
            x_off = offset * (i % 2)
            for j in range(M):
                xs[i, j] = j * dx_col + x_off
                ys[i, j] = i * dy_row
        thetas = noise * np.random.randn(N, M)
        return xs, ys, thetas, d

    # ── grade honeycomb (colmeia) ──────────────────────────────────────────
    elif geometry == 'honeycomb':
        # lado do hexágono = 2R (= distância entre vizinhos adjacentes)
        # Linha tipo A (par):   offset_x = 0,      colunas a 2R·√3
        # Linha tipo B (ímpar): offset_x = R·√3,   colunas a 2R·√3
        # Δy(A→B) = R,  Δy(B→A) = 2R
        #
        # Verificação: A(0,0) → B(R·√3, R):
        #   d = √((R√3)² + R²) = √(4R²) = 2R  ✓
        d      = 2.0 * R          # distância entre primeiros vizinhos
        dx_col = d * s3           # espaçamento horizontal entre colunas = 2R√3
        xs_list, ys_list = [], []
        y = 0.0
        for row in range(N):
            if row % 2 == 0:      # linha tipo A
                x_offset = 0.0
                dy_next  = R      # Δy até a próxima linha (B)
            else:                  # linha tipo B
                x_offset = R * s3
                dy_next  = d      # Δy = 2R até a próxima linha (A)
            for col in range(M):
                xs_list.append(col * dx_col + x_offset)
                ys_list.append(y)
            y += dy_next

        xs     = np.array(xs_list).reshape(N, M)
        ys     = np.array(ys_list).reshape(N, M)
        thetas = noise * np.random.randn(N, M)
        return xs, ys, thetas, d

    else:
        raise ValueError(f"Geometria desconhecida: '{geometry}'. "
                         "Use 'square', 'triangular' ou 'honeycomb'.")


# ══════════════════════════════════════════════════════════════════════════════
# 6. DESENHO DE UMA AGULHA DE BÚSSOLA
# ══════════════════════════════════════════════════════════════════════════════

def draw_compass(ax, x, y, theta, length=0.42, width=0.10,
                 color_n='#FFFFFF', color_s='#2E6DB4',
                 edge='#1a1a1a', zorder=4):
    """
    Desenha uma agulha de bússola tradicional: losango bicolor.

    Geometria do losango (coordenadas locais, eixo longo ao longo de +x):
        vértice 0  (+half,  0)      → ponta norte  (polo +)
        vértice 1  (0,  +half_w)    → largura superior
        vértice 2  (−half,  0)      → ponta sul    (polo −)
        vértice 3  (0,  −half_w)    → largura inferior

    O losango é dividido em duas metades pelo centro (x, y):
        norte (branca) : triângulo [0, 1, centro, 3]
        sul   (azul)   : triângulo [2, 1, centro, 3]

    Após construir em coordenadas locais, aplica rotação por theta e
    translação para (x, y).

    Parâmetros
    ----------
    ax               : eixo matplotlib onde desenhar
    x, y             : posição do centro (pivô) da agulha
    theta            : ângulo de orientação (rad), medido do eixo +x
    length           : comprimento total da agulha
    width            : largura máxima (no centro do losango)
    color_n, color_s : cores das metades norte e sul
    edge             : cor da borda do losango
    zorder           : ordem de empilhamento para sobreposição correta
    """
    half   = length / 2.0
    half_w = width  / 2.0

    # vértices em coordenadas locais (não rotacionados)
    pts_local = np.array([
        [ half,     0.0    ],   # 0: ponta norte
        [ 0.0,      half_w ],   # 1: largura superior
        [-half,     0.0    ],   # 2: ponta sul
        [ 0.0,     -half_w ],   # 3: largura inferior
    ])

    # matriz de rotação 2×2 pelo ângulo theta
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s],
                  [s,  c]])
    # aplica rotação e translação
    pts = (R @ pts_local.T).T + np.array([x, y])

    # metade norte: polígono [ponta_norte, largura_sup, centro, largura_inf]
    north = plt.Polygon(
        [pts[0], pts[1], [x, y], pts[3]],
        closed=True, facecolor=color_n,
        edgecolor=edge, linewidth=0.5, zorder=zorder)

    # metade sul: polígono [ponta_sul, largura_sup, centro, largura_inf]
    south = plt.Polygon(
        [pts[2], pts[1], [x, y], pts[3]],
        closed=True, facecolor=color_s,
        edgecolor=edge, linewidth=0.5, zorder=zorder)

    ax.add_patch(south)
    ax.add_patch(north)

    # pino central (eixo de rotação)
    ax.plot(x, y, 'o', ms=2.0, color='#555555',
            markeredgecolor='#222222', markeredgewidth=0.4,
            zorder=zorder + 1)


# ══════════════════════════════════════════════════════════════════════════════
# 7. VISUALIZAÇÃO DO ESTADO DA REDE
# ══════════════════════════════════════════════════════════════════════════════

def _fmt_B(B):
    """Formata intensidade de campo em unidade legível (T, mT ou µT)."""
    if B == 0:      return "0 T"
    if B >= 0.1:    return f"{B:.4f} T"
    if B >= 1e-4:   return f"{B*1e3:.4f} mT"
    return              f"{B*1e6:.2f} µT"


def draw_ext_field_on_lattice(ax, xs, ys, B_ext, phi_ext_deg,
                              needle_len, color='#FFD700'):
    """
    Representa o campo externo uniforme na imagem de duas formas:

    1. SETAS NOS SÍTIOS — uma seta dourada centrada em cada agulha,
       apontando na direção do campo. Comprimento fixo = 35% de needle_len,
       sempre visível independente de B_ext. Renderizadas abaixo das
       agulhas (zorder=2) para servir de campo de fundo.

    2. PAINEL DE CAMPO — caixa no canto inferior direito da figura com:
         • Seta grande mostrando a direção (comprimento fixo = needle_len)
         • Texto: intensidade em SI + ângulo em graus
       O painel usa coordenadas de eixo (transform=ax.transAxes) para
       ficar sempre dentro da figura independente da escala da rede.

    Se B_ext = 0, nada é desenhado.

    Parâmetros
    ----------
    ax          : eixo matplotlib
    xs, ys      : arrays 2D de posições das agulhas  [m]
    B_ext       : intensidade do campo externo  [T]
    phi_ext_deg : direção do campo  [graus; 0=+x, 90=+y, anti-horário]
    needle_len  : comprimento de referência das agulhas  [m]
    color       : cor das setas e do texto do campo (padrão: dourado)
    """
    if B_ext <= 0:
        return

    phi      = np.deg2rad(phi_ext_deg)
    cos_phi  = np.cos(phi)
    sin_phi  = np.sin(phi)

    # ── 1. setas pequenas em cada sítio (campo de fundo) ──────────────────
    # Comprimento fixo = 35% de needle_len — sempre visível
    alen  = needle_len * 0.35
    dx_a  = alen * cos_phi
    dy_a  = alen * sin_phi
    N, M  = xs.shape

    for i in range(N):
        for j in range(M):
            x0 = xs[i, j] - dx_a / 2.0
            y0 = ys[i, j] - dy_a / 2.0
            ax.annotate(
                '', xy=(x0 + dx_a, y0 + dy_a), xytext=(x0, y0),
                arrowprops=dict(
                    arrowstyle='->', color=color,
                    lw=1.0, mutation_scale=6, alpha=0.50),
                zorder=2)

    # ── 2. painel de campo no canto inferior direito ───────────────────────
    # Usa coordenadas de eixo (0–1) para posição independente da escala.
    # A seta é desenhada em coordenadas de dados, mas ancorada em coords de eixo
    # usando um ponto de dados calculado a partir dos limites do eixo.

    # posição do centro da seta em coordenadas de dados
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    xspan = xlim[1] - xlim[0]
    yspan = ylim[1] - ylim[0]

    # canto inferior direito: 88% em x, 9% em y (em fração dos limites)
    cx_d = xlim[0] + 0.88 * xspan
    cy_d = ylim[0] + 0.09 * yspan

    # comprimento da seta de referência = needle_len (tamanho de uma agulha)
    slen = needle_len * 1.1
    x_tail = cx_d - slen * cos_phi / 2.0
    y_tail = cy_d - slen * sin_phi / 2.0
    x_head = cx_d + slen * cos_phi / 2.0
    y_head = cy_d + slen * sin_phi / 2.0

    # seta da direção do campo
    ax.annotate(
        '', xy=(x_head, y_head), xytext=(x_tail, y_tail),
        arrowprops=dict(
            arrowstyle='->', color=color, lw=2.5,
            mutation_scale=18),
        zorder=20)

    # caixa de fundo semitransparente atrás do painel
    pad   = needle_len * 0.55
    rect  = plt.Rectangle(
        (cx_d - slen * 0.7, cy_d - pad * 0.5),
        slen * 1.4, pad * 3.8,
        facecolor='#0A0A1A', edgecolor=color,
        linewidth=1.0, alpha=0.75, zorder=19,
        transform=ax.transData)
    ax.add_patch(rect)

    # ── texto: intensidade e direção ──────────────────────────────────────
    # linha 1: intensidade em SI
    # linha 2: ângulo em graus
    b_str  = _fmt_B(B_ext)
    phi_str = f"φ = {phi_ext_deg:.1f}°"

    # posição do texto: acima da seta dentro da caixa
    tx = cx_d
    ty = cy_d + pad * 1.4

    ax.text(tx, ty, b_str,
            color=color, fontsize=9, fontweight='bold',
            fontfamily='monospace', ha='center', va='bottom',
            zorder=21)
    ax.text(tx, ty - pad * 0.85, phi_str,
            color=color, fontsize=8,
            fontfamily='monospace', ha='center', va='bottom',
            zorder=21)

    # rótulo "B ext" acima
    ax.text(tx, cy_d - pad * 0.3, "B ext",
            color=color, fontsize=7, alpha=0.8,
            fontfamily='monospace', ha='center', va='top',
            zorder=21)


def plot_state(thetas, xs, ys, title="Rede de bussolas", show_order=True,
               needle_len=0.42, needle_width=0.10, r_halo=None,
               B_ext=0.0, phi_ext_deg=0.0):
    """
    Gera uma figura com o estado instantâneo da rede de agulhas.

    Elementos visuais
    -----------------
    - Halo colorido em torno de cada agulha (verde = vizinhos alinhados,
      vermelho = vizinhos anti-alinhados ou frustrados) — parâmetro de
      ordem LOCAL, calculado como média de cos(Δθ) com os vizinhos em
      grade (deslocamentos ±1 em i e j).
    - Agulhas desenhadas como losangos bicolores (branco = norte, azul = sul).
    - Setas douradas em cada sítio mostrando a direção do campo externo.
    - Painel no canto inferior direito com seta de direção, intensidade
      em SI e ângulo em graus.
    - Legenda com a convenção de cores.

    Parâmetros
    ----------
    thetas       : array N×M de ângulos das agulhas  [rad]
    xs, ys       : arrays N×M de posições  [m]
    title        : título da figura
    show_order   : se True, desenha os halos de parâmetro de ordem local
    needle_len   : comprimento das agulhas  [m]
    needle_width : largura das agulhas  [m]
    B_ext        : intensidade do campo externo  [T]
    phi_ext_deg  : direção do campo externo  [graus]

    Retorna
    -------
    fig, ax : objetos matplotlib
    """
    fig, ax = plt.subplots(figsize=(8, 8), facecolor='#1A1A2E')
    ax.set_facecolor('#16213E')
    N, M = thetas.shape

    # ── halos de parâmetro de ordem local ────────────────────────────────
    if show_order:
        for i in range(N):
            for j in range(M):
                align, count = 0.0, 0
                for di, dj in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                    ni_, nj_ = i + di, j + dj
                    if 0 <= ni_ < N and 0 <= nj_ < M:
                        align += np.cos(thetas[i, j] - thetas[ni_, nj_])
                        count += 1
                if count:
                    align /= count
                col    = cm.RdYlGn((align + 1.0) / 2.0)
                _r     = r_halo if r_halo is not None else needle_len * 0.58
                ax.add_patch(plt.Circle(
                    (xs[i, j], ys[i, j]), _r,
                    color=col, alpha=0.20, zorder=1))

    # ── agulhas ───────────────────────────────────────────────────────────
    for i in range(N):
        for j in range(M):
            draw_compass(ax, xs[i, j], ys[i, j], thetas[i, j],
                         length=needle_len, width=needle_width)

    # ── formatação dos eixos (deve vir ANTES do painel de campo) ──────────
    ax.set_aspect('equal')
    margin = needle_len * 1.6
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + margin * 1.5)
    ax.set_title(title, color='#ECF0F1', fontsize=11, pad=10,
                 fontfamily='monospace')
    ax.tick_params(colors='#7F8C8D')
    for spine in ax.spines.values():
        spine.set_edgecolor('#2C3E50')

    # ── campo externo (setas + painel) — depois de set_xlim/ylim ─────────
    draw_ext_field_on_lattice(ax, xs, ys, B_ext, phi_ext_deg, needle_len)

    # ── legenda ───────────────────────────────────────────────────────────
    patch_n = mpatches.Patch(facecolor='#FFFFFF', edgecolor='#555',
                              label='Polo Norte (branco)')
    patch_s = mpatches.Patch(facecolor='#2E6DB4', edgecolor='#555',
                              label='Polo Sul (azul)')
    handles = [patch_n, patch_s]
    if B_ext > 0:
        patch_b = mpatches.Patch(facecolor='#FFD700', edgecolor='#555',
                                  label=f'Campo ext.  {_fmt_B(B_ext)}  '
                                        f'φ={phi_ext_deg:.1f}°')
        handles.append(patch_b)
    ax.legend(handles=handles, loc='lower left',
              facecolor='#1A1A2E', edgecolor='#2C3E50',
              labelcolor='#ECF0F1', fontsize=8)

    return fig, ax


# ══════════════════════════════════════════════════════════════════════════════
# 8. ANIMAÇÃO DA RELAXAÇÃO
# ══════════════════════════════════════════════════════════════════════════════

def animate_relaxation(thetas_hist, xs, ys,
                       needle_len=0.42, needle_width=0.10,
                       B_ext=0.0, phi_ext_deg=0.0,
                       interval=80, save_gif=None):
    """
    Gera uma animação quadro a quadro da evolução temporal da rede.

    Cada quadro corresponde a um snapshot salvo durante a relaxação
    (um a cada 20 passos de integração).

    Parâmetros
    ----------
    thetas_hist  : lista de arrays N×M (saída do campo `hist` de `relax`)
    xs, ys       : posições das agulhas
    needle_len   : tamanho das agulhas
    needle_width : largura das agulhas
    B_ext        : intensidade do campo externo (para seta)
    phi_ext_deg  : direção do campo externo (graus)
    interval     : tempo entre quadros em ms (menor = mais rápido)
    save_gif     : caminho para salvar o GIF; None = exibe na tela

    Retorna
    -------
    fig, ani : objetos matplotlib (FuncAnimation)
    """
    fig, ax = plt.subplots(figsize=(7, 7), facecolor='#1A1A2E')
    ax.set_facecolor('#16213E')
    ax.set_aspect('equal')
    margin = needle_len * 1.6
    ax.set_xlim(xs.min() - margin, xs.max() + margin)
    ax.set_ylim(ys.min() - margin, ys.max() + margin * 2.5)

    N, M = thetas_hist[0].shape

    def update(frame):
        # limpa todos os artistas do quadro anterior
        for art in ax.lines + ax.patches + ax.collections:
            art.remove()
        thetas = thetas_hist[frame]
        for i in range(N):
            for j in range(M):
                draw_compass(ax, xs[i, j], ys[i, j], thetas[i, j],
                             length=needle_len, width=needle_width)
        draw_ext_field_on_lattice(ax, xs, ys, B_ext, phi_ext_deg, needle_len)
        ax.set_title(f"Relaxação dipolar — passo {frame * 20}",
                     color='#ECF0F1', fontsize=11, fontfamily='monospace')
        return []

    ani = FuncAnimation(fig, update, frames=len(thetas_hist),
                        interval=interval, blit=False)
    if save_gif:
        ani.save(save_gif, writer='pillow', fps=12, dpi=90)
        print(f"GIF salvo em: {save_gif}")
    return fig, ani


# ══════════════════════════════════════════════════════════════════════════════
# 9. GRÁFICO DO PARÂMETRO DE ORDEM GLOBAL
# ══════════════════════════════════════════════════════════════════════════════

def plot_order_parameter(thetas_hist, outpath, dt=None):
    """
    Plota a evolução temporal do parâmetro de ordem magnético global S(t).

        S(t) = |〈e^{iθ}〉|  ∈ [0, 1]

    S ≈ 1 → agulhas alinhadas;  S ≈ 0 → orientações aleatórias.

    Parâmetros
    ----------
    thetas_hist : lista de arrays N×M (histórico de `relax`, cada 20 passos)
    outpath     : caminho do arquivo de saída (PNG)
    dt          : passo de tempo [s] — se fornecido, eixo x em segundos reais;
                  se None, eixo x em índice de snapshot
    """
    order_params = [
        np.abs(np.mean(np.exp(1j * th)))
        for th in thetas_hist
    ]

    if dt is not None:
        # cada snapshot corresponde a 20 passos
        time_ax = np.arange(len(order_params)) * 20 * dt
        xlabel  = "Tempo  t  [s]"
    else:
        time_ax = np.arange(len(order_params))
        xlabel  = "Snapshot (a cada 20 passos)"

    fig, ax = plt.subplots(figsize=(8, 3.5), facecolor='#1A1A2E')
    ax.set_facecolor('#0F3460')
    ax.plot(time_ax, order_params, color='#E94560', lw=2)
    ax.set_xlabel(xlabel, color='#BDC3C7')
    ax.set_ylabel(r"Parâmetro de ordem $S = |\langle e^{i\theta}\rangle|$",
                  color='#BDC3C7')
    ax.set_title("Evolução do parâmetro de ordem magnético global",
                 color='#ECF0F1', fontfamily='monospace')
    ax.tick_params(colors='#7F8C8D')
    ax.set_ylim(0, 1.05)
    ax.grid(True, color='#2C3E50', alpha=0.5)
    for spine in ax.spines.values():
        spine.set_edgecolor('#2C3E50')
    plt.tight_layout()
    plt.savefig(outpath, dpi=130, bbox_inches='tight', facecolor='#1A1A2E')
    print(f"Parâmetro de ordem salvo em: {outpath}")


# ══════════════════════════════════════════════════════════════════════════════
# 10. INTERFACE DE LINHA DE COMANDO (main)
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """
    Ponto de entrada da simulação via linha de comando.

    Fluxo de execução
    -----------------
    1. Parseia os argumentos (geometria, espaçamento, campo externo, etc.)
    2. Converte (B_ext, phi_ext) → (Bext_x, Bext_y) em coordenadas cartesianas
    3. Gera a grade de posições e ângulos iniciais (make_grid)
    4. Calcula o cutoff e tamanho das agulhas de acordo com a geometria
    5. Salva figura do estado inicial
    6. Executa a relaxação (relax)
    7. Salva figura do estado de equilíbrio
    8. Salva figura comparativa inicial vs. equilíbrio
    9. Plota o parâmetro de ordem ao longo da relaxação
    10. (Opcional) Gera animação da relaxação
    """
    parser = argparse.ArgumentParser(
        description="Simulacao de rede de bussolas — campo dipolar 2D",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # ── grupo 1: geometria da rede ────────────────────────────────────────
    grp = parser.add_argument_group("Geometria da rede")
    grp.add_argument('--N', type=int, default=8,
                     help='Numero de LINHAS de agulhas')
    grp.add_argument('--M', type=int, default=8,
                     help='Numero de COLUNAS de agulhas')
    grp.add_argument('--R', type=float, default=0.025,
                     help='Raio do circulo que envolve cada agulha [m]. '
                          'A distancia entre centros de vizinhos = 2R. '
                          'Padrao: 0.025 m = 2.5 cm')
    grp.add_argument('--needle_frac', type=float, default=0.80,
                     help='Comprimento da agulha como fracao do diametro 2R '
                          '(0.0–1.0). Padrao: 0.80  →  agulha = 0.80 * 2R')
    grp.add_argument('--geometry',
                     choices=['square', 'triangular', 'honeycomb'],
                     default='square',
                     help='Tipo de rede')

    # ── grupo 2: parâmetros físicos e de simulação ────────────────────────
    grp2 = parser.add_argument_group("Fisica e simulacao")
    grp2.add_argument('--moment', type=float, default=MOMENT_DEFAULT,
                      help='Momento magnetico de cada agulha [A·m²]. '
                           'Tipico: 0.01 (pequena) a 1.0 (grande). '
                           f'Padrao: {MOMENT_DEFAULT} A·m²')
    grp2.add_argument('--inertia', type=float, default=INERTIA_DEFAULT,
                      help='Momento de inercia de cada agulha [kg·m²]. '
                           'Barra fina: I = (1/12)·massa·comprimento². '
                           f'Padrao: {INERTIA_DEFAULT:.2e} kg·m² '
                           '(agulha 0.5 g × 5 cm)')
    grp2.add_argument('--damping', type=float, default=DAMPING_DEFAULT,
                      help='Amortecimento viscoso do ar b [N·m·s/rad]. '
                           'Zero = sem amortecimento (oscila infinitamente). '
                           'Grande = sem oscilacoes visíveis. '
                           f'Padrao: {DAMPING_DEFAULT:.2e} N·m·s/rad '
                           '(sub-amortecido, Q >> 1)')
    grp2.add_argument('--t_sim', type=float, default=2.0,
                      help='Tempo TOTAL de simulacao [s]. O numero de passos '
                           'e calculado automaticamente a partir de dt. '
                           'A simulacao pode terminar antes se convergir.')
    grp2.add_argument('--dt_factor', type=float, default=0.05,
                      help='Fracao do periodo natural T0 usada como passo dt '
                           '(0.02-0.10; menor = mais preciso, mais lento)')
    grp2.add_argument('--noise', type=float, default=1.5,
                      help='Amplitude do ruido inicial nos angulos [rad]; '
                           '0=todas para +x, 3.14=aleatorio total')
    grp2.add_argument('--seed', type=int, default=42,
                      help='Semente do gerador aleatorio (reprodutibilidade)')

    # ── grupo 3: campo externo ────────────────────────────────────────────
    grp3 = parser.add_argument_group(
        "Campo externo uniforme (SI)",
        "Intensidade em Tesla. Exemplos: campo terrestre = 50e-6 T; "
        "ima de geladeira a 5 cm ≈ 1e-3 T.")
    grp3.add_argument('--B_ext', type=float, default=0.0,
                      help='Intensidade do campo externo [T]. '
                           'Ex: 50e-6 (campo terrestre), 1e-3 (1 mT)')
    grp3.add_argument('--phi_ext', type=float, default=0.0,
                      help='Direcao do campo externo [graus]: '
                           '0=direita (+x), 90=cima (+y), anti-horario')
    grp3.add_argument('--ext_Bx', type=float, default=None,
                      help='Componente Bx do campo externo [T] '
                           '(sobrescreve --B_ext/--phi_ext se informado)')
    grp3.add_argument('--ext_By', type=float, default=None,
                      help='Componente By do campo externo [T] '
                           '(sobrescreve --B_ext/--phi_ext se informado)')

    # ── grupo 4: saída ────────────────────────────────────────────────────
    grp4 = parser.add_argument_group("Saida")
    grp4.add_argument('--video', type=str, default=None,
                      metavar='ARQUIVO.mp4',
                      help='Gera vídeo MP4 da relaxação e salva neste caminho '
                           '(requer ffmpeg instalado)')
    grp4.add_argument('--frame_every', type=int, default=5,
                      help='Salva um frame a cada N passos (menor = mais suave, '
                           'mais lento)')
    grp4.add_argument('--fps', type=int, default=24,
                      help='Quadros por segundo do vídeo MP4')
    grp4.add_argument('--keep_frames', action='store_true',
                      help='Mantém a pasta de PNGs intermediários após gerar o MP4')
    args = parser.parse_args()

    np.random.seed(args.seed)

    # ── converte campo externo para coordenadas cartesianas ───────────────
    # Prioridade: --ext_Bx/--ext_By (cartesiano) sobrescreve --B_ext/--phi_ext
    if args.ext_Bx is not None or args.ext_By is not None:
        Bext_x = args.ext_Bx if args.ext_Bx is not None else 0.0
        Bext_y = args.ext_By if args.ext_By is not None else 0.0
        B_ext_mag = np.sqrt(Bext_x**2 + Bext_y**2)
        phi_ext_deg = np.degrees(np.arctan2(Bext_y, Bext_x))
    else:
        # converte (intensidade, ângulo) → (Bx, By)
        phi_rad = np.deg2rad(args.phi_ext)
        Bext_x  = args.B_ext * np.cos(phi_rad)
        Bext_y  = args.B_ext * np.sin(phi_rad)
        B_ext_mag   = args.B_ext
        phi_ext_deg = args.phi_ext

    # ── resumo dos parâmetros ─────────────────────────────────────────────
    # formata campo para exibição em unidade mais legível
    def fmt_field(B):
        if B == 0:      return "0 T"
        if B >= 0.1:    return f"{B:.4f} T"
        if B >= 1e-4:   return f"{B*1e3:.4f} mT"
        return              f"{B*1e6:.2f} µT"

    print(f"\n{'═'*62}")
    print(f"  Rede         : {args.geometry}  {args.N}×{args.M} agulhas")
    print(f"  Raio R       : {args.R*100:.2f} cm  "
          f"(dist. entre vizinhos = 2R = {2*args.R*100:.2f} cm)")
    print(f"  Agulha       : {args.needle_frac*100:.0f}% de 2R "
          f"= {args.needle_frac*2*args.R*100:.2f} cm de comprimento")
    print(f"  Momento mag. : {args.moment:.4g} A·m²  por agulha")
    print(f"  Inércia      : {args.inertia:.3e} kg·m²  por agulha")
    print(f"  Amortecimento: {args.damping:.3e} N·m·s/rad  (ar)")
    print(f"  Campo externo: {fmt_field(B_ext_mag)}  φ={phi_ext_deg:.1f}°"
          f"  → Bx={fmt_field(Bext_x)}  By={fmt_field(Bext_y)}")
    print(f"  Simulação    : t_sim={args.t_sim:.3f} s  dt_factor={args.dt_factor}"
          f"  ruído={args.noise:.2f} rad  seed={args.seed}")
    print(f"{'═'*62}\n")

    # ── gera grade de posições e ângulos iniciais ─────────────────────────
    # make_grid usa R como parâmetro único; retorna também nn_dist = 2R
    xs, ys, thetas_init, nn_dist = make_grid(
        args.N, args.M,
        geometry=args.geometry,
        noise=args.noise,
        R=args.R,
    )

    # ── tamanhos derivados de R ───────────────────────────────────────────
    # R  = raio do círculo visual de cada agulha
    # 2R = distância entre vizinhos (nn_dist) em qualquer geometria
    # needle_len = fração de 2R controlada por --needle_frac
    # needle_width = 22% de needle_len (proporção visual de bússola)
    # r_halo = R (exatamente o raio do círculo) — círculos se tocam tangencialmente
    # cutoff = 2.6 * 2R — cobre 1ª e 2ª camadas de vizinhos
    R          = args.R
    needle_len = args.needle_frac * 2.0 * R   # comprimento da agulha [m]
    needle_width = needle_len * 0.22           # largura do losango [m]
    r_halo     = R * 0.98                      # ligeiramente menor que R para gap visível
    cutoff     = nn_dist * 2.6                 # raio de corte da interação dipolar [m]

    # campo dipolar de referência entre primeiros vizinhos (para exibição)
    B_ref = MU0_OVER_4PI * 2.0 * args.moment / nn_dist**3
    print(f"  Campo dipolar de referência (vizinho mais próximo): {fmt_field(B_ref)}")
    if B_ext_mag > 0:
        ratio = B_ext_mag / B_ref
        print(f"  B_ext / B_ref = {ratio:.3f}  "
              + ("(campo externo DOMINANTE)" if ratio > 1
                 else "(campo externo fraco — interações dominam)"))
    print()

    ext_kwargs = dict(B_ext=B_ext_mag, phi_ext_deg=phi_ext_deg,
                      r_halo=r_halo)

    # ── figura do estado inicial ──────────────────────────────────────────
    fig0, _ = plot_state(thetas_init, xs, ys,
                         title="Estado inicial (aleatório)",
                         needle_len=needle_len, needle_width=needle_width,
                         **ext_kwargs)
    plt.tight_layout()
    plt.savefig("compass_initial.png",
                dpi=130, bbox_inches='tight', facecolor='#1A1A2E')
    plt.close(fig0)
    print("Estado inicial salvo.")

    # ── integração dinâmica inercial ──────────────────────────────────────
    frame_dir   = None
    final_video = None
    if args.video:
        import os
        # resolve o nome final do vídeo ANTES de criar a pasta de frames,
        # garantindo que pasta e arquivo usem o mesmo nome base
        final_video = next_available_path(args.video)
        if final_video != args.video:
            print(f"  '{args.video}' já existe → será salvo como '{final_video}'")
        base      = os.path.splitext(final_video)[0]
        frame_dir = base + "_frames"
        print(f"Integrando e gravando frames em '{frame_dir}/'...")
    else:
        print("Integrando dinâmica inercial...")

    thetas_eq, omegas_eq, thetas_hist, n_frames, sim_dt, stop_reason = relax(
        thetas_init.copy(), xs, ys,
        t_sim=args.t_sim,
        dt_factor=args.dt_factor,
        inertia=args.inertia,
        damping=args.damping,
        cutoff=cutoff,
        ext_field=(Bext_x, Bext_y),
        moment=args.moment,
        frame_dir=frame_dir,
        frame_every=args.frame_every,
        needle_len=needle_len,
        needle_width=needle_width,
        r_halo=r_halo,
        B_ext=B_ext_mag,
        phi_ext_deg=phi_ext_deg,
    )
    print(f"Integração concluída — {stop_reason}"
          + (f"  ({n_frames} frames salvos)" if n_frames else ""))

    # ── figura do estado de equilíbrio ────────────────────────────────────
    title_eq = (f"Equilíbrio — {args.geometry} {args.N}×{args.M}"
                + (f"  |  B={fmt_field(B_ext_mag)} @ {phi_ext_deg:.0f}°"
                   if B_ext_mag > 0 else ""))
    fig1, _ = plot_state(thetas_eq, xs, ys,
                         title=title_eq,
                         needle_len=needle_len, needle_width=needle_width,
                         **ext_kwargs)
    plt.tight_layout()
    plt.savefig("compass_equilibrium.png",
                dpi=130, bbox_inches='tight', facecolor='#1A1A2E')
    plt.close(fig1)
    print("Estado de equilíbrio salvo.")

    # ── figura comparativa lado a lado ────────────────────────────────────
    fig2, axes = plt.subplots(1, 2, figsize=(15, 7), facecolor='#1A1A2E')
    margin = needle_len * 1.6
    for ax_ in axes:
        ax_.set_facecolor('#16213E')
        ax_.set_aspect('equal')
        ax_.set_xlim(xs.min() - margin, xs.max() + margin)
        ax_.set_ylim(ys.min() - margin, ys.max() + margin * 2.5)
        ax_.tick_params(left=False, bottom=False,
                        labelleft=False, labelbottom=False)
        for sp in ax_.spines.values():
            sp.set_edgecolor('#2C3E50')

    axes[0].set_title("Estado inicial", color='#ECF0F1',
                      fontsize=12, fontfamily='monospace')
    axes[1].set_title("Equilíbrio dipolar", color='#ECF0F1',
                      fontsize=12, fontfamily='monospace')

    N_g, M_g = thetas_init.shape
    for i in range(N_g):
        for j in range(M_g):
            draw_compass(axes[0], xs[i, j], ys[i, j], thetas_init[i, j],
                         length=needle_len, width=needle_width)
            draw_compass(axes[1], xs[i, j], ys[i, j], thetas_eq[i, j],
                         length=needle_len, width=needle_width)

    # setas de campo externo em ambos os painéis
    draw_ext_field_on_lattice(axes[0], xs, ys, B_ext_mag, phi_ext_deg, needle_len)
    draw_ext_field_on_lattice(axes[1], xs, ys, B_ext_mag, phi_ext_deg, needle_len)

    bfield_str = (f"  |  B_ext={B_ext_mag:.2f} @ {phi_ext_deg:.0f}°"
                  if B_ext_mag > 0 else "")
    fig2.suptitle(
        f"Rede {args.geometry} {args.N}×{args.M}"
        f"  |  R={args.R*100:.1f} cm  2R={2*args.R*100:.1f} cm"
        f"{bfield_str}  |  interação dipolar 2D",
        color='#BDC3C7', fontsize=11, fontfamily='monospace')
    plt.tight_layout()
    plt.savefig("compass_comparison.png",
                dpi=130, bbox_inches='tight', facecolor='#1A1A2E')
    plt.close(fig2)
    print("Comparação salva.")

    # ── parâmetro de ordem global ─────────────────────────────────────────
    # hist contém tuplas (thetas, omegas); extraímos só os thetas
    plot_order_parameter([th for th, _ in thetas_hist],
                         "compass_order_param.png",
                         dt=sim_dt)

    # ── geração do vídeo MP4 ──────────────────────────────────────────────
    if final_video and frame_dir and n_frames > 0:
        import shutil
        ok = render_video(frame_dir, final_video, fps=args.fps)
        if ok and not args.keep_frames:
            shutil.rmtree(frame_dir)
            print(f"  Pasta de frames removida: {frame_dir}/")
        elif not ok:
            print(f"  Os frames PNG foram mantidos em: {frame_dir}/")

    print("\nConcluído.")


if __name__ == '__main__':
    main()
