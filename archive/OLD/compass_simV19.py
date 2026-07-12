"""
============================================================
compass_sim.py — Simulação de rede de agulhas de bússola
============================================================

Modela uma grade 2D de dipolos magnéticos clássicos (agulhas de
bússola) que interagem pelo campo que cada uma gera sobre as vizinhas.
A dinâmica é inercial (2ª Lei de Newton rotacional) sem atrito no pino,
com amortecimento viscoso do ar. Integrador Velocity-Verlet.
Todas as grandezas físicas em unidades SI.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PARÂMETROS DE LINHA DE COMANDO
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

┌─ GEOMETRIA DA REDE ────────────────────────────────────────┐
│                                                            │
│  --geometry  square | triangular | honeycomb               │
│              Tipo de rede. Padrão: square                  │
│              · square     : grade retangular               │
│              · triangular : grade triangular equilateral   │
│              · honeycomb  : colmeia (buracos hexagonais)   │
│                                                            │
│  --N         int   Número de LINHAS de agulhas. Padrão: 8  │
│  --M         int   Número de COLUNAS de agulhas. Padrão: 8 │
│                                                            │
│  --R         float Raio do círculo de cada agulha [m].     │
│                    Distância entre vizinhos = 2R.          │
│                    Padrão: 0.025 m (2.5 cm)                │
│                                                            │
│  --needle_frac float Comprimento da agulha como fração     │
│                    do diâmetro 2R (0.0–1.0).               │
│                    Padrão: 0.80  → agulha = 80% de 2R     │
└────────────────────────────────────────────────────────────┘

┌─ FÍSICA E SIMULAÇÃO ───────────────────────────────────────┐
│                                                            │
│  --moment    float Momento magnético de cada agulha [A·m²] │
│                    Padrão: 0.1  (bússola de mesa, ~5 cm)   │
│                    Ref: bolso ≈ 0.01 | náutica ≈ 1.0       │
│                                                            │
│  --inertia   float Momento de inércia [kg·m²].             │
│                    Barra fina: I = (1/12)·massa·L²         │
│                    Padrão: 1e-7  (0.5 g × 5 cm)           │
│                                                            │
│  --damping   float Amortecimento viscoso do ar [N·m·s/rad] │
│                    Controla o fator de qualidade Q:        │
│                    Q = omega_0·I/b  (Q alto → mais oscilações) │
│                    Padrão: 5e-8 (Q≈25, bússola realista)   │
│                    Para B_ext=0.1T suave: use 8e-6 (Q≈4)  │
│                                                            │
│  --t_sim     float Tempo físico total de simulação [s].     │
│                    É a soma de todos os passos dt integrados │
│                    — equivalente ao que um cronômetro real   │
│                    marcaria observando as agulhas se mover.  │
│                    O vídeo exibe esse tempo físico.          │
│                    A sim. para antes se S=1.00 ou repouso.  │
│                    Padrão: 2.0 s                            │
│                                                            │
│  --dt_factor float Fração do período natural T₀ usada      │
│                    como passo de integração (0.02–0.10).   │
│                    Menor = mais preciso, mais lento.        │
│                    Padrão: 0.05                            │
│                                                            │
│  --noise     float Amplitude do ruído inicial [rad].       │
│                    0 = todas apontam para +x               │
│                    π ≈ 3.14 = orientação totalmente aleat. │
│                    Padrão: 1.5                             │
│                                                            │
│  --seed      int   Semente do gerador aleatório.           │
│                    Garante reprodutibilidade dos resultados.│
│                    Padrão: 42                              │
└────────────────────────────────────────────────────────────┘

┌─ CAMPO EXTERNO UNIFORME (unidades SI) ─────────────────────┐
│  Duas formas de especificar — não usar as duas juntas.     │
│                                                            │
│  Forma A — intensidade + ângulo (recomendada):             │
│  --B_ext     float Intensidade do campo [T].               │
│                    0.0      = sem campo (padrão)           │
│                    50e-6    = campo terrestre (≈50 µT)     │
│                    1e-3     = ímã de geladeira a 5 cm      │
│                    0.1      = campo forte (alinha tudo)    │
│                                                            │
│  --phi_ext   float Direção do campo [graus].               │
│                    0   = direita (+x)  ← padrão            │
│                    90  = cima (+y)                         │
│                    180 = esquerda (−x)                     │
│                    270 = baixo (−y)                        │
│                    Sentido anti-horário. Padrão: 0.0       │
│                                                            │
│  Forma B — componentes cartesianas (sobrescreve A):        │
│  --ext_Bx    float Componente x do campo [T]               │
│  --ext_By    float Componente y do campo [T]               │
└────────────────────────────────────────────────────────────┘

┌─ SAÍDA ─────────────────────────────────────────────────────┐
│  Arquivos PNG sempre gerados no diretório atual:           │
│    compass_initial.png      estado inicial                  │
│    compass_equilibrium.png  estado final                    │
│    compass_comparison.png   comparação lado a lado          │
│    compass_order_param.png  parâmetro de ordem S(t)         │
│                                                            │
│  --video     str   Caminho do vídeo MP4 a gerar.           │
│                    Requer ffmpeg instalado.                 │
│                    Se o arquivo já existir, salva como     │
│                    nome0001.mp4, nome0002.mp4, etc.        │
│                    Ex: --video simulacao.mp4               │
│                                                            │
│  --frame_every int Salva um frame a cada N passos.         │
│                    Menor = vídeo mais suave, mais lento.   │
│                    Padrão: 5                               │
│                                                            │
│  --fps       int   Quadros por segundo do vídeo MP4.       │
│                    Padrão: 24                              │
│                                                            │
│  --keep_frames     Se presente, mantém a pasta de PNGs     │
│                    intermediários após gerar o MP4.        │
└────────────────────────────────────────────────────────────┘

┌─ CONTROLES DURANTE A SIMULAÇÃO ────────────────────────────┐
│  A simulação para automaticamente quando:                  │
│    · S = 1.00  (todas as agulhas alinhadas)                │
│    · Rede em repouso (ω_max → 0)                           │
│    · Tempo t_sim atingido                                  │
│    · Ctrl+I (Tab) pressionado no terminal                  │
│      → interrompe e salva o vídeo imediatamente            │
└────────────────────────────────────────────────────────────┘

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
EXEMPLOS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  # Padrões — rede quadrada 8×8 sem campo
  python compass_sim.py

  # Honeycomb com campo terrestre
  python compass_sim.py --geometry honeycomb --N 10 --M 10 --B_ext 50e-6

  # Triangular com campo de 0.1 T a 45°, movimento suave
  python compass_sim.py --geometry triangular --N 10 --M 10 \\
      --B_ext 0.1 --phi_ext 45 --damping 8e-6 --t_sim 2.0

  # Vídeo com agulhas maiores e muitas oscilações
  python compass_sim.py --R 0.03 --needle_frac 0.85 --damping 1e-9 \\
      --t_sim 5.0 --frame_every 2 --fps 30 --video sim.mp4

  # Campo via componentes cartesianas
  python compass_sim.py --ext_Bx 0.05 --ext_By -0.05

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Dependências: numpy, matplotlib, ffmpeg (para vídeo)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.animation import FuncAnimation
import argparse
import time
import sys

# ── configura saída para flush imediato e line endings corretos ──────────────
# No macOS com zsh, o modo raw do teclado (usado pelo Ctrl+I) pode desabilitar
# o processamento de \n → \r\n no terminal. Forçar line_buffering e flush
# garante que cada print() apareça em sua própria linha corretamente.
sys.stdout.reconfigure(line_buffering=True)

def _print(*args, **kwargs):
    """
    Substituto de print() que garante flush imediato e \r\n correto
    mesmo quando o terminal está em modo raw (macOS/zsh com Ctrl+I ativo).
    """
    import sys
    msg = ' '.join(str(a) for a in args)
    sys.stdout.write(msg + '\r\n')
    sys.stdout.flush()

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
#   Q = omega_0·I/b   onde  omega_0 = sqrt(m·B_eff/I)
#
# Com B_ext = 0.1 T:   omega_0 ≈ 316 rad/s  →  b = 8e-6 dá Q ≈ 4  (suave)
# Com B_ext = 0  :     omega_0 ≈  13 rad/s  →  b = 5e-8 dá Q ≈ 25 (bússola)
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

def total_field_on(i, j, thetas, xs, ys, cutoff, moment,
                   pbc=False, Lx=None, Ly=None, n_images=1):
    """
    Soma o campo dipolar de todas as agulhas vizinhas sobre a agulha (i, j).

    Para eficiência computacional, apenas agulhas dentro do raio `cutoff`
    são consideradas. Isso é fisicamente justificado porque o campo dipolar
    decai como 1/r³ — contribuições de agulhas distantes são desprezíveis.

    Condições periódicas de contorno (PBC)
    ---------------------------------------
    Quando pbc=True, a rede é tratada como uma célula unitária que se repete
    infinitamente em x e em y (estrutura periódica, sem bordas). Para cada
    par de agulhas (i,j)-(ni,nj), o campo é somado sobre TODAS as réplicas
    periódicas dentro de ±n_images células em x e em y — não apenas a mais
    próxima. Com n_images=1 (padrão), cada vizinha contribui através de
    (2·1+1)² = 9 réplicas (a célula original + 1 réplica de cada lado
    em x e em y — equivalente à convenção de imagem mínima clássica).
    Aumentar n_images soma mais réplicas distantes, útil quando o cutoff
    é grande em relação ao período da rede.

    Isso elimina os efeitos de borda: agulhas na margem da rede passam a
    "ver" vizinhas do lado oposto (e suas réplicas mais distantes), como
    se a estrutura fosse infinita.

    Parâmetros
    ----------
    i, j     : índices da agulha receptora na grade (N×M)
    thetas   : array 2D dos ângulos atuais de todas as agulhas  [rad]
    xs, ys   : arrays 2D das posições fixas das agulhas  [m]
    cutoff   : raio máximo de interação  [m]
    moment   : momento magnético de cada agulha  [A·m²]
    pbc      : se True, aplica condições periódicas de contorno  [bool]
    Lx, Ly   : dimensões do período da rede em x e y  [m]
               (necessários apenas se pbc=True; calculados em make_grid)
    n_images : número de réplicas periódicas a somar de cada lado, em
               cada direção (padrão: 1). Total de réplicas consideradas
               por vizinha = (2·n_images+1)² células.

    Retorna
    -------
    (Bx_tot, By_tot) : campo total na posição da agulha (i, j)  [T]
    """
    Bx_tot, By_tot = 0.0, 0.0
    N, M = thetas.shape
    xi, yi = xs[i, j], ys[i, j]

    # deslocamentos das réplicas periódicas a considerar (0 = célula original)
    if pbc and Lx and Ly:
        img_range = range(-n_images, n_images + 1)
        x_shifts = [k * Lx for k in img_range]
        y_shifts = [k * Ly for k in img_range]
    else:
        x_shifts = [0.0]
        y_shifts = [0.0]

    for ni in range(N):
        for nj in range(M):
            if ni == i and nj == j and not pbc:
                continue  # sem PBC: a agulha não interage com ela mesma

            rx0 = xi - xs[ni, nj]   # [m]
            ry0 = yi - ys[ni, nj]   # [m]

            for dx_img in x_shifts:
                for dy_img in y_shifts:
                    # pula a própria agulha na célula original (distância zero)
                    if ni == i and nj == j and dx_img == 0.0 and dy_img == 0.0:
                        continue

                    rx = rx0 + dx_img
                    ry = ry0 + dy_img
                    dist = np.sqrt(rx*rx + ry*ry)

                    if dist > cutoff:
                        continue  # fora do raio de corte: ignora

                    bx, by = dipole_field_2d(rx, ry, thetas[ni, nj], moment)
                    Bx_tot += bx
                    By_tot += by

    return Bx_tot, By_tot   # [T]

    return Bx_tot, By_tot   # [T]


# ══════════════════════════════════════════════════════════════════════════════
# 4. DINÂMICA INERCIAL (2ª Lei de Newton para rotação — sem atrito no pino)
# ══════════════════════════════════════════════════════════════════════════════

def _plot_hysteresis(log):
    """
    Plota a curva de histerese M_proj(B) e salva como PNG e CSV.

    O gráfico mostra a magnetização projetada na direção do campo em função
    da intensidade do campo aplicado — a clássica curva de histerese M×H.

    Parâmetros
    ----------
    log : lista de tuplas (t, B_scalar, M_proj, S)
          gerada durante a integração com field_mode='hysteresis'
    """
    log = np.array(log)
    t_arr  = log[:, 0]
    B_arr  = log[:, 1] * 1e3    # converte T → mT para legibilidade
    M_arr  = log[:, 2]          # magnetização projetada (adimensional, ∈ [−1,1])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), facecolor='#1A1A2E')

    # ── painel esquerdo: M×B (curva de histerese) ─────────────────────────
    ax1 = axes[0]
    ax1.set_facecolor('#0D1B2A')
    ax1.plot(B_arr, M_arr, color='#E94560', lw=1.2, alpha=0.85)
    ax1.axhline(0, color='#4A4A6A', lw=0.6, ls='--')
    ax1.axvline(0, color='#4A4A6A', lw=0.6, ls='--')
    ax1.set_xlabel("Campo B  [mT]",    color='#BDC3C7')
    ax1.set_ylabel("Magnetização  M (projeção)",  color='#BDC3C7')
    ax1.set_title("Curva de Histerese  M × B", color='#ECF0F1',
                  fontfamily='monospace')
    ax1.tick_params(colors='#7F8C8D')
    ax1.grid(True, color='#2C3E50', alpha=0.5)
    for sp in ax1.spines.values():
        sp.set_edgecolor('#2C3E50')

    # ── painel direito: M(t) e B(t) versus tempo ──────────────────────────
    ax2 = axes[1]
    ax2.set_facecolor('#0D1B2A')
    l1, = ax2.plot(t_arr, M_arr, color='#E94560', lw=1.2, label='M (proj.)')
    ax2b = ax2.twinx()
    ax2b.plot(t_arr, B_arr, color='#FFD700', lw=1.0, alpha=0.7,
              ls='--', label='B (mT)')
    ax2b.set_ylabel("Campo B  [mT]", color='#FFD700')
    ax2b.tick_params(colors='#FFD700')
    ax2.set_xlabel("Tempo físico  t  [s]", color='#BDC3C7')
    ax2.set_ylabel("Magnetização  M", color='#E94560')
    ax2.set_title("Evolução temporal", color='#ECF0F1', fontfamily='monospace')
    ax2.tick_params(colors='#7F8C8D')
    ax2.grid(True, color='#2C3E50', alpha=0.4)
    for sp in ax2.spines.values():
        sp.set_edgecolor('#2C3E50')
    lines = [l1, plt.Line2D([0],[0], color='#FFD700', ls='--')]
    ax2.legend(lines, ['M (proj.)', 'B (mT)'], facecolor='#1A1A2E',
               edgecolor='#2C3E50', labelcolor='#ECF0F1', fontsize=8)

    fig.suptitle("Simulação de Histerese Magnética — Rede de Bússolas",
                 color='#BDC3C7', fontsize=12, fontfamily='monospace')
    plt.tight_layout()
    plt.savefig("hysteresis_loop.png", dpi=130, bbox_inches='tight',
                facecolor='#1A1A2E')
    plt.close(fig)
    _print("  Gráfico de histerese salvo: hysteresis_loop.png")


def _plot_sine(log, freq):
    """
    Plota M(t) e B(t) para o modo de campo senoidal e salva como PNG.

    Parâmetros
    ----------
    log  : lista de tuplas (t, B_scalar, M_proj, S)
    freq : frequência do campo senoidal [Hz]
    """
    log   = np.array(log)
    t_arr = log[:, 0]
    B_arr = log[:, 1] * 1e3    # T → mT
    M_arr = log[:, 2]

    fig, ax = plt.subplots(figsize=(10, 4), facecolor='#1A1A2E')
    ax.set_facecolor('#0D1B2A')
    ax.plot(t_arr, M_arr, color='#E94560', lw=1.2, label='M (proj.)')
    ax2 = ax.twinx()
    ax2.plot(t_arr, B_arr, color='#FFD700', lw=1.0, alpha=0.7,
             ls='--', label=f'B (mT)  f={freq:.2f} Hz')
    ax.set_xlabel("Tempo físico  t  [s]", color='#BDC3C7')
    ax.set_ylabel("Magnetização  M", color='#E94560')
    ax2.set_ylabel("Campo B  [mT]", color='#FFD700')
    ax.set_title(f"Campo Senoidal — f = {freq:.3f} Hz",
                 color='#ECF0F1', fontfamily='monospace')
    ax.tick_params(colors='#7F8C8D')
    ax2.tick_params(colors='#FFD700')
    ax.grid(True, color='#2C3E50', alpha=0.4)
    for sp in ax.spines.values():
        sp.set_edgecolor('#2C3E50')
    lines = [plt.Line2D([0],[0], color='#E94560'),
             plt.Line2D([0],[0], color='#FFD700', ls='--')]
    ax.legend(lines, ['M (proj.)', f'B (mT)  f={freq:.2f}Hz'],
              facecolor='#1A1A2E', edgecolor='#2C3E50',
              labelcolor='#ECF0F1', fontsize=8)
    plt.tight_layout()
    plt.savefig("sine_field.png", dpi=130, bbox_inches='tight',
                facecolor='#1A1A2E')
    plt.close(fig)
    _print("  Gráfico senoidal salvo: sine_field.png")


def relax(thetas, xs, ys, t_sim=2.0, dt_factor=0.05,
          inertia=INERTIA_DEFAULT, damping=DAMPING_DEFAULT,
          cutoff=3.5, ext_field=(0.0, 0.0), moment=MOMENT_DEFAULT,
          field_mode='static', field_freq=1.0,
          callback=None,
          frame_dir=None, frame_every=10,
          needle_len=0.042, needle_width=0.010, r_halo=None,
          frame_dpi=120, figsize_inches=None,
          pbc=False, Lx=None, Ly=None, n_images=1,
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

    Modos de campo externo (field_mode)
    ------------------------------------
    'static'
        Campo constante (Bext_x, Bext_y) em todo o intervalo.

    'hysteresis'
        Campo varia linearmente no tempo ao longo de 3 rampas:
          0 → t_sim/3  : B sobe de 0 até +B_max
          t_sim/3 → 2t_sim/3 : B desce de +B_max até −B_max
          2t_sim/3 → t_sim   : B sobe de −B_max até +B_max
        A direção é sempre a de phi_ext. A amplitude B_max é |ext_field|.
        Desabilita a parada por S=1.00 (a magnetização oscila).

    'sine'
        Campo senoidal: B(t) = B_max · sin(2π·f·t)
        Direção: phi_ext. Amplitude: B_max = |ext_field|. Freq: field_freq [Hz].
        Desabilita a parada por S=1.00.

    Parâmetros
    ----------
    thetas      : array N×M de ângulos iniciais  [rad]
    xs, ys      : posições fixas das agulhas  [m]
    t_sim       : tempo físico total da simulação  [s]
                  (= soma dos passos dt, não tempo de CPU)
    dt_factor   : fração do período natural T₀ usada como passo de tempo
    inertia     : momento de inércia de cada agulha  [kg·m²]
    damping     : coeficiente de amortecimento viscoso b  [N·m·s/rad]
    cutoff      : raio máximo de interação dipolar  [m]
    ext_field   : tupla (Bext_x, Bext_y) — campo base [T]
                  Em 'static': campo constante.
                  Em 'hysteresis'/'sine': define direção e amplitude máxima.
    moment      : momento magnético de cada agulha  [A·m²]
    field_mode  : 'static' | 'hysteresis' | 'sine'
    field_freq  : frequência do campo senoidal [Hz] (só para field_mode='sine')
    callback    : função opcional callback(step, thetas, omegas)
    frame_dir   : diretório onde salvar PNGs de frames (ou None)
    frame_every : intervalo em passos entre frames salvos
    needle_len  : tamanho das agulhas para renderização  [m]
    needle_width: largura das agulhas  [m]
    r_halo      : raio dos halos de parâmetro de ordem  [m]
    pbc         : se True, aplica condições periódicas de contorno em x e y
                  (a rede se comporta como uma estrutura infinita periódica;
                  agulhas nas bordas interagem com réplicas do lado oposto)
    Lx, Ly      : dimensões do período da rede em x e y  [m]
                  (necessários apenas se pbc=True; calculados em make_grid)
    n_images    : número de réplicas periódicas a somar de cada lado em
                  cada direção quando pbc=True (padrão: 1 → soma sobre
                  uma grade de réplicas (2·1+1)² = 9 células)
    B_ext       : intensidade do campo externo para rótulo nos frames  [T]
    phi_ext_deg : direção do campo externo  [graus]

    Retorna
    -------
    theta_cur   : array N×M — ângulos no estado final  [rad]
    omega_cur   : array N×M — velocidades angulares  [rad/s]
    hist        : lista de tuplas (thetas, omegas) a cada 20 passos
    n_frames    : número de frames PNG salvos
    dt          : passo de tempo usado  [s]
    stop_reason : string descrevendo como a simulação terminou
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

    _print(f"  Dinamica inercial:")
    _print(f"    r_nn    = {r_nn*100:.2f} cm")
    _print(f"    B_ref   = {B_ref*1e3:.4f} mT  (dipolar entre vizinhos)")
    if B_ext_mag > 0:
        _print(f"    B_ext   = {B_ext_mag*1e3:.4f} mT  (campo externo)")
    _print(f"    B_eff   = {B_eff*1e3:.4f} mT  (campo dominante)")
    _print(f"    omega_0 = {omega0:.2f} rad/s   T0 = {T0:.5f} s")
    _print(f"    dt      = {dt:.6f} s  ({dt_factor:.0%} de T0)")
    _print(f"    t_sim   = {t_sim:.3f} s  -> {n_steps} passos")
    if Q > 2:
        q_desc = "sub-amortecido (oscila)"
    elif Q > 0.5:
        q_desc = "criticamente amortecido"
    else:
        q_desc = "super-amortecido"
    _print(f"    Q       = {Q:.1f}  ({q_desc})")
    _print()

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

    def _draw_clock(ax, t_phys, needle_len, stop_label=None):
        """
        Desenha um cronômetro no canto superior esquerdo mostrando o
        TEMPO FÍSICO DO SISTEMA — a soma dos passos de integração dt,
        equivalente ao tempo que um cronômetro real marcaria observando
        as agulhas se mover.

        t_phys = step * dt   [s]

        Elementos:
          - Texto grande: t = X.XXXX s  (tempo físico atual)
          - Subtexto: / X.XXXX s        (tempo físico total t_sim)
          - Barra de progresso: t_phys / t_sim
          - Rótulo de parada quando a simulação termina antes de t_sim
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
        bar_y = py - 0.095 * yspan

        extra_h = 0.030 * yspan if stop_label else 0.0

        # caixa de fundo semitransparente
        pad = needle_len * 0.25
        ax.add_patch(plt.Rectangle(
            (px - pad * 0.3, bar_y - pad * 0.6 - extra_h),
            bar_w + pad, 0.140 * yspan + pad + extra_h,
            facecolor='#080818', edgecolor='#3A3A6A',
            linewidth=0.8, alpha=0.82, zorder=19,
            transform=ax.transData))

        # ── barra de progresso: t_phys / t_sim ────────────────────────────
        frac = min(t_phys / t_sim, 1.0) if t_sim > 0 else 0.0

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

        # marcador branco na posição atual
        if 0 < frac < 1.0:
            mx = px + bar_w * frac
            ax.plot([mx, mx], [bar_y, bar_y + bar_h],
                    color='white', lw=1.2, zorder=24,
                    transform=ax.transData)

        # ── texto: tempo físico atual ──────────────────────────────────────
        # formata em mm:ss.ss se >= 60 s, senão em s com 4 casas decimais
        if t_phys >= 60.0:
            mins = int(t_phys // 60)
            secs = t_phys - mins * 60
            time_str = f"t = {mins:02d}:{secs:05.2f}"
        else:
            time_str = f"t = {t_phys:.4f} s"

        ax.text(px + pad * 0.2, py,
                time_str,
                color='#E8E8FF', fontsize=11, fontweight='bold',
                fontfamily='monospace', va='top', ha='left',
                zorder=25, transform=ax.transData)

        # subtexto: tempo total solicitado
        ax.text(px + pad * 0.2, py - 0.038 * yspan,
                f"/ {t_sim:.4f} s  (tempo físico)",
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

    def _save_frame(step, th, om, stop_label=None,
                    B_ext_inst=None, phi_ext_inst=None):
        """
        Renderiza e salva um frame PNG.

        B_ext_inst : valor COM SINAL do campo projetado na direção phi [T].
                     Negativo na fase inversa da histerese.
                     O painel mostra a seta invertida e o valor com sinal.
        phi_ext_inst: direção base do campo [graus].
        B_ext (closure) : amplitude máxima — usada para escalar a barra.
        """
        nonlocal n_frames
        t_phys   = step * dt
        S        = np.abs(np.mean(np.exp(1j * th)))
        om_max   = np.max(np.abs(om))

        # campo com sinal: B_now pode ser negativo (campo invertido)
        b_now    = B_ext_inst   if B_ext_inst   is not None else B_ext
        phi_now  = phi_ext_inst if phi_ext_inst is not None else phi_ext_deg

        # texto com sinal: "+60.8 mT" ou "-60.8 mT" ou "0 T"
        if b_now is not None and abs(b_now) > 1e-12:
            sign_str = "+" if b_now > 0 else "-"
            b_str = sign_str + _fmt_B(abs(b_now))
        else:
            b_str = "0 T"

        title = f"S = {S:.4f}   w_max = {om_max:.2f} rad/s   B = {b_str}"

        # quando B < 0: a seta aponta na direção oposta (phi + 180°)
        phi_display = phi_now if (b_now is None or b_now >= 0) else (phi_now + 180.0) % 360.0

        fig, ax = plot_state(th, xs, ys, title=title,
                             needle_len=needle_len, needle_width=needle_width,
                             r_halo=r_halo,
                             B_ext=abs(b_now) if b_now is not None else 0.0,
                             phi_ext_deg=phi_display,
                             B_ext_max=B_ext,
                             B_signed=b_now,        # valor com sinal para o texto do painel
                             figsize_inches=figsize_inches)
        plt.tight_layout()
        _draw_clock(ax, t_phys, needle_len, stop_label=stop_label)

        fpath = os.path.join(frame_dir, f"frame_{n_frames:05d}.png")
        plt.savefig(fpath, dpi=frame_dpi, bbox_inches='tight', facecolor='#1A1A2E')
        plt.close(fig)
        n_frames += 1

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

    # ── campo externo base: direção e amplitude máxima ────────────────────
    Bext_x0, Bext_y0 = ext_field        # campo base (direção + amplitude)
    B_max = np.sqrt(Bext_x0**2 + Bext_y0**2)   # amplitude máxima [T]
    phi_rad = np.arctan2(Bext_y0, Bext_x0)      # direção [rad]
    cos_phi = np.cos(phi_rad)
    sin_phi = np.sin(phi_rad)

    def field_at(t):
        """
        Retorna (Bx, By) do campo externo no instante físico t [s].

        Modos
        -----
        'static'
            Campo constante = (Bext_x0, Bext_y0) durante toda a simulação.

        'hysteresis'
            Ciclo completo em 5 segmentos iguais de t_sim/5 cada:
              segmento 1: t ∈ [0,       t5]   →  B: 0 → +B_max
              segmento 2: t ∈ [t5,    2*t5]   →  B: +B_max → 0
              segmento 3: t ∈ [2*t5,  3*t5]   →  B: 0 → -B_max
              segmento 4: t ∈ [3*t5,  4*t5]   →  B: -B_max → 0
              segmento 5: t ∈ [4*t5,  t_sim]  →  B: 0 → +B_max
            Direção: phi_ext (constante). Amplitude máxima: B_max = |ext_field|.

        'sine'
            B(t) = B_max · sin(2π · field_freq · t)
            Direção: phi_ext. Frequência: field_freq [Hz].

        'pulse'
            Fase 1 (campo aplicado): B = B_max (constante) até S ≥ 0.99.
            Fase 2 (relaxação): B = 0 após S ≥ 0.99, sistema relaxa livre.
            A transição entre fases é controlada pela flag _pulse_relaxing,
            atualizada no loop principal. field_at() consulta essa flag.
        """
        if field_mode == 'static':
            return Bext_x0, Bext_y0

        elif field_mode == 'hysteresis':
            # Ciclo completo: 0→+Bmax→0→-Bmax→0→+Bmax
            # 5 segmentos lineares iguais de t5 = t_sim/5
            T  = t_sim if t_sim > 0 else 1.0
            t5 = T / 5.0
            if t <= t5:
                # segmento 1: rampa 0 → +B_max
                B_scalar = B_max * (t / t5)
            elif t <= 2.0 * t5:
                # segmento 2: rampa +B_max → 0
                B_scalar = B_max * (1.0 - (t - t5) / t5)
            elif t <= 3.0 * t5:
                # segmento 3: rampa 0 → -B_max
                B_scalar = -B_max * ((t - 2.0 * t5) / t5)
            elif t <= 4.0 * t5:
                # segmento 4: rampa -B_max → 0
                B_scalar = -B_max * (1.0 - (t - 3.0 * t5) / t5)
            else:
                # segmento 5: rampa 0 → +B_max
                B_scalar = B_max * ((t - 4.0 * t5) / t5)
            return B_scalar * cos_phi, B_scalar * sin_phi

        elif field_mode == 'sine':
            B_scalar = B_max * np.sin(2.0 * np.pi * field_freq * t)
            return B_scalar * cos_phi, B_scalar * sin_phi

        elif field_mode == 'pulse':
            # fase 1: campo ligado; fase 2: campo desligado (relaxação)
            if _pulse_relaxing[0]:
                return 0.0, 0.0
            else:
                return Bext_x0, Bext_y0

        else:
            return Bext_x0, Bext_y0

    # flag mutável para o modo 'pulse' (lista de 1 elemento para ser
    # modificável dentro de _save_frame e do loop sem 'nonlocal')
    _pulse_relaxing = [False]   # False = campo ligado; True = relaxando

    # imprime instrução de uso
    _print("  [Ctrl+I para interromper e salvar o vídeo agora]")

    # ── função auxiliar: calcula torques SI para estado θ e campo B(t) ────
    def _torques(th, bx_ext, by_ext):
        """
        Retorna array N×M de torques [N·m] para o estado angular th.

        Parâmetros
        ----------
        th     : array N×M de ângulos [rad]
        bx_ext : componente x do campo externo instantâneo [T]
        by_ext : componente y do campo externo instantâneo [T]
        """
        tau = np.zeros((N, M))
        for i in range(N):
            for j in range(M):
                Bx, By = total_field_on(i, j, th, xs, ys, cutoff, moment,
                                        pbc=pbc, Lx=Lx, Ly=Ly, n_images=n_images)
                Bx += bx_ext
                By += by_ext
                mx = moment * np.cos(th[i, j])
                my = moment * np.sin(th[i, j])
                tau[i, j] = mx * By - my * Bx
        return tau

    # ── histórico de histerese: (t, B_scalar, M_x, M_y) ──────────────────
    # Gravado a cada passo; exportado para CSV ao final se field_mode=hysteresis
    field_log      = []   # log universal: (t, B, M_proj, S) para todos os modos
    hysteresis_log = [] if field_mode == 'hysteresis' else None
    sine_log       = [] if field_mode == 'sine'       else None

    # ── loop de integração Velocity-Verlet ────────────────────────────────
    bx_cur, by_cur    = field_at(0.0)           # campo no instante inicial
    tau_cur           = _torques(theta_cur, bx_cur, by_cur)
    _converged_count  = 0
    _S1_count         = 0   # passos consecutivos com S >= 0.9999 (parada em static)
    # S=1 como critério de parada só faz sentido em campo estático.
    # Em hysteresis e sine o campo oscila; em pulse queremos continuar
    # após S=1 para observar a relaxação com campo zerado.
    _allow_S1_stop     = (field_mode == 'static')
    # em hysteresis e sine o campo muda continuamente: parada por repouso
    # desabilitada (a rede pode estar momentaneamente parada no cruzamento B=0)
    _allow_rest_stop   = (field_mode in ('static', 'pulse'))
    _stop_reason       = "tempo total atingido"
    _pulse_phase   = "field_on"   # 'field_on' → 'relaxing' → 'done'
    _S99_count        = 0            # passos consecutivos com S ≥ 0.99 (pulse)

    for step in range(1, n_steps + 1):

        t_now = step * dt          # tempo físico atual [s]

        # ── verifica interrupção por teclado (Ctrl+I) ──────────────────────
        if _stop_flag.is_set():
            _print(f"\n  Interrompido (Ctrl+I) em t={t_now:.4f}s  (passo {step}/{n_steps})")
            _stop_reason = "interrompido pelo usuário (Ctrl+I)"
            if frame_dir is not None:
                _save_frame(step, theta_cur, omega_cur,
                            stop_label="■ interrompido (Ctrl+I)")
            break

        # ── campo externo no instante t (pode variar com o tempo) ──────────
        bx_new, by_new = field_at(t_now)

        # ── passo 1: atualiza θ usando θ(t), ω(t), τ(t) ───────────────────
        accel_cur  = (tau_cur - damping * omega_cur) / inertia
        theta_new  = theta_cur + omega_cur * dt + 0.5 * accel_cur * dt**2
        theta_new  = (theta_new + np.pi) % (2.0 * np.pi) - np.pi

        # ── passo 2: calcula τ(t+dt) com o novo θ e o novo campo ──────────
        tau_new = _torques(theta_new, bx_new, by_new)

        # ── passo 3: atualiza ω (Velocity-Verlet implícito) ────────────────
        b_half    = damping * dt / (2.0 * inertia)
        omega_new = (omega_cur * (1.0 - b_half)
                     + dt * (tau_cur + tau_new) / (2.0 * inertia)) \
                    / (1.0 + b_half)

        theta_cur = theta_new
        omega_cur = omega_new
        tau_cur   = tau_new
        bx_cur    = bx_new
        by_cur    = by_new

        if callback:
            callback(step, theta_cur.copy(), omega_cur.copy())

        if step % 20 == 0:
            hist.append((theta_cur.copy(), omega_cur.copy()))

        # ── calcula S e ω_max ──────────────────────────────────────────────
        S_now     = np.abs(np.mean(np.exp(1j * theta_cur)))
        omega_max = np.max(np.abs(omega_cur))

        # ── registra log de histerese / senoidal ──────────────────────────
        # M_proj: componente da magnetização média na direção do campo
        # ── log universal: campo e magnetização a cada passo ────────────────
        mx_mean  = float(np.mean(np.cos(theta_cur)))
        my_mean  = float(np.mean(np.sin(theta_cur)))
        M_proj   = mx_mean * cos_phi + my_mean * sin_phi
        B_scalar = bx_cur * cos_phi + by_cur * sin_phi
        entry    = (t_now, B_scalar, M_proj, S_now)
        field_log.append(entry)
        if hysteresis_log is not None:
            hysteresis_log.append(entry)
        if sine_log is not None:
            sine_log.append(entry)

        # ── salva frame periódico ──────────────────────────────────────────
        if frame_dir is not None and step % frame_every == 0:
            # B_signed: projeção com sinal na direção do campo
            # (negativo quando o campo aponta contra phi_ext, ex: hysteresis)
            B_signed = bx_cur * cos_phi + by_cur * sin_phi   # [T], com sinal
            _save_frame(step, theta_cur, omega_cur,
                        B_ext_inst=B_signed, phi_ext_inst=phi_ext_deg)
            # imprime progresso a cada 50 frames salvos
            if n_frames % 50 == 0:
                _print(f"  frame {n_frames-1:4d}  t={t_now:.4f}s  B={B_signed*1e3:+.3f}mT  S={S_now:.4f}  w={omega_max:.2f}rad/s")

        # ── condição 1: S = 1.00 (só em campo estático) ──────────────────
        if _allow_S1_stop:
            if S_now >= 0.9999:
                _S1_count += 1
            else:
                _S1_count = 0

            if _S1_count >= 30:
                _print(f"\n  S = 1.00 em t={t_now:.4f}s  (passo {step}/{n_steps})")
                _stop_reason = f"S = 1.00 atingido em t = {t_now:.4f} s"
                _stop_flag.set()
                if frame_dir is not None:
                    _save_frame(step, theta_cur, omega_cur,
                                stop_label="★ S = 1.00  alinhamento total")
                break

        # ── modo pulse: detecta S ≥ 0.99 e desliga o campo ─────────────
        if field_mode == 'pulse' and _pulse_phase == 'field_on':
            if S_now >= 0.99:
                _S99_count += 1
            else:
                _S99_count = 0

            if _S99_count >= 20:
                # alinhamento atingido → desliga o campo e entra em relaxação
                _pulse_relaxing[0] = True
                _pulse_phase = 'relaxing'
                _print(f"\n  Pulso: S>=0.99 em t={t_now:.4f}s  campo zerado, relaxando")
                # atualiza tau com campo=0 para o próximo passo
                tau_cur = _torques(theta_cur, 0.0, 0.0)
                if frame_dir is not None:
                    # B_ext_inst=0 mostra campo desligado no painel,
                    # mas phi_ext_deg é mantido para mostrar a direção
                    _save_frame(step, theta_cur, omega_cur,
                                stop_label="campo zerado - relaxando",
                                B_ext_inst=0.0)

        # ── condição 2: rede em repouso (ω_max → 0) ───────────────────────
        # Só para em modos static e pulse; em hysteresis/sine o campo muda
        # e a rede pode ficar momentaneamente parada no cruzamento por zero.
        if _allow_rest_stop:
            if omega_max < omega0 * 1e-3:
                _converged_count += 1
            else:
                _converged_count = 0

            if _converged_count >= 50:
                _print(f"\n  Rede em repouso em t={t_now:.4f}s  (passo {step}/{n_steps})  S={S_now:.4f}")
                _stop_reason = f"rede em repouso em t = {t_now:.4f} s"
                _stop_flag.set()
                if frame_dir is not None:
                    _save_frame(step, theta_cur, omega_cur,
                                stop_label="● rede em repouso")
                break

    _stop_flag.set()   # garante que a thread de teclado encerra

    # ── exporta dados de histerese / senoidal para CSV ─────────────────────
    if hysteresis_log:
        import csv
        csv_path = "hysteresis_loop.csv"
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['t_s', 'B_T', 'M_proj', 'S'])
            w.writerows(hysteresis_log)
        _print(f"  Dados de histerese salvos: {csv_path}")
        _plot_hysteresis(hysteresis_log)

    if sine_log:
        import csv
        csv_path = "sine_field.csv"
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['t_s', 'B_T', 'M_proj', 'S'])
            w.writerows(sine_log)
        _print(f"  Dados de campo senoidal salvos: {csv_path}")
        _plot_sine(sine_log, field_freq)

    return theta_cur, omega_cur, hist, n_frames, dt, _stop_reason, field_log

def next_available_path(path):
    """
    Sempre usa sufixo numérico de 4 dígitos, nunca cria o arquivo sem sufixo.
    Encontra o próximo número disponível na sequência verificando apenas
    o arquivo de vídeo (.mp4). O CSV correspondente sempre usa o mesmo número.

    Exemplos:
        sim.mp4         →  sim0000.mp4  (primeira vez)
                        →  sim0001.mp4  (segunda vez)
        sim0003.mp4     →  mesma série "sim", encontra próximo disponível
        dominios.mp4    →  dominios0000.mp4

    Parâmetros
    ----------
    path : str — nome desejado (com ou sem sufixo numérico)

    Retorna
    -------
    str — próximo caminho disponível na série nome0000, nome0001, ...
    """
    import os, re

    base, ext = os.path.splitext(path)

    # remove sufixo numérico final (4 dígitos) do base, se houver
    # "sim0008" → "sim",  "dominios0003" → "dominios",  "sim" → "sim"
    prefix = re.sub(r'\d{4}$', '', base)

    # encontra o próximo número na sequência verificando só o arquivo de vídeo
    n = 0
    while True:
        candidate = f"{prefix}{n:04d}{ext}"
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
        _print("AVISO: ffmpeg não encontrado. Instale com:")
        _print("  Ubuntu/Debian : sudo apt install ffmpeg")
        _print("  macOS (brew)  : brew install ffmpeg")
        _print("  conda         : conda install -c conda-forge ffmpeg")
        _print("  Windows       : https://ffmpeg.org/download.html")
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

    _print(f"\nMontando vídeo MP4: {output_path}")
    _print(f"  Fonte : {frame_dir}/frame_*.png   FPS={fps}   CRF={crf}")

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
            _print(f"  Codec : {desc}")
            _print(f"  Salvo : {output_path}  ({size_mb:.1f} MB)")
            return True
        else:
            last_err = result.stderr.strip().splitlines()
            short_err = last_err[-1] if last_err else "(sem mensagem)"
            _print(f"  [{desc}] falhou: {short_err}")

    _print(f"\nERRO: ffmpeg não conseguiu gerar o vídeo. Erro completo:")
    _print(result.stderr[-600:])
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
        _print("AVISO: ffmpeg não encontrado. Instale com:")
        _print("  Ubuntu/Debian : sudo apt install ffmpeg")
        _print("  macOS (brew)  : brew install ffmpeg")
        _print("  conda         : conda install -c conda-forge ffmpeg")
        _print("  Windows       : https://ffmpeg.org/download.html")
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

    _print(f"\nMontando vídeo MP4: {output_path}")
    _print(f"  Fonte : {frame_dir}/frame_*.png   FPS={fps}   CRF={crf}")

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
            _print(f"  Codec : {desc}")
            _print(f"  Salvo : {output_path}  ({size_mb:.1f} MB)")
            return True
        else:
            # extrai apenas a última linha do erro para não poluir o terminal
            last_err = result.stderr.strip().splitlines()
            short_err = last_err[-1] if last_err else "(sem mensagem)"
            _print(f"  [{desc}] falhou: {short_err}")

    # todas as estratégias falharam — mostra erro completo da última tentativa
    _print(f"\nERRO: ffmpeg não conseguiu gerar o vídeo. Erro completo:")
    _print(result.stderr[-600:])
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
    Lx, Ly  : dimensões do período da rede em x e y  [m]
              (usadas para condições periódicas de contorno, --pbc)
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
        Lx, Ly = M * d, N * d   # período da rede para PBC [m]
        return xs, ys, thetas, d, Lx, Ly

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
        Lx, Ly = M * dx_col, N * dy_row   # período da rede para PBC [m]
        return xs, ys, thetas, d, Lx, Ly

    # ── grade honeycomb (colmeia) ──────────────────────────────────────────
    elif geometry == 'honeycomb':
        # ── Geometria da imagem ──────────────────────────────────────────
        # Linhas COMPLETAS alternando com linhas de MEIA DENSIDADE,
        # criando buracos hexagonais visíveis.
        #
        # Linha par  (row=0,2,4,...): M agulhas, passo 2R, offset 0
        #   x = 0, 2R, 4R, ..., (M-1)*2R
        #
        # Linha ímpar(row=1,3,5,...): ~M/2 agulhas, passo 4R, offset R
        #   x = R, 5R, 9R, ...
        #
        # Δy entre linhas consecutivas = R*√3
        #
        # Cada agulha da linha completa toca:
        #   - 2 vizinhos na mesma linha (a dist 2R)
        #   - 1 vizinho na linha ímpar acima  (a dist 2R)  ← B(-R, +dy)
        #   - 1 vizinho na linha ímpar abaixo (a dist 2R)
        #
        # Verificação:  A(0,0) → B(-R, dy):
        #   d = √(R² + 3R²) = 2R  ✓
        #
        # Estratégia de preenchimento:
        #   Gera rede maior (N+4 pares de linhas) e recorta o retângulo W×H.
        # ─────────────────────────────────────────────────────────────────
        dy  = R * np.sqrt(3.0)   # espaçamento vertical entre linhas [m]
        d   = 2.0 * R            # distância entre vizinhos [m]

        # Dimensões alvo do retângulo
        W = (M - 1) * 2.0 * R   # largura: M agulhas na linha completa
        H = (N - 1) * dy         # altura:  N linhas com espaçamento dy

        # Gera rede maior (padding de 2 linhas/colunas em cada lado)
        N_rows = (N + 4) * 2     # * 2 porque alternamos par e ímpar
        x_start = -2.0 * 2 * R   # começa 2 agulhas à esquerda
        y_start = -2.0 * dy       # começa 2 linhas abaixo

        xs_list, ys_list = [], []
        for row in range(N_rows):
            y = y_start + row * dy
            if row % 2 == 0:
                # linha completa: x = x_start, x_start+2R, x_start+4R, ...
                x = x_start
                while x <= W + 2 * 2 * R:
                    xs_list.append(x)
                    ys_list.append(y)
                    x += 2.0 * R
            else:
                # linha meia: x = x_start+R, x_start+5R, x_start+9R, ...
                # offset de R em relação à linha par, passo 4R
                x = x_start + R
                while x <= W + 2 * 2 * R:
                    xs_list.append(x)
                    ys_list.append(y)
                    x += 4.0 * R

        all_x = np.array(xs_list)
        all_y = np.array(ys_list)

        # Recorta: mantém pontos dentro do retângulo alvo com margem R
        margin = R * 0.99
        mask = ((all_x >= -margin) & (all_x <= W + margin) &
                (all_y >= -margin) & (all_y <= H + margin))
        clipped_x = all_x[mask]
        clipped_y = all_y[mask]

        n_pts = len(clipped_x)
        xs     = clipped_x.reshape(n_pts, 1)
        ys     = clipped_y.reshape(n_pts, 1)
        thetas = noise * np.random.randn(n_pts, 1)
        Lx, Ly = W, H   # período aproximado da rede para PBC [m]
        return xs, ys, thetas, d, Lx, Ly

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
                              needle_len, B_ext_max=None, B_signed=None,
                              color='#FFD700'):
    """
    Representa o campo externo na imagem.  SEMPRE desenhado, mesmo quando
    B_ext = 0 (campo desligado após pulso, etc.).

    Elementos
    ---------
    1. SETAS NOS SÍTIOS — setas douradas em cada agulha apontando na
       direção do campo. Quando B_ext = 0, as setas ficam invisíveis
       (alpha=0) mas o painel de referência continua visível.

    2. PAINEL DE CAMPO — caixa no canto inferior direito com:
         • Seta mostrando a direção (comprimento proporcional a B_ext/B_max)
         • Texto: intensidade atual em SI  e  ângulo em graus
         • Barra de intensidade (verde cheia = B_max, vazia = 0)
       Quando B_ext = 0: seta apagada, texto "0 T", barra vazia.
       A direção exibida é sempre phi_ext_deg (direção do campo máximo).

    Parâmetros
    ----------
    ax          : eixo matplotlib
    xs, ys      : arrays 2D de posições  [m]
    B_ext       : intensidade ATUAL do campo  [T]  (pode ser 0)
    phi_ext_deg : direção do campo  [graus]
    needle_len  : comprimento de referência das agulhas  [m]
    B_ext_max   : intensidade máxima do campo  [T]  (para escalar a seta)
                  Se None, usa B_ext como máximo.
    color       : cor das setas e do painel (padrão: dourado)
    """
    phi      = np.deg2rad(phi_ext_deg)
    cos_phi  = np.cos(phi)
    sin_phi  = np.sin(phi)
    B_max_ref = B_ext_max if (B_ext_max and B_ext_max > 0) else max(B_ext, 1e-30)
    frac      = min(abs(B_ext) / B_max_ref, 1.0)   # 0→1, proporção do campo máximo
    field_on  = (B_ext > 1e-30)                     # True se campo não-nulo

    # ── 1. setas nos sítios (só quando campo está ligado) ─────────────────
    if field_on:
        alen = needle_len * 0.35
        dx_a = alen * cos_phi
        dy_a = alen * sin_phi
        N, M = xs.shape
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

    # ── 2. painel de campo (SEMPRE visível) ───────────────────────────────
    xlim  = ax.get_xlim()
    ylim  = ax.get_ylim()
    xspan = xlim[1] - xlim[0]
    yspan = ylim[1] - ylim[0]

    cx_d = xlim[0] + 0.88 * xspan
    cy_d = ylim[0] + 0.09 * yspan
    slen = needle_len * 1.1
    pad  = needle_len * 0.55

    # caixa de fundo
    rect = plt.Rectangle(
        (cx_d - slen * 0.7, cy_d - pad * 0.5),
        slen * 1.4, pad * 4.2,
        facecolor='#0A0A1A', edgecolor=color,
        linewidth=1.0, alpha=0.80, zorder=19,
        transform=ax.transData)
    ax.add_patch(rect)

    # seta de direção (comprimento proporcional à intensidade atual)
    arrow_len = slen * max(frac, 0.0)
    if arrow_len > slen * 0.05:   # só desenha se visível
        x_tail = cx_d - arrow_len * cos_phi / 2.0
        y_tail = cy_d - arrow_len * sin_phi / 2.0
        x_head = cx_d + arrow_len * cos_phi / 2.0
        y_head = cy_d + arrow_len * sin_phi / 2.0
        arrow_color = color if field_on else '#555555'
        ax.annotate(
            '', xy=(x_head, y_head), xytext=(x_tail, y_tail),
            arrowprops=dict(
                arrowstyle='->', color=arrow_color, lw=2.5,
                mutation_scale=18),
            zorder=20)
    else:
        # campo zero: desenha um ponto no centro
        ax.plot(cx_d, cy_d, 'o', ms=4, color='#555555', zorder=20)

    # barra de intensidade (fração do campo máximo)
    bar_w = slen * 1.2
    bar_h = pad * 0.22
    bar_y = cy_d - pad * 0.3
    # trilha
    ax.add_patch(plt.Rectangle(
        (cx_d - bar_w/2, bar_y), bar_w, bar_h,
        facecolor='#252540', edgecolor='none', zorder=20,
        transform=ax.transData))
    # preenchimento proporcional
    if frac > 0.01:
        ax.add_patch(plt.Rectangle(
            (cx_d - bar_w/2, bar_y), bar_w * frac, bar_h,
            facecolor=color if field_on else '#555555',
            edgecolor='none', alpha=0.8, zorder=21,
            transform=ax.transData))
    # borda da barra
    ax.add_patch(plt.Rectangle(
        (cx_d - bar_w/2, bar_y), bar_w, bar_h,
        facecolor='none', edgecolor='#5A5A8A',
        linewidth=0.7, zorder=22, transform=ax.transData))

    # ── textos ────────────────────────────────────────────────────────────
    tx = cx_d
    ty = cy_d + pad * 1.5

    # usa B_signed se disponível (mostra sinal), senão usa B_ext (positivo)
    b_val = B_signed if B_signed is not None else (B_ext if field_on else 0.0)
    if abs(b_val) > 1e-12:
        sign_str = "+" if b_val >= 0 else "-"
        b_str = sign_str + _fmt_B(abs(b_val))
    else:
        b_str = "0 T  (desligado)" if not field_on else "0 T"

    phi_str = f"dir: {phi_ext_deg:.1f} graus"

    ax.text(tx, ty, b_str,
            color=color if field_on else '#777777',
            fontsize=9, fontweight='bold',
            fontfamily='monospace', ha='center', va='bottom', zorder=23)
    ax.text(tx, ty - pad * 0.85, phi_str,
            color=color if field_on else '#777777',
            fontsize=7, fontfamily='monospace',
            ha='center', va='bottom', zorder=23)
    ax.text(tx, bar_y - pad * 0.15, "B ext",
            color=color, fontsize=7, alpha=0.8,
            fontfamily='monospace', ha='center', va='top', zorder=23)


def plot_state(thetas, xs, ys, title="Rede de bussolas", show_order=True,
               needle_len=0.42, needle_width=0.10, r_halo=None,
               B_ext=0.0, phi_ext_deg=0.0, B_ext_max=None,
               B_signed=None, figsize_inches=None):
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
    thetas        : array N×M de ângulos das agulhas  [rad]
    xs, ys        : arrays N×M de posições  [m]
    title         : título da figura
    show_order    : se True, desenha os halos de parâmetro de ordem local
    needle_len    : comprimento das agulhas  [m]
    needle_width  : largura das agulhas  [m]
    r_halo        : raio dos halos  [m]
    B_ext         : intensidade do campo externo  [T]
    phi_ext_deg   : direção do campo externo  [graus]
    figsize_inches: tupla (largura, altura) em polegadas, ou None para
                    calcular automaticamente em função do número de agulhas

    Retorna
    -------
    fig, ax : objetos matplotlib
    """
    # ── tamanho da figura: escala com a extensão física da rede ──────────
    # Garante que agulhas tenham tamanho visual adequado independente de
    # quantas há na rede. A figura cresce com a rede mas é limitada.
    if figsize_inches is None:
        x_span = xs.max() - xs.min() + 2 * needle_len * 2.0
        y_span = ys.max() - ys.min() + 2 * needle_len * 4.0
        aspect = x_span / y_span if y_span > 0 else 1.0
        # base: 8 polegadas no lado menor, máximo 20 polegadas
        base = 8.0
        if aspect >= 1.0:
            fig_w = min(base * aspect, 20.0)
            fig_h = fig_w / aspect
        else:
            fig_h = min(base / aspect, 20.0)
            fig_w = fig_h * aspect
        figsize_inches = (max(fig_w, 6.0), max(fig_h, 6.0))

    fig, ax = plt.subplots(figsize=figsize_inches, facecolor='#1A1A2E')
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
    draw_ext_field_on_lattice(ax, xs, ys, B_ext, phi_ext_deg, needle_len,
                              B_ext_max=B_ext_max, B_signed=B_signed)

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
        _print(f"GIF salvo em: {save_gif}")
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
    _print(f"Parâmetro de ordem salvo em: {outpath}")


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
                          '(0.0 a 0.8). Padrao: 0.80  ->  agulha = 0.80 * 2R. '
                          'Valores fora do intervalo serao limitados (clamp).')
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
                      help='Tempo fisico total da simulacao [s] '
                           '(= soma dos passos dt integrados). '
                           'Em hysteresis: cobre 1 ciclo completo (0->Hmax->-Hmax->Hmax). '
                           'Em sine: deve cobrir varios periodos (>= 3/field_freq). '
                           'Padrao: 2.0 s')
    grp2.add_argument('--field_mode',
                      choices=['static', 'hysteresis', 'sine', 'pulse'],
                      default='static',
                      help='Modo do campo externo: '
                           'static=constante (padrao), '
                           'hysteresis=ciclo completo 0->+Bmax->0->-Bmax->0->+Bmax (5 rampas), '
                           'sine=senoidal Hmax*sin(2pi*f*t), '
                           'pulse=aplica campo ate S>=0.99 depois zera e relaxa. '
                           'Amplitude=--B_ext, direcao=--phi_ext.')
    grp2.add_argument('--field_freq', type=float, default=1.0,
                      help='Frequencia do campo senoidal [Hz]. '
                           'Apenas para --field_mode sine. Padrao: 1.0 Hz')
    grp2.add_argument('--dt_factor', type=float, default=0.05,
                      help='Fracao do periodo natural T0 usada como passo dt '
                           '(0.02-0.10; menor = mais preciso, mais lento)')
    grp2.add_argument('--noise', type=float, default=1.5,
                      help='Amplitude do ruido inicial nos angulos [rad]; '
                           '0=todas para +x, 3.14=aleatorio total')
    grp2.add_argument('--seed', type=int, default=42,
                      help='Semente do gerador aleatorio (reprodutibilidade)')
    grp2.add_argument('--pbc', type=int, choices=[0, 1], default=0,
                      help='Condicoes periodicas de contorno (PBC): '
                           '0=desligado (padrao, rede finita com bordas), '
                           '1=ligado (rede tratada como estrutura periodica '
                           'infinita em x e y; agulhas nas bordas interagem '
                           'com replicas do lado oposto via convencao de '
                           'imagem minima)')
    grp2.add_argument('--pbc_images', type=int, default=1,
                      help='Numero de replicas periodicas somadas de cada '
                           'lado, em cada direcao, quando --pbc 1. Padrao: 1 '
                           '(soma sobre uma grade de (2*1+1)^2=9 celulas). '
                           'Usado apenas se --pbc 1.')

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
                      metavar='NOME',
                      help='Gera video MP4 com este nome base (requer ffmpeg). '
                           'Extensao .mp4 e opcional e adicionada automaticamente. '
                           'Sempre recebe sufixo numerico: nome0000.mp4, '
                           'nome0001.mp4, etc. Ex: --video simples '
                           '-> simples0000.mp4')
    grp4.add_argument('--frame_every', type=int, default=5,
                      help='Salva um frame a cada N passos (menor = mais suave, '
                           'mais lento). Padrao: 5')
    grp4.add_argument('--fps', type=int, default=24,
                      help='Quadros por segundo do video MP4. Padrao: 24')
    grp4.add_argument('--dpi', type=int, default=120,
                      help='Resolucao dos frames em pontos por polegada. '
                           'Padrao: 120. Para redes grandes (>15x15) use '
                           '150-200 para melhor qualidade de imagem.')
    grp4.add_argument('--keep_frames', action='store_true',
                      help='Mantém a pasta de PNGs intermediários após gerar o MP4')
    grp4.add_argument('--csv_order', choices=['t', 'B'], default='t',
                      help='Ordem das colunas no CSV exportado: '
                           't (padrao) = tempo na 1a coluna: t,B,M_proj,S. '
                           'B = tempo na ultima coluna: B,M_proj,S,t.')
    args = parser.parse_args()

    # ── valida e limita needle_frac ao intervalo [0, 0.8] ──────────────────
    # acima de 0.8 as agulhas começam a se tocar/sobrepor visualmente
    if args.needle_frac < 0.0 or args.needle_frac > 0.8:
        clamped = max(0.0, min(args.needle_frac, 0.8))
        _print(f"  Aviso: --needle_frac {args.needle_frac} fora do intervalo "
                f"[0, 0.8]; ajustado para {clamped}")
        args.needle_frac = clamped

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

    _print(f"\n{'═'*62}")
    _print(f"  Rede         : {args.geometry}  {args.N}x{args.M} agulhas")
    pbc_str = (f"ligado  (n_images={args.pbc_images}, soma sobre "
               f"{(2*args.pbc_images+1)**2} celulas)") if args.pbc else "desligado"
    _print(f"  PBC          : {pbc_str}")
    _print(f"  Raio R       : {args.R*100:.2f} cm  (2R = {2*args.R*100:.2f} cm)")
    _print(f"  Agulha       : {args.needle_frac*100:.0f}% de 2R = {args.needle_frac*2*args.R*100:.2f} cm")
    _print(f"  Momento mag. : {args.moment:.4g} A.m2  por agulha")
    _print(f"  Inercia      : {args.inertia:.3e} kg.m2  por agulha")
    _print(f"  Amortecimento: {args.damping:.3e} N.m.s/rad  (ar)")
    _print(f"  Campo externo: {fmt_field(B_ext_mag)}  phi={phi_ext_deg:.1f} graus")
    _print(f"  Componentes  : Bx={fmt_field(abs(Bext_x))}  By={fmt_field(abs(Bext_y))}")
    field_mode_str = args.field_mode
    if args.field_mode == 'sine':
        field_mode_str = f"sine  f={args.field_freq:.3f} Hz"
    elif args.field_mode == 'hysteresis':
        field_mode_str = "hysteresis  (0->+Bmax->0->-Bmax->0->+Bmax)"
    elif args.field_mode == 'pulse':
        field_mode_str = "pulse  (campo -> S>=0.99 -> zera -> relaxa)"
    _print(f"  Modo campo   : {field_mode_str}")
    _print(f"  t_sim        : {args.t_sim:.3f} s  dt_factor={args.dt_factor}")
    _print(f"  Ruido/seed   : {args.noise:.2f} rad  seed={args.seed}")
    _print(f"{'═'*62}\n")

    # ── gera grade de posições e ângulos iniciais ─────────────────────────
    # make_grid usa R como parâmetro único; retorna também nn_dist = 2R
    # e Lx, Ly = período da rede para condições periódicas de contorno
    xs, ys, thetas_init, nn_dist, Lx_period, Ly_period = make_grid(
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

    # ── PBC: limita cutoff a min(Lx,Ly)/2 para evitar contagem dupla ───────
    # Acima desse limite, a convenção de imagem mínima pode contar a mesma
    # réplica mais de uma vez, gerando resultados fisicamente incorretos.
    if args.pbc:
        max_cutoff_pbc = min(Lx_period, Ly_period) / 2.0
        if cutoff > max_cutoff_pbc:
            _print(f"  PBC: cutoff reduzido de {cutoff*100:.2f}cm para "
                    f"{max_cutoff_pbc*100:.2f}cm (limite min(Lx,Ly)/2)")
            cutoff = max_cutoff_pbc

    # campo dipolar de referência entre primeiros vizinhos (para exibição)
    B_ref = MU0_OVER_4PI * 2.0 * args.moment / nn_dist**3
    _print(f"  B_dipolar ref: {fmt_field(B_ref)}  (entre vizinhos)")
    if B_ext_mag > 0:
        ratio = B_ext_mag / B_ref
        dom = "DOMINANTE" if ratio > 1 else "fraco"
        _print(f"  B_ext/B_ref  : {ratio:.3f}  (campo externo {dom})")
    _print()

    ext_kwargs = dict(B_ext=B_ext_mag, phi_ext_deg=phi_ext_deg,
                      r_halo=r_halo)

    # ── tamanho da figura proporcional à extensão da rede ─────────────────
    # Calculado uma vez e reutilizado em todos os frames e figuras estáticas.
    # A lógica espelha a de plot_state mas usa xs/ys já conhecidos.
    _x_span = (xs.max() - xs.min()) + 4 * needle_len
    _y_span = (ys.max() - ys.min()) + 6 * needle_len
    _aspect = _x_span / _y_span if _y_span > 0 else 1.0
    _base   = 8.0
    if _aspect >= 1.0:
        _fw = min(_base * _aspect, 20.0)
        _fh = _fw / _aspect
    else:
        _fh = min(_base / _aspect, 20.0)
        _fw = _fh * _aspect
    figsize_main = (max(_fw, 6.0), max(_fh, 6.0))

    ext_kwargs['figsize_inches'] = figsize_main

    # ── figura do estado inicial ──────────────────────────────────────────
    fig0, _ = plot_state(thetas_init, xs, ys,
                         title="Estado inicial (aleatório)",
                         needle_len=needle_len, needle_width=needle_width,
                         **ext_kwargs)
    plt.tight_layout()
    plt.savefig("compass_initial.png",
                dpi=args.dpi, bbox_inches='tight', facecolor='#1A1A2E')
    plt.close(fig0)
    _print("Estado inicial salvo.")

    # ── integração dinâmica inercial ──────────────────────────────────────
    frame_dir   = None
    final_video = None
    if args.video:
        import os
        # se o usuário não informou extensão, adiciona .mp4 automaticamente
        # ex: "--video simples" -> "simples.mp4" -> "simples0000.mp4"
        if not os.path.splitext(args.video)[1]:
            args.video = args.video + ".mp4"
        # resolve o nome final do vídeo ANTES de criar a pasta de frames,
        # garantindo que pasta e arquivo usem o mesmo nome base
        final_video = next_available_path(args.video)
        _print(f"  Video sera salvo como '{final_video}'")
        base      = os.path.splitext(final_video)[0]
        frame_dir = base + "_frames"
        _print(f"Integrando e gravando frames em '{frame_dir}/'...")
    else:
        _print("Integrando dinâmica inercial...")

    thetas_eq, omegas_eq, thetas_hist, n_frames, sim_dt, stop_reason, field_log = relax(
        thetas_init.copy(), xs, ys,
        t_sim=args.t_sim,
        dt_factor=args.dt_factor,
        inertia=args.inertia,
        damping=args.damping,
        cutoff=cutoff,
        ext_field=(Bext_x, Bext_y),
        moment=args.moment,
        field_mode=args.field_mode,
        field_freq=args.field_freq,
        frame_dir=frame_dir,
        frame_every=args.frame_every,
        needle_len=needle_len,
        needle_width=needle_width,
        r_halo=r_halo,
        frame_dpi=args.dpi,
        figsize_inches=figsize_main,
        pbc=bool(args.pbc),
        Lx=Lx_period,
        Ly=Ly_period,
        n_images=args.pbc_images,
        B_ext=B_ext_mag,
        phi_ext_deg=phi_ext_deg,
    )
    frames_str = f"  ({n_frames} frames salvos)" if n_frames else ""
    _print(f"Integracao concluida - {stop_reason}{frames_str}")

    # ── exporta CSV universal com campo e magnetização ────────────────────
    # Nome = mesmo nome do vídeo (se houver), senão "compass_field_log.csv"
    import csv as _csv, os as _os
    if args.video and final_video:
        csv_path = _os.path.splitext(final_video)[0] + ".csv"
    else:
        csv_path = "compass_field_log.csv"

    # ── ordem das colunas conforme --csv_order ────────────────────────────
    # field_log armazena tuplas (t, B, M_proj, S)
    # 't' (padrao): tempo na 1a coluna  -> t, B, M_proj, S
    # 'B'         : tempo na ultima coluna -> B, M_proj, S, t
    if args.csv_order == 'B':
        header = ['B_aplicado_T', 'M_proj', 'S', 't_s']
        rows   = [(B, Mp, S, t) for (t, B, Mp, S) in field_log]
    else:
        header = ['t_s', 'B_aplicado_T', 'M_proj', 'S']
        rows   = field_log

    with open(csv_path, 'w', newline='') as f:
        w = _csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    _print(f"  CSV salvo: {csv_path}  ({len(field_log)} pontos)")

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
                dpi=args.dpi, bbox_inches='tight', facecolor='#1A1A2E')
    plt.close(fig1)
    _print("Estado de equilíbrio salvo.")

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
                dpi=args.dpi, bbox_inches='tight', facecolor='#1A1A2E')
    plt.close(fig2)
    _print("Comparação salva.")

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
            _print(f"  Pasta de frames removida: {frame_dir}/")
        elif not ok:
            _print(f"  Os frames PNG foram mantidos em: {frame_dir}/")

    _print("\nConcluído.")


if __name__ == '__main__':
    main()
