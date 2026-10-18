# Cronologia do projeto: do primeiro script ao Relatório V60

Data: 18/10/2026. Guia reconstruído a partir de quatro fontes:

1. **Código arquivado:** 76 versões em `archive/OLD/` e em "simula agulhas/OLD" no Google Drive (`compass_simT*.py`, `compass_simV11–V77.py`), mais `archive/compass.py.bkp` (V79) e `archive/compassV02*.py` / `compassV2_1 (1).py` (V2.x). As mudanças de cada versão foram **extraídas automaticamente**: flags de linha de comando adicionadas e removidas, modos de campo, funções novas e tamanho do código. Nos poucos casos sem registro escrito, foram conferidas por `diff`.
2. **Registro de versões do código V1–V65** e **narrativa das sessões de junho**, no relatório `relatorio_projeto_V22_EN.tex` (05/07).
3. **Relatórios e guias de 08–11/07** na raiz da pasta do Drive, com datas de modificação dos arquivos.
4. **Histórico do git**, de 12/07 em diante.

As datas vêm da modificação dos arquivos ou dos commits. Onde uma afirmação é inferência, isso está indicado.

---

## 0. Atenção: há três numerações diferentes

| Linha | Numeração | Exemplo de ambiguidade |
|---|---|---|
| **Código do simulador** | T1–T10 → V11–V79 → V2.0.1, V2.0.2, V2.1 → `compass.py` 2.1.0 com correções (git) | "**V60**" do código é uma versão de **05/07**: a tradução da interface de linha de comando para o inglês. |
| **Relatório** | arquivos `relatorio_projeto_V09`…`V22` (04–05/07, com versões internas V01–V32 na tabela do próprio relatório) → `project_report` **V51, V52, V53** (10–12/07) → **V60** (17/10) | O **Relatório V60** atual não tem relação com o código V60. Não existiram relatórios V54–V59: o salto de V53 para V60 marcou a reestruturação de outubro. Entre o relatório interno V32 (05/07) e o V51 (10/07) não há arquivos intermediários; a numeração foi provavelmente avançada para casar com a linha do código. |
| **Git** | commits e PRs #4–#8 | Os PRs são a referência mais confiável a partir de 12/07. |

---

## 1. Linha do tempo resumida

| Período | Versões | O que mudou |
|---|---|---|
| **15/06** | T1–T2 | Núcleo: agulhas em losango com integração por velocity Verlet e momento de inércia, interação dipolar completa, geometrias quadrada, triangular e honeycomb, campo externo, vídeo MP4. Comentários em português; o cabeçalho ainda descrevia a dinâmica como "inspirada em LLG". |
| **16/06** | T3–T10 | Rede parametrizada por R (`--R`, `--needle_frac`); a partir da T4, dinâmica **newtoniana inercial** explícita; modos de campo `static`, `hysteresis`, `sine`, `saturate`, `pulse`. |
| **16–17/06** | V11–V39 | Registro formal de versões a partir da V20. CSV, PBC, GPU (CuPy), I e m automáticos pela geometria, renderização vetorizada (17–30×), domínios magnéticos, `--make_images`. |
| **18/06** | — | Revisão de literatura e estratégia de publicação (relatório V22, Fase 6). |
| **24/06** | `compass_sim_gpu.py/.ipynb` | Variante para GPU, em notebook. |
| **04/07** | V40–V44 | Tensor dipolar pré-calculado (35–45×), novos protocolos (histerese log, pulso temporizado, degraus), `--t_sim_full`, painel de progresso. Relatórios internos V01–V09. |
| **05/07** | V45–V65 | Vídeo com painel duplo M×H; tradução para o inglês (V47–V65, com queda de 3945 para 2430 linhas); interface (HUD); `--damping_noise` (V59); ajustes dos scripts de varredura. Relatórios V10–V32 (arquivos V12–V22). |
| **06–07/07** | V66–V77 | **FORC** (V67–V71, com taxa constante na V71); **desmagnetização** (V70–V76: rotacional, `--demag_temp`, `--demag_delay`); tempo de execução no log (V77). |
| **08–09/07** | V78–V79 | Correção do pivô em I e m (`--pivot_*`). Resumos para o LAW3M 2026; manual de varreduras; guia do usuário. |
| **10/07** | V2.0.1 → V2.0.2 | **Ramo de pesquisa reescrito** (1292 linhas): cutoff explícito, observáveis S1/S2, contadores de flips, energias, metadados JSON, estados NPZ, PNGs automáticos. Relatórios V51 e V52; `damping_sweepV03.py`. |
| **11/07** | V2.1 | Contadores de avalanche robustos a oscilação (Schmitt + dwell), monitor de estabilidade, `--event_log`. |
| **12/07** | git (`ad865d7`…`cb14300`) | Repositório criado; auditoria de código; `config.yaml`; testes de regressão; motor único `compass.py`; janela lenta na histerese; Relatório V53. |
| **16/10** | git, PRs #4–#6 | **Auditoria física** e **Fase 1**: soma completa, dt pela cota do campo local, Q₀ intrínseco, livro-caixa de energia, taxa adimensional r, testes físicos, README, licença MIT. |
| **16–17/10** | git, PR #7 (aberto) | **Fase 2**: 690 ciclos com taxa controlada (Mac + mumax + GPU); protocolo quase-estático por eventos; manual com PDF gerado; **Relatório V60**. |
| **18/10** | git, PR #8 + `21d8ec2` | Descoberta do **artefato de pontos de sela** na primeira grade quase-estática; correção por chute e nova relaxação; grade corrigida em execução. |

---

## 2. Detalhe por fase

### 2.1 Fase 1 de junho: núcleo da simulação (15/06, T1–T2)
Pedidos sucessivos estabeleceram:
- agulhas em losango (metade norte azul, metade sul branca), vistas de cima;
- interação dipolar entre **todas** as agulhas e busca do equilíbrio;
- parâmetros da rede (`--dx`, `--dy`, `--N`, `--M`, `--geometry`) e campo uniforme (`--B_ext`, `--phi_ext`, `--ext_Bx/By`);
- vídeo (`--video`, `--fps`, `--frame_every`).

Já no primeiro arquivo há integração por velocity Verlet e `--inertia`. Duas decisões físicas desse dia definiram o projeto:
- **unidades SI** em todo o código;
- **dinâmica inercial**, de início por realismo visual: a agulha oscila antes de parar. Mais tarde isso virou o diferencial científico do projeto.

### 2.2 Fase 2: controle, tempo físico e honeycomb (16–17/06, T3–V19)
- **T3:** a rede passa a ser parametrizada pelo raio R da circunferência da agulha (`--R`, `--needle_frac`, substituindo `--dx/--dy`). Isso garante r_nn = 2R igual nas três redes.
- **T4:** o cabeçalho deixa de mencionar LLG; a dinâmica newtoniana inercial passa a ser a descrição oficial.
- **T5–T9:** modos de campo `static`, `hysteresis`, `sine` (T5), `saturate` (T6) e `pulse` (T9). O ciclo de histerese foi corrigido para **cinco** segmentos.
- **T10:** `--dpi`.
- **V13:** `--csv_order`, exportação CSV de (t, B, M, S), que virou a espinha dorsal das análises.
- **V16–V17:** condições periódicas (`--pbc`, `--pbc_images`).
- **Outras correções do período:**
  - `t_sim` passou a ser o tempo físico simulado, não o tempo de relógio;
  - a honeycomb foi refeita como rede triangular sem o sítio central de cada hexágono;
  - interrupção com Ctrl+I.

### 2.3 Fases 3 a 5: instrumento de medida (17–18/06, V20–V39)
- **V20:** marco do **registro formal de versões** (V1–V19 consolidadas retroativamente).
- **V21–V24:** critério de platô no pulso; indicadores na tela.
- **V25–V30:** GPU via CuPy (V26: torque vetorizado; V30: `--gpu`).
- **V29–V33:**
  - imagens periódicas no campo dipolar;
  - barra de progresso (V32);
  - **momento de inércia automático** pela geometria real (V33: `--needle_thickness`, `--steel_density`; chapa de 0,4 mm);
  - Ctrl+C e Ctrl+I corrigidos.
- **V35:** renderização vetorizada (17–30×).
- **V36:** `--progress_bar`.
- **V37:** **momento magnético automático**, m = M_s·V (`--steel_Bsat`).
- **V38:** **domínios magnéticos** (`--halo_mode`, `--domain_tol`, union-find), mais tarde usados como medida.
- **V39:** `--make_images` (3× mais rápido em redes pequenas).

Nesse período também nasceram os scripts de campanha (`damping_sweep.py`, `damping_sweep_analysis.py`) e o estudo piloto de amortecimento.

### 2.4 Fase 6: de ferramenta a pesquisa (18/06 e 04–07/07, V40–V77)
- **V40 (04/07):** **tensor dipolar pré-calculado**, com ganho medido de 35–45× por corrida. O `damping_sweepV02.py` ganhou paralelização e retomada.
- **V41:** histerese com espaçamento log (`--hyst_spacing`, `--hyst_log_k`); pulso com `--field_delay`, `--t_pulse` e `--torque_tol`; modos `step_pos` e `step_neg`.
- **V42:** `--t_sim_full`.
- **V43:** halos também na figura de comparação.
- **V44:** painel de progresso de duas linhas.
- **V45:** vídeo de histerese com painel duplo (agulhas e curva M×H).
- **V46–V65 (05/07):** alinhamento da ajuda (V46); `--hyst_autoscale` e tradução dos comentários (V47); ajustes do HUD, mostrador e seta de direção (V49–V58); **`--damping_noise`** (V59); tradução da interface de linha de comando (**V60**), dos gráficos (V61) e dos avisos (V62); constantes pré-calculadas no laço (V63); imports e tradução dos scripts de varredura (V64–V65).
- **V66 (06/07):** só ajustes de formatação (alinhamento da ajuda e da assinatura de funções).
- **V67–V69:** **FORC**. Na V67, `--forc_Br`, `--forc_t_rev` e `--forc_t_sat`. Na V68, a família de curvas (`--forc_n_curves`, `--forc_Br_min`, `--forc_t_ramp_down/up`). Na V69, `--forc_sweep`.
- **V70:** **desmagnetização** (`--demag`, `--demag_cycles`, `--demag_freq`).
- **V71:** FORC com **taxa constante** (`--forc_rate`), em resposta à "armadilha da taxa variável" documentada em `FORC.md`.
- **V72:** desmagnetização rotacional com envelope linear, no lugar do decaimento geométrico 0,7ⁿ.
- **V73:** `--hyst_autoscale` ligado por padrão.
- **V74:** `--demag_temp` (opção de recozimento térmico).
- **V75:** barra de progresso com etapas ("demag" × "integrando").
- **V76:** `--demag_delay`.
- **V77 (07/07):** tempo de relógio no log final.

### 2.5 Consolidação e primeira divulgação (08–09/07, V78–V79)
- **V79:** correção do pivô no momento de inércia e no volume magnético (`--pivot_radius`, `--pivot_thickness`, `--pivot_density`, `--pivot_mass`). A V78 não foi preservada como arquivo separado.
- Resumos para o **LAW3M 2026** (`Abstract.tex` de 08/07; `LAW3M2026_compass_needle_abstract*.tex` de 09/07), manual das varreduras (`damping_sweep_manual.tex`) e guia do usuário.
- `calc_defaults_temp.py` usa a espessura de 0,26 mm, origem da inconsistência com os 0,4 mm do simulador que ainda está pendente.

### 2.6 Ramo de pesquisa V2 (10–11/07)
- **V2.0.1 (`compassV02.py`):** **reescrita completa**, de 3058 para 1292 linhas e 26 funções a menos. Ficaram de fora vídeo, halos, HUD e as opções de desmagnetização "térmica". Entraram:
  - cutoff explícito (`--cutoff_shells`, `--cutoff_m`);
  - dimensões físicas da agulha independentes do espaçamento (`--needle_len`, `--needle_width`);
  - observáveis S1, S2, diretor, q_axis, contadores de flips e energias;
  - metadados JSON e estados NPZ;
  - `--float32` e `--use_gpu`.
- **V2.0.2 (`compassV02_corrected.py`):** PNGs automáticos dos estados inicial e final (`--png_*`) e carimbo de data na fonte.
- **V2.1 (`compassV2_1.py`, 11/07):** contadores com histerese de Schmitt e tempo de permanência (`--flip_band_deg`, `--flip_dwell_T0`, `--flip_settle_frac`), monitor de estabilidade (`--dt_guard_alpha`, `--dt_guard_substep`) e `--event_log`.
- **Relatórios:** `project_report_V51_compassV2.tex` (10/07, 08h36) reorganizou o relatório em torno dos objetivos acadêmicos; o **V52** (10/07, 12h19) alinhou-o à V2.0.2. Guias do usuário V2.0.2. `damping_sweepV03.py` compatível com a V2.

### 2.7 Repositório git e auditoria de código (12/07)
Repositório criado (`ad865d7`), com auditoria de código (`docs/AUDIT.md`), parâmetros extraídos para `config.yaml`, testes de regressão e arquivamento das linhagens mortas em `archive/`. **A V2.1 virou o `compass.py` canônico** (`9a2cc63`). Ganhou ainda log adaptativo e janela de varredura lenta na histerese (`--log_adaptive`, `--hyst_slow_window`, `--hyst_slow_factor`). O `USER_GUIDE.md` foi reescrito, e o **Relatório V53** saiu na mesma data.

### 2.8 Auditoria física e Fase 1 (16/10; PRs #4, #5 e #6, em `main`)
- **Auditoria** (`docs/AUDIT_PHYSICS_ROADMAP.md`). O cutoff enviesava o campo em 15–17 %; o dt estava sub-resolvido (erro de energia de 1,4–1,8 %); Q dependia do protocolo (fator √8); a área do laço era dominada pela taxa (∝ taxa^~0,6). As figuras antigas de "área × Q" e "expoente × Q" não mediam o que diziam.
- **Correções**, cada uma com teste físico:
  - soma dipolar completa (`de455f6`);
  - dt pela cota do campo local (`31f435b`);
  - Q₀ intrínseco (`db21d33`);
  - livro-caixa W_field / E_diss (`814744f`);
  - taxa adimensional `--sweep_rate_dimless` (`7528531`, PR #5).
- **Documentação:** README, `requirements.txt`, manual revisado e licença MIT (PR #6).

### 2.9 Fase 2: histerese dinâmica e limite quase-estático (16–18/10; PRs #7 e #8)
- **Campanha com taxa controlada** (`fase2_campaign.py`), com divisão entre Mac e mumax e as redes grandes na GPU RTX 3090: 690 ciclos. Para Q₀ ≤ 1 há colapso com a taxa em unidades viscosas; em taxa finita, a inércia estreita o laço em 20–50 %.
- **Protocolo quase-estático por eventos** (`--field_mode hysteresis_qs`, `6312d4e`) e fluxo gradiente como referência Q₀ → 0. Manual com fonte única e PDF gerado por pandoc (`0d82bcf`).
- **Relatório V60** (`c06e880`, 17/10): reestruturado como relatório de métodos e resultados, com três seções de resultados. Sem comandos nem tabela de flags, que ficam no manual.
- **18/10 (PR #8):** a análise de seleção de caminho mostrou que boa parte do "efeito inercial" da primeira grade quase-estática era **artefato de pontos de sela**: o fluxo gradiente preserva a simetria da rede e fica no ramo instável. Em `21d8ec2`, o critério de repouso passou a exigir um **chute de 10⁻⁴ rad e nova relaxação**. A grade corrigida está em execução, e as seções de resultados III e o resumo do V60 serão reescritos com ela.

---

## 3. Onde está cada coisa

| Item | Local |
|---|---|
| Código T1–V77 | `archive/OLD/compass_sim*.py` (git) e "simula agulhas/OLD" (Drive) |
| V79, V2.0.1, V2.0.2, V2.1 | `archive/compass.py.bkp`, `archive/compassV02.py`, `archive/compassV02_corrected.py`, `archive/compassV2_1 (1).py` |
| Relatórios de julho (V09–V22, internos V01–V32) | `archive/OLD/relatorio_projeto_*` (git) e Drive/OLD |
| Relatórios V51–V53 | `archive/project_report_V51_compassV2.*`, `archive/project_report.*` (V51/V52) e o histórico do git de `docs/project_report.tex` (V53) |
| Relatório atual (V60) | `docs/project_report.tex` e `docs/project_report.pdf` |
| Transcrição das conversas de junho | `archive/OLD/compass_sim_conversation.md` (2,1 MB; não usada aqui, mas é a fonte primária das sessões de 15–18/06) |
| Histórico detalhado de 12/07 em diante | `git log` |
