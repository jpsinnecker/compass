# Cronologia do projeto: do primeiro script ao Relatório V60

Data: 18/10/2026 (revisão 6, com as fontes primárias, a pasta original do Dropbox, os dados da IDE Antigravity, o histórico local de arquivos do editor e o pacote `OLD/ALLfiles`). Guia reconstruído a partir de nove fontes:

1. **Código arquivado:** 76 versões em `archive/OLD/` e em "simula agulhas/OLD" no Google Drive (`compass_simT*.py`, `compass_simV11–V77.py`), mais `archive/compass.py.bkp` (V79) e `archive/compassV02*.py` / `compassV2_1 (1).py` (V2.x). As mudanças de cada versão foram **extraídas automaticamente**: flags de linha de comando adicionadas e removidas, modos de campo, funções novas e tamanho do código. Nos poucos casos sem registro escrito, foram conferidas por `diff`.
2. **Registro de versões do código V1–V65** e **narrativa das sessões de junho**, no relatório `relatorio_projeto_V22_EN.tex` (05/07).
3. **Relatórios e guias de 08–11/07** na raiz da pasta do Drive, com datas de modificação dos arquivos.
4. **Histórico do git**, de 12/07 em diante.
5. **Transcrições das sessões de 15–18/06 no claude.ai:** 99 pedidos datados, listados em `docs/history/requests_2026-06_to_07-05.csv`. A fonte completa está em `archive/OLD/compass_sim_conversation.md`, que reúne as mesmas 7 sessões.
6. **Sessões de 04–10/07 na IDE Antigravity** (`~/.gemini/antigravity-ide`, um link para `Dropbox/.antigravity-ide`): **13 conversas do compass, com 184 pedidos datados de 04/07 14:05 a 10/07 07:31** (`docs/history/requests_antigravity_2026-07.csv`). As transcrições completas estão em `archive/transcripts/antigravity_2026-07/conversations/`, ao lado do log `chat_antiIDE.md` e dos relatórios de tarefa (`*_walkthrough.md`, `*_task.md`). A pasta `scratch/texbuild` da IDE guardava as últimas compilações: dela vieram o **Relatório V47** e um rascunho do guia da V78.

7. **A pasta original do Dropbox** (`~/Library/CloudStorage/Dropbox/simula agulhas`), onde o repositório git foi criado em 12/07. O histórico do git é o mesmo do GitHub (30 commits), mas os arquivos apagados nele e os `.zip` de julho devolveram os fontes dos **relatórios V10 e V11**. Também havia arquivos não versionados de 13/07 (seção 2.7).

8. **Histórico local de arquivos do editor** (`~/Library/Application Support/Antigravity IDE/User/History`): duas cópias salvas do `compass.py`, de 05/07 21:31 (idêntica à V62) e de **08/07 07:56, que é a V78**, agora em `archive/compass_simV78.py`.

9. **O pacote `OLD/ALLfiles`** (`ALLfiles.zip`, montado em 06/07, no Drive): os **fontes dos relatórios V01–V09** de 04/07, que não tinham arquivo (com os PDFs V01–V04), e as **primeiras versões dos scripts de campanha** (`damping_sweep.py`, `damping_sweepV02.py`, `damping_sweep_analysis.py`, `sweep_damping_hysteresis.py`), ainda em português e importando `compass_sim`. Estão em `archive/OLD/relatorio_projeto_V01…V11.tex` e `archive/OLD/campaign_early/`. As figuras e os vídeos de junho do pacote ficaram fora do git (são resultados).

**Como navegar:** cada versão arquivada é um commit da branch `history/pre-git`, com a data original e uma tag. Ver a seção 4.

As datas vêm da modificação dos arquivos ou dos commits. Onde uma afirmação é inferência, isso está indicado.

---

## 0. Atenção: há três numerações diferentes

| Linha | Numeração | Exemplo de ambiguidade |
|---|---|---|
| **Código do simulador** | T1–T10 → V11–V79 → V2.0.1, V2.0.2, V2.1 → `compass.py` 2.1.0 com correções (git) | "**V60**" do código é uma versão de **05/07**: a tradução da interface de linha de comando para o inglês. |
| **Relatório** | numeração **contínua**: V01–V32 (04–05/07, arquivos `relatorio_projeto_V01`…`V22`) → V33–V50 (06–09/07, `project_report.tex` nas sessões da Antigravity; o V33 documentou o código V66 e o V47, a V77) → **V51, V52** (10/07) → **V53** (12/07) → **V60** (17/10) | O **Relatório V60** não tem relação com o código V60. Não existiram V54–V59: o salto de V53 para V60 marcou a reestruturação de outubro. Dos arquivos V33–V50, só o V47 foi preservado. |
| **Git** | commits, PRs #4–#8 e tags (`v2.2`, `v2.3.0`, `v2.4.0`, `report-Vnn`) | Os PRs são a referência mais confiável a partir de 12/07. A política de versões daqui em diante está no `CHANGELOG.md`. |

**Ferramentas e locais de armazenamento por período:**
- **15–18/06 e manhã de 04/07** (até a V44): claude.ai (chat), com os arquivos no Dropbox.
- **04/07 14:05 a 10/07:** IDE Antigravity, ainda no Dropbox (`~/Library/CloudStorage/Dropbox/simula agulhas`), cobrindo V45–V79, os relatórios V10–V50 e os resumos. Por dia: 04/07 (2 pedidos), 05/07 (55), 06/07 (50), 07/07 (31), 08/07 (12), 09/07 (22), 10/07 (12).
- **10–11/07:** pasta "simula agulhas" no Google Drive.
- **12/07 em diante:** repositório git; desde 16/10, Claude Code. O motor tem `__version__` próprio desde a 2.4.0.

---

## 1. Linha do tempo resumida

| Período | Versões | O que mudou |
|---|---|---|
| **15/06** | T1–T2 | Núcleo: agulhas em losango com integração por velocity Verlet e momento de inércia, interação dipolar completa, geometrias quadrada, triangular e honeycomb, campo externo, vídeo MP4. Comentários em português; o cabeçalho ainda descrevia a dinâmica como "inspirada em LLG". |
| **16/06** | T3–T10 | Rede parametrizada por R (`--R`, `--needle_frac`); a partir da T4, dinâmica **newtoniana inercial** explícita; modos de campo `static`, `hysteresis`, `sine`, `saturate`, `pulse`. |
| **16–17/06** | V11–V39 | Registro formal de versões a partir da V20. CSV, PBC, GPU (CuPy), I e m automáticos pela geometria, renderização vetorizada (17–30×), domínios magnéticos, `--make_images`. |
| **18/06** | — | Revisão de literatura e estratégia de publicação (relatório V22, Fase 6). |
| **24/06** | `compass_sim_gpu.py/.ipynb` | Variante para GPU, em notebook. |
| **04/07** | V40–V44 | Tensor dipolar pré-calculado (35–45×), novos protocolos (histerese log, pulso temporizado, degraus), `--t_sim_full`, painel de progresso. Relatórios V01–V09, um por mudança do código (fontes recuperados do `ALLfiles`). |
| **05/07** | V45–V65 | Vídeo com painel duplo M×H; tradução para o inglês (V47–V65, com queda de 3945 para 2430 linhas); interface (HUD); `--damping_noise` (V59); ajustes dos scripts de varredura. Relatórios V10–V32 (arquivos V12–V22). |
| **06–07/07** | V66–V77 | Barra de campo bidirecional no HUD (V66); **FORC** (V67–V71, com taxa constante na V71); **desmagnetização** (V70–V76: rotacional, `--demag_temp`, `--demag_delay`); tempo de execução no log (V77). Relatórios V33–V47. Sessões na Antigravity. |
| **08–09/07** | V78–V79 | **Geometria física da agulha real** (V79: R = 7,5 mm, espessura de 0,26 mm, pivô de latão com o furo descontado de I e m). Guia do usuário da V78, resumos para o LAW3M 2026, manual de varreduras, FORC nas varreduras. |
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

Nesse período também nasceram os scripts de campanha (`damping_sweep.py`, `damping_sweep_analysis.py`) e o estudo piloto de amortecimento. As versões mais antigas que sobreviveram, em português e importando `compass_sim`, estão em `archive/OLD/campaign_early/`; são anteriores à V64 (05/07), que passou os scripts a importar `compass.py`.

### 2.4 Fase 6: de ferramenta a pesquisa (18/06 e 04–07/07, V40–V77)
- **Relatórios V01–V09 (04/07, 08:21–13:44):** o primeiro relatório de projeto, em português, nasceu nessa manhã e foi revisado a cada passo. O V01 tinha 31 kB; o V02 acrescentou a tabela de parâmetros da V39; o V03, o roteiro de estudos e as melhorias propostas; o V04 documentou a V40 e o `damping_sweepV02`; o V05, exemplos de uso; o V06–V09 documentaram as V41–V44. As figuras incluídas já eram `area_vs_Q.png` e `tau_vs_Q.png`, as mesmas que a auditoria de outubro mostrou serem dominadas pela taxa.
- **V40 (04/07):** **tensor dipolar pré-calculado**, com ganho medido de 35–45× por corrida. O `damping_sweepV02.py` ganhou paralelização e retomada.
- **V41:** histerese com espaçamento log (`--hyst_spacing`, `--hyst_log_k`); pulso com `--field_delay`, `--t_pulse` e `--torque_tol`; modos `step_pos` e `step_neg`.
- **V42:** `--t_sim_full`.
- **V43:** halos também na figura de comparação.
- **V44:** painel de progresso de duas linhas.
- **V45:** vídeo de histerese com painel duplo (agulhas e curva M×H).
- **V64:** mudou só os scripts de varredura (imports). O código do simulador é igual ao da V63/V65, que são idênticas fora o número de versão; por isso não há arquivo próprio.
- **V46–V65 (05/07):** alinhamento da ajuda (V46); `--hyst_autoscale` e tradução dos comentários (V47); ajustes do HUD, mostrador e seta de direção (V49–V58); **`--damping_noise`** (V59); tradução da interface de linha de comando (**V60**), dos gráficos (V61) e dos avisos (V62); constantes pré-calculadas no laço (V63); imports e tradução dos scripts de varredura (V64–V65).
- **V66 (06/07):** barra horizontal de intensidade do campo no HUD **bidirecional**, crescendo a partir do zero, com marca vermelha em zero; tradução da ajuda de `--damping_noise`. Junto veio o Relatório **V33**, o primeiro `project_report.tex`, derivado do `relatorio_projeto_V22_EN.tex`. A comparação automática de flags não detecta mudanças visuais como essa; a fonte é o relatório de tarefa da Antigravity `df4c3b25`.
- **V67–V69:** **FORC**. Na V67, `--forc_Br`, `--forc_t_rev` e `--forc_t_sat`. Na V68, a família de curvas (`--forc_n_curves`, `--forc_Br_min`, `--forc_t_ramp_down/up`). Na V69, `--forc_sweep`.
- **V70:** **desmagnetização** (`--demag`, `--demag_cycles`, `--demag_freq`).
- **V71:** FORC com **taxa constante** (`--forc_rate`), em resposta à "armadilha da taxa variável" documentada em `FORC.md`.
- **V72:** desmagnetização rotacional com envelope linear, no lugar do decaimento geométrico 0,7ⁿ.
- **V73:** `--hyst_autoscale` ligado por padrão.
- **V74:** `--demag_temp` (opção de recozimento térmico).
- **V75:** barra de progresso com etapas ("demag" × "integrando").
- **V76:** `--demag_delay`.
- **V77 (07/07):** tempo de relógio no log final ("Completed in … sec"), documentado no Relatório **V47**.

### 2.5 Consolidação e primeira divulgação (08–09/07, V78–V79)
- **V78 (08/07 07:56), recuperada do histórico local do editor:**
  - caixa de **magnetização projetada M** no HUD dos vídeos de histerese e FORC, entre as caixas de tempo e de campo;
  - Q em notação científica abaixo de 0,1;
  - comentários do cabeçalho reorganizados em blocos (propriedades da agulha × parâmetros da simulação);
  - semente aleatória tirada do relógio em nanossegundos por padrão;
  - **contagem regressiva de 10 s** com confirmação por tecla antes de iniciar.

  O diff V77 → V78 tem +417/−158 linhas. Os guias de 08/07 (rascunho) e de 09/07 documentam a interface completa.
- **V79:** **valores padrão tirados da agulha real** (relatório de tarefa `44a7391c`): R = 7,5 mm e `needle_frac` = 2/3 (L = 10 mm), razão largura/comprimento de 0,30 (W = 3 mm), espessura de **0,26 mm**, pivô de latão (raio de 1 mm, altura de 2 mm, 8500 kg/m³) e o furo do pivô descontado da massa, de I e de m. Resultado: I ≈ 1,63×10⁻¹⁰ kg·m² e m ≈ 4,9×10⁻³ A·m². Novas flags `--pivot_radius`, `--pivot_thickness`, `--pivot_density` e `--pivot_mass`. É a origem dos 0,26 mm; a reescrita V2 (10/07) voltou ao padrão de 0,4 mm e a um espaçamento de 13 mm, e daí vem a inconsistência de espessura que ainda está pendente.
- Na Antigravity (pedidos de 08–10/07 em `docs/history/requests_antigravity_2026-07.csv`): guia detalhado da V78 em Markdown e em LaTeX, revisões do `project_report.tex`, coerência dos scripts de varredura com o novo `compass.py`, barra de progresso por worker, manual das varreduras, gráficos da análise em inglês, **FORC nas varreduras** (relatório de tarefa `4e2c5a1a`) e o caso de estudo do resumo do LAW3M.
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
Repositório criado (`ad865d7`), com auditoria de código (`docs/AUDIT.md`), parâmetros extraídos para `config.yaml`, testes de regressão e arquivamento das linhagens mortas em `archive/`. **A V2.1 virou o `compass.py` canônico** (`9a2cc63`). Ganhou ainda log adaptativo e janela de varredura lenta na histerese (`--log_adaptive`, `--hyst_slow_window`, `--hyst_slow_factor`). O `USER_GUIDE.md` foi reescrito, e o **Relatório V53** saiu na mesma data. Em **13/07**, uma tentativa de auditoria automática do código com um agente da Antigravity (`agy.py`, modelo Gemini) produziu um resumo curto (`audit_summary.json`). Os dois estão em `archive/2026-07-13_antigravity_agent_audit/`.

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
| Relatórios de julho (V01–V32; arquivos V01–V22) | `archive/OLD/relatorio_projeto_*` (git); V01–V09 também em Drive/OLD/ALLfiles (zip com figuras) |
| Relatórios V51–V53 | `archive/project_report_V51_compassV2.*`, `archive/project_report.*` (V51/V52) e o histórico do git de `docs/project_report.tex` (V53) |
| Relatório atual (V60) | `docs/project_report.tex` e `docs/project_report.pdf` |
| Transcrição das conversas de junho | `archive/OLD/compass_sim_conversation.md` (2,1 MB; fonte primária das sessões de 15–18/06, resumida na seção 5) |
| Histórico detalhado de 12/07 em diante | `git log` |
| Versões anteriores ao git, como commits datados | branch `history/pre-git` (seção 4) |
| Pedidos de junho (99) | `docs/history/requests_2026-06_to_07-05.csv` |
| Sessões da Antigravity (04–10/07) | `archive/transcripts/antigravity_2026-07/` (13 transcrições, log de conversa e relatórios de tarefa); 184 pedidos datados em `docs/history/requests_antigravity_2026-07.csv` |
| Primeiras versões dos scripts de campanha | `archive/OLD/campaign_early/` |
| Política de versões e mudanças | `CHANGELOG.md` |

---

## 4. Como navegar pelo histórico no git

A branch `history/pre-git` tem **um commit por versão arquivada** (123 commits, de 15/06 a 11/07), com a data original. O código fica sempre em `compass.py`, e o nome original do arquivo aparece na mensagem de cada commit. Os relatórios ficam em `docs/report.tex`, os guias em `docs/user_guide.tex` e os resumos em `docs/abstract.tex`. Datas aproximadas (V28, V79, V2.0 e versões salvas de novo depois da seguinte) estão marcadas na mensagem do commit.

```bash
git log --oneline history/pre-git
```

```bash
git diff code-V40 code-V41 -- compass.py
```

```bash
git show report-V22-en:docs/report.tex
```

**Tags:**
- código: `code-T1`…`code-T10`, `code-V11`…`code-V79`, `code-V2.0`, `code-V2.0.2`, `code-V2.1`;
- relatórios: `report-V01-pt` … `report-V11-pt`, `report-V12-pt`, `report-V12-en` … `report-V22-en`, `report-V47`, `report-V51`, `report-V51b`, `report-V52`;
- guias e resumos: `guide-V78-draft`, `guide-V78`, `guide-V2.0.2`, `abstract-LAW3M2026`;
- período do git: `code-v2.1.0-canonical`, `v2.2` (sua tag de julho), `report-V53`, `phase1-done`, `v2.3.0`, `report-V60`, `qs-saddle-fix`, `v2.4.0`.

Para ver tudo como uma única linha contínua (só localmente, sem alterar o GitHub):

```bash
git replace --graft ad865d7 history/pre-git
```

---

## 5. Sessões de junho (claude.ai): pedidos por sessão

Horário de Brasília. Os nomes dos arquivos de transcrição são as datas de exportação, não as das sessões.

| Sessão (arquivo) | Pedidos | Período | Conteúdo principal |
|---|---|---|---|
| compass-needle-simulation | 13 | 15/06 16:56 → 17:56 | Pedido inicial ("simular uma rede de agulhas… buscar o equilíbrio… vista por cima"), agulha em losango azul e branca, parâmetros da rede, honeycomb, campo externo, vídeo, unidades SI. |
| compass-sim-pbc-feature | 22 | 15/06 18:06 → 16/06 07:45 | Indicação do campo na imagem, numeração de vídeos, interrupção, Q, honeycomb corrigida, PBC. |
| compass-sim-gpu-perf | 21 | 16/06 07:52 → 17/06 07:02 | `t_sim` como tempo físico, modos de campo, histerese em 5 segmentos, CSV, GPU. |
| compass-sim-gpu-ctrlc | 22 | 17/06 07:12 → 13:34 | PBC no campo dipolar, registro formal de versões (V20), GPU vetorizada, Ctrl+C e Ctrl+I, inércia automática. |
| compass-sim-gpu-domains | 9 | 17/06 13:39 → 20:08 | Renderização vetorizada (V35), momento automático (V37), domínios (V38). |
| compass-sim-publication-strategy | 5 | 17/06 20:11 → 21:18 | `--make_images` (V39), primeiros scripts de campanha. |
| compass-sim-english-translation | 7 | 17/06 21:50 → 18/06 05:55 | Revisão de literatura e de arXiv, nichos de publicação, revistas. |

A lista completa, com data e hora de cada pedido, está em `docs/history/requests_2026-06_to_07-05.csv`.

---

## 6. O que ainda falta

**Não houve salto de numeração** (exceto o salto deliberado do relatório, de V53 para V60). As versões sem arquivo existiram; os arquivos é que se perderam:
- **Código V57:** o registro de versões do relatório V22 (`archive/OLD/relatorio_projeto_V22_EN.tex`, tabela de versões) descreve a sequência V56 (comprimento mínimo da seta do mostrador), **V57** (seta fina trocada por uma seta grossa preenchida, `FancyArrow`) e V58 (a mesma seta, escalonada pela intensidade do campo). A V58 substituiu o arquivo logo em seguida.
- **Relatórios V33–V46:** a tabela de histórico do relatório V47 (`git show report-V47:docs/report.tex`) tem uma linha datada para cada um, de V33 (06/07, documenta a V66) a V46 (07/07, documenta o `--demag_delay`). Perderam-se porque, desde 06/07, a Antigravity editava sempre o mesmo `project_report.tex`, e cada revisão sobrescrevia a anterior; o V47 sobreviveu só na pasta de compilação da IDE.
- **Relatórios V48–V50:** evidência indireta. O V51 foi reescrito e já não tem a tabela de histórico, mas os pedidos de 07–08/07 na Antigravity pedem para atualizar o relatório a cada versão do código (por exemplo, 08/07 07:32: "change the actual version number in the code and update the report"), o que corresponde às V78 e V79 e fecha a sequência até o V51 (10/07).

- **Código sem arquivo:** só a V57 (o mostrador com seta grossa, que é a V58 sem o escalonamento proporcional). A V64 é igual à V63/V65 no simulador. A V78 foi recuperada.
- **Relatórios sem arquivo:** V33–V46 e V48–V50. O V47 e os V01–V09 foram recuperados.
- **Já verificados (18/10):** a pasta original do Dropbox e o histórico do git dela; `~/.gemini/antigravity-ide`, que é o mesmo que `~/Dropbox/.antigravity-ide`; `~/.antigravity`, `~/.antigravity-ide` e `Application Support/Antigravity`, que só têm a instalação, as extensões e históricos de outros arquivos; o histórico local de arquivos da 'Antigravity IDE', de onde veio a V78; e `OLD/ALLfiles`, de onde vieram os relatórios V01–V09. Ainda poderiam estar no histórico de versões do próprio Dropbox (pelo site; depende do prazo de retenção do plano), no MacBook usado em julho ou num backup do NAS ou do Time Machine.
- As conversas em formato binário da Antigravity (`conversations/`, 274 MB) não foram decodificadas, mas as transcrições JSONL equivalentes (`brain/*/transcript_full.jsonl`) foram lidas e arquivadas.
- A tradução para o inglês começou em 05/07 09:34, com o pedido de criar o `compass_sim_ingV01.py` a partir da V46 (conversa `304910a1` da Antigravity).

