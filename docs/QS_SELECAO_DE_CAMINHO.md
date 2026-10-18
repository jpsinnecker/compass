# Seleção de caminho inercial nos laços quase-estáticos (Fase 2)

Data: 18/10/2026. Branch `fase2/qs-path-selection`. Script: `qs_path_selection.py`.
Dados (fora do git): `~/compass_runs/fase2_qs_paths` (tabelas, estados por passo, figuras em `figs/`).

## 1. Pergunta

No protocolo quase-estático (QS), o laço do fluxo gradiente ("od", Q₀ → 0) e o laço inercial com Q₀ = 10 coincidem
na maioria das condições. Em algumas, a área muda muito: honeycomb 10×10 (−43 %), honeycomb 6×6 (−28 %),
triangular 32×32 (−27 %), quadrada 32×32 (−15 %), honeycomb 24×24 (−10 %). A hipótese da Seção "Results III" do
relatório é que, em poucas avalanches grandes, o *overshoot* dos rotores leva o arranjo a outra bacia.
Esta nota localiza essas avalanches e testa a hipótese.

## 2. Método

1. **Alinhamento (`scan`).** As duas corridas (od e Q10, mesma semente) estão no mesmo estado em +B_max: a saturação
   apaga o estado inicial, e as sementes 0 e 1 dão o mesmo laço bit a bit (exceto honeycomb 6×6, cuja corrida Q10 s0
   foi feita no Linux da mumax). A partir daí, as linhas do CSV (uma por passo de campo, em repouso) são comparadas
   uma a uma. Com d = K (M^Q10 − M^od), em "agulhas equivalentes" (uma reversão completa muda K·M de 2):
   *início* = primeiro passo com |d| > 0,1; *divergência* = |d| > 1; *salto* = |Δd| > 0,5 em um passo.
2. **Reexecução (`rerun`).** Os pares relevantes foram rodados de novo com os argumentos exatos da campanha, mais
   duas opções novas do `compass.py`: `--qs_state_log` (ângulos de todas as agulhas em cada linha do CSV) e
   `--qs_stop_after` (para no passo 1800 nos casos grandes). As reexecuções reproduzem os CSVs da campanha
   **bit a bit** (diferença máxima de M_proj = 0), exceto honeycomb 6×6 Q10 (rodada original no Linux; ver §3).
   Também foi gravado o log de eventos (`--event_log`).
3. **Mapas (`map`).** Por passo: agulhas revertidas em uma corrida e não na outra (sinal de cos(θ−φ) diferente e
   rotação > 30°), agulhas que reverteram no passo (pegada da avalanche), borda × interior (coordenação menor que a do
   interior), reversões múltiplas da mesma agulha no mesmo passo (log de eventos) e **estabilidade do estado de
   repouso**: menor autovalor λ_min da Hessiana d²E/dθ_i dθ_j, em unidades de m·B_ref.
4. **Controle (`control`).** O mesmo fluxo gradiente, mas com um chute aleatório de ε = 10⁻³ rad em cada agulha antes
   de cada relaxação ("od com ruído"). Um estado que seja ponto de sela é então abandonado, como em qualquer arranjo
   real. A área é calculada pelo trabalho dos saltos de campo, que reproduz exatamente a área da campanha.

## 3. Resultado principal: nos honeycomb pequenos, o "efeito inercial" é um artefato de pontos de sela

Em honeycomb 6×6 e 10×10, o fluxo gradiente **para em pontos de sela**. Entre B ≈ +1,0 B_ref e o colapso do domínio
(B ≈ −0,2 B_ref), os estados de repouso od têm de 1 a 4 modos instáveis, com λ_min até −1,9 m·B_ref. O critério de
repouso do protocolo QS só testa |τ| e |ω|, e um ponto de sela tem torque nulo.

Por que o fluxo gradiente fica lá: com o campo ao longo de x, a energia é invariante pela operação
S = (espelho x → −x) ∘ (inversão de todas as agulhas). O estado saturado é simétrico, e o fluxo gradiente preserva S
até o arredondamento (|θ − S(θ)| ≤ 3×10⁻¹³ grau, do início do laço até o colapso do domínio). Quando o estado simétrico perde a
estabilidade (bifurcação em B = +1,03 B_ref no 10×10 e +0,95 B_ref no 6×6), o fluxo gradiente segue no ramo
simétrico, agora instável. A dinâmica inercial com Q₀ = 10 quebra S na bifurcação e fica em mínimos verdadeiros.

| honeycomb | passos em sela (de 3200): od | Q1 | Q3 | Q10 | od com ruído |
|---|---|---|---|---|---|
| 6×6 | 506 | 611 | 467 | 80 | 2 |
| 10×10 | 518 | 506 | 427 | 6 | 5 |

| honeycomb | área od | Q1 | Q3 | Q10 | od com ruído (3 sementes do chute) |
|---|---|---|---|---|---|
| 6×6 | 1,019 | 1,242 | 0,882 | 0,735 | 0,7337–0,7338 |
| 10×10 | 0,793 | 0,760 | 0,632 | 0,454 | 0,4556 |

As áreas são as da campanha. As contagens de sela do 6×6 Q10 vêm da reexecução no Mac, cujo caminho Q10 difere do
original rodado no Linux (área 0,730 contra 0,735). Até a corrida Q10 do 6×6 depende do arredondamento: ela também
passa 80 passos em selas.

Com o chute, o fluxo gradiente segue **o mesmo caminho** que Q₀ = 10. No 10×10, nenhuma agulha difere de Q10 em
253 passos, do início ao reencontro. A área fica a 0,3 % (10×10) e 0,1 % (6×6) da área Q10. O resultado não depende
de ε (10⁻⁴ a 10⁻² rad dão a mesma área). Q₀ = 1 e Q₀ = 3 também ficam presos em selas. Portanto a sequência
"área diminui com Q₀" nesses dois arranjos mede **quanto tempo cada dinâmica leva para sair da sela**, e não um
efeito físico da inércia. O laço quase-estático correto, para qualquer Q₀, é o laço Q10.

Como a divergência aparece (10×10, Figs. `figs/qs_paths/timeline_honeycomb_N10.pdf` e
`figs/qs_paths/map_honeycomb_N10_max_split_row1636.pdf`):
- B = +1,03 B_ref: começa uma diferença contínua de inclinação (até 55 agulhas inclinadas > 10° de forma diferente),
  sem nenhuma reversão. É a quebra de simetria da Q10.
- B = −0,03 a −0,18 B_ref: Q10 reverte em avalanches pequenas (1 a 10 agulhas por passo, quase todas no interior).
  O od continua no ramo simétrico, cada vez mais instável (4 modos negativos).
- B = −0,215 B_ref: o od reverte 36 agulhas de uma vez (7 na borda, dissipando 19,8 m·B_ref) e cai exatamente no
  estado da Q10 (d = 0). O mesmo se repete, por simetria, no ramo de volta.

No 6×6 o quadro é o mesmo (7 agulhas revertidas pela Q10 em B = −0,105 B_ref, 2 na borda; até 22 das 27 agulhas
diferentes; reencontro em B = −1,09 B_ref).

## 4. Seleção de caminho genuína: honeycomb 16×16 e 24×24

Nesses arranjos, as duas dinâmicas preservam S, e todos os estados de repouso são mínimos (λ_min > 0, com uma única
exceção marginal no od em 24×24, λ = −0,007).

- **16×16 (sem efeito na área).** A Q10 faz uma avalanche de 16 agulhas (2 na borda) em B = −0,145 B_ref. O od faz a
  **mesma** avalanche um passo depois (−0,150 B_ref). Nesse passo, o mínimo do od já estava quase instável
  (λ_min = 0,026, contra 0,42 da Q10 depois de reverter). O od com ruído também espera até −0,150. Este é o
  mecanismo da hipótese, mas em escala de um passo: o overshoot dispara uma avalanche iminente um pouco antes.
  Não há efeito na área (−0,2 %).
- **24×24 (−10 % na campanha).** Na bifurcação em B = +1,25 B_ref (λ_min → 0 nas duas corridas), od e Q10 escolhem
  ramos diferentes, ambos simétricos e estáveis. Em 40 passos, todas as 432 agulhas ficam inclinadas de forma
  diferente, sem nenhuma reversão. Perto de B = 0, as avalanches são **fileiras inteiras de 24 agulhas** (cadeias ao
  longo do campo, 2 na borda). Elas ocorrem nas duas dinâmicas, mas em ordem e posição diferentes, com até 256
  agulhas diferentes ao mesmo tempo (B = −0,14 B_ref). Áreas: od 0,387, Q10 0,348, od com ruído 0,367 (igual para 3
  sementes do chute). Portanto **cerca de metade do efeito de −10 % é inercial** (Q10 −5,3 % em relação ao od com
  ruído). A outra metade é do od sem ruído, que passa um passo em sela na bifurcação e escolhe o outro ramo.

## 5. Arranjos quadrado e triangular 32×32

- **Quadrada 32×32 (−15 % na campanha): artefato.** Em B ≈ +0,5 B_ref, o estado uniforme fica marginalmente instável
  (λ_min ≈ −0,003 a −0,014). A Q10 sai em B = +0,39 B_ref com uma reorganização contínua de inclinação: 632 agulhas
  mudam mais de 10°, sem reversão, dissipando 3,8 m·B_ref. O od fica no estado instável até B ≈ +0,10 B_ref, e daí em
  diante as avalanches (até 288 agulhas, quase todas no interior) acontecem em ordem diferente até B = +0,06 B_ref.
  O od com ruído segue a Q10 **exatamente** (0 agulhas diferentes de Q10 de B = +0,40 a +0,06). Área no laço
  completo: od 0,2257, Q10 0,1912, od com ruído 0,1897. Diferença Q10 × od com ruído: +0,8 %.
- **Triangular 32×32 (−27 % na campanha): misto, e sensível.** Não há selas relevantes (λ_min ≥ −0,006). A primeira
  divergência é **inercial**: em B = +0,04 B_ref, as duas dinâmicas nucleiam domínios revertidos nos mesmos dois
  cantos (a uma ou duas fileiras da borda), mas a Q10 reverte 46 agulhas e o od 40. A avalanche inercial é um pouco
  maior, e uma faixa diagonal de agulhas inclinadas de forma diferente liga os dois núcleos
  (Fig. `figs/qs_paths/map_triangular_N32_first_split_row1592.pdf`). O od com ruído segue o od nesse passo, ou seja,
  é a inércia que muda o resultado. Depois, perto de B = 0, há avalanches de até 286 agulhas em sequências diferentes,
  com até 212 agulhas diferentes. O od com ruído se afasta das duas corridas e, em B = −0,35 B_ref, faz a mesma
  avalanche de 150 agulhas que a Q10 (50 na borda), que o od só faz dois passos depois. Áreas: od 0,437, Q10 0,317,
  od com ruído 0,389. Portanto cerca de 11 pontos dos −27 % vêm da referência sem ruído, e cerca de −18 % (Q10 em
  relação ao od com ruído) seriam inerciais. Como o caminho depende do chute, esse número ainda precisa da dispersão
  entre sementes do chute (em execução; ver §8).

### Resumo dos casos

| caso | efeito Q10 na campanha (vs od) | Q10 vs od com ruído | mecanismo |
|---|---|---|---|
| honeycomb 6×6 | −28 % | +0,1 % | od preso em selas simétricas (artefato) |
| honeycomb 10×10 | −43 % | −0,3 % | od preso em selas simétricas (artefato) |
| honeycomb 16×16 | −0,2 % | −0,2 % | avalanche de 16 agulhas um passo antes (inercial, sem efeito) |
| honeycomb 24×24 | −10 % | −5,3 % | ramo diferente numa bifurcação, depois fileiras de 24 em outra ordem |
| quadrada 32×32 | −15 % | +0,8 % | od preso num estado marginalmente instável (artefato) |
| triangular 32×32 | −27 % | −18 % (1 semente do chute) | núcleos de canto maiores na Q10 (46 × 40), depois caminho sensível |

## 6. Que tipo de avalanche dispara a divergência

- **Não é uma avalanche grande.** A divergência começa com uma diferença contínua de inclinação, na bifurcação de um
  estado simétrico (honeycomb) ou perto dela, sem reversões. As primeiras reversões diferentes são pequenas
  (1 a 46 agulhas) e quase todas no interior.
- **A avalanche grande vem depois e é do od.** Nos casos de artefato, o od acumula a instabilidade e reverte de uma vez
  (36 de 75 agulhas no 10×10), caindo no estado da Q10.
- **Não há reversões múltiplas por overshoot.** Nos logs de eventos das reexecuções, a reversão dupla de uma agulha no
  mesmo passo de campo é rara e não é própria da inércia: 1 agulha no honeycomb 6×6 Q10, 2 no triangular 32×32 **od**,
  nenhuma nos demais casos. O efeito inercial visível é outro: a avalanche da Q10 é um pouco maior (46 × 40 agulhas no
  triangular 32×32) ou acontece um passo antes (honeycomb 16×16). As avalanches com s > K da
  Fig. QS3 do relatório estão **todas no primeiro passo de campo**: são a relaxação do estado inicial aleatório, não
  avalanches do laço. Dentro do laço, o maior s é 204 em todas as dinâmicas.
- **Borda × interior:** nas avalanches que diferem entre od e Q10, a fração de agulhas de borda não passa da fração
  de sítios de borda (quadrada 32×32: 6 % contra 12 % dos sítios; triangular 32×32: 9 % contra 12 %;
  honeycomb 10×10: 21 % contra 31 %; honeycomb 24×24: 19 % contra 13 %). As divergências são eventos de interior.
  A nucleação do triangular 32×32 é perto dos cantos, mas a uma ou duas fileiras da borda. Há alguns eventos tardios
  só de borda (por exemplo, 30 agulhas de borda na quadrada 32×32 em B = −0,885 B_ref, um passo antes na Q10).
- **Tamanho:** a primeira reversão diferente tem de 1 a 46 agulhas (de 0,1 % a 5 % de K). As avalanches de
  centenas de agulhas vêm depois, quando os caminhos já se separaram.
- **Divergências são comuns, efeitos na área não.** od e Q10 divergem em 13 das 15 combinações geometria/tamanho
  (todas, exceto triangular 6×6 e 10×10), mas quase sempre se reencontram em poucos passos.

## 7. Consequências

1. **Relatório (Results III).** A frase "Q₀ = 10 muda a área em 10–43 %" não se sustenta como efeito inercial nos
   honeycomb 6×6 e 10×10 nem na quadrada 32×32: lá, o fluxo gradiente de referência está errado. No honeycomb 24×24
   e no triangular 32×32, parte do efeito é inercial (cerca de −5 % e −18 % em relação ao od com ruído). O mesmo
   problema afeta Q₀ = 1 e 3 (a comparação "Q₀ = 1 é fluxo gradiente" dos honeycomb pequenos compara duas dinâmicas presas em selas). A
   legenda da Fig. QS3 ("s > K, assinatura direta de overshoot") deve ser corrigida: esses eventos são a relaxação
   inicial.
2. **Protocolo QS.** O critério de repouso precisa testar estabilidade: por exemplo, um chute de 10⁻³ rad seguido de
   nova relaxação antes de registrar o estado, ou o sinal de λ_min. Depois disso, a grade od/Q1/Q3 precisa ser
   recalculada antes de qualquer afirmação sobre o efeito de Q₀ no laço quase-estático.
3. **Fase 3 (memória).** A seleção de caminho genuína existe (16×16, 24×24), mas é decidida perto de bifurcações e
   instabilidades, onde perturbações de 10⁻³ rad também decidem o caminho. Para estudar memória, é preciso tratar
   isso estatisticamente: conjuntos de perturbações (ruído ou desordem fraca de amortecimento) para cada Q₀.

## 8. Reprodução

```bash
python3 qs_path_selection.py scan
python3 qs_path_selection.py rerun --cases honeycomb_N06 honeycomb_N10 honeycomb_N16 --jobs 6
python3 qs_path_selection.py rerun --cases honeycomb_N24 square_N32 triangular_N32 --stop_after 1800 --jobs 6
python3 qs_path_selection.py rerun --cases honeycomb_N06 honeycomb_N10 --labels Q1 Q3
python3 qs_path_selection.py map --cases honeycomb_N06 honeycomb_N10 honeycomb_N16 honeycomb_N24
OPENBLAS_NUM_THREADS=3 python3 qs_path_selection.py map --cases square_N32 triangular_N32
python3 qs_path_selection.py control --cases honeycomb_N06 honeycomb_N10 honeycomb_N16 --kick_seeds 0 1 2
python3 qs_path_selection.py control --cases honeycomb_N24 --stride 5 --kick_seeds 0 1 2
OPENBLAS_NUM_THREADS=2 python3 qs_path_selection.py control --cases square_N32 triangular_N32 --stride 50
```

O controle de laço completo em N = 32 leva cerca de 3 h por semente do chute (K = 1024). As sementes 1 e 2 do chute
no triangular 32×32 estavam rodando quando esta nota foi escrita
(`control_kicked_loop_tri32k{1,2}.csv`).

Saídas: `path_selection_scan.csv` (todas as combinações), `path_selection_jumps.csv`, `path_selection_spatial.csv`
(eventos-chave), `timeline_<caso>.csv` (por passo, com λ_min), `control_noisy_od_*.csv` (janela),
`control_kicked_loop*.csv` (laço completo), figuras `figs/map_*.png` (od, Q10 e |Δθ| no passo) e
`figs/timeline_*.png`. Cópias em PDF das figuras principais estão em `docs/figs/qs_paths/`.
