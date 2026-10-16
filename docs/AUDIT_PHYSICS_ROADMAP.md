# Auditoria do projeto "simula agulhas" e roteiro para artigo

Data: 16/10/2026. A auditoria foi feita primeiro na cópia do Google Drive (`compassV2_1.py`, de 11/07) e depois rechecada contra o `compass.py` canônico deste repositório (commit `cb14300`, de 13/07).
Os achados de física das seções 2.1 a 2.4 continuam válidos no `compass.py` canônico: o código de cutoff, dt, Q e B_eff é o mesmo.
Este documento complementa `docs/AUDIT.md`, que trata da organização do código; aqui o foco é a validade física e o roteiro para o artigo.
Os testes que codificam estes achados estão em `tests/test_physics.py`. Os marcados `xfail(strict=True)` são defeitos conhecidos: quando um deles for corrigido, o teste passa a falhar como XPASS, e é preciso remover o marcador.

---

## 1. Diagnóstico em uma frase

O código está bem construído (integrador correto, metadados, contadores de avalanche robustos). Porém **os resultados preliminares (área do laço × Q, expoente τ × Q) ainda não medem o que dizem medir**: estão dominados pela taxa de varredura do campo, por um cutoff que corta ~20 % do campo dipolar e por uma definição de Q que muda com o protocolo. A pergunta científica também precisa ser reenquadrada para sobreviver a um revisor.

---

## 2. Problemas encontrados (em ordem de gravidade)

### 2.1 Taxa de varredura não controlada: o efeito de Q na área do laço é sobretudo atraso viscoso (CRÍTICO)
Teste: rede quadrada 12×12, sem cutoff, mesma semente, variando apenas a duração da rampa:

| Q | rampa = 71 T0 | 283 T0 | 1134 T0 |
|---|---|---|---|
| 0,06 | 16,1 mT | 6,9 mT | 2,8 mT |
| 0,6 | 4,0 mT | 1,5 mT | 0,68 mT |
| 6 | 1,27 mT | 0,67 mT | 0,39 mT |

- A área cai aproximadamente como (taxa)^0,6 em **todos** os Q e não converge para um laço quase-estático.
- No regime superamortecido, o tempo de relaxação é τ ≈ b/(mB), que cresce ∝ 1/Q. Então "área maior a Q baixo" é, em boa parte, histerese dinâmica trivial (o sistema não acompanha o campo). O pico em Q≈0,1 e a queda em Q=0,05 da figura `figs/area_vs_Q.png` têm a forma típica de laço de relaxação (tipo Debye), e não de física coletiva.
- **Consequência:** a figura área×Q e a frase do resumo do LAW3M ("damping changes the shape of the hysteresis loop") não se sustentam como estão.
- **Correção:** tratar a taxa adimensional r = (dB/dt)·τ(Q)/B_ref como eixo explícito; medir A(r, Q) e extrair o limite quase-estático A₀(Q) = lim_{r→0} A. Só a diferença entre A₀(Q) para Q diferentes é efeito inercial genuíno.

### 2.2 Cutoff de 3,5 r_nn: erro grande e sem ganho de velocidade (ALTO)
Campo no centro de um array 36×36 uniformemente magnetizado, em unidades de B_ref:

| cutoff | 1 r_nn | 2,5 | 3,5 (padrão) | 5 | 10 | soma completa |
|---|---|---|---|---|---|---|
| B/B_ref | 1,00 | 1,66 | 1,80 | 1,92 | 2,10 | 2,18 |

- Em 2D, a soma de 1/r³ converge, mas devagar (a cauda vale ~π m/(a² r_c)). O cutoff padrão subestima o campo interno em cerca de 17 % e **remove os efeitos de forma e de borda**, que são justamente a física de um array finito.
- O tensor é uma matriz densa K×K. Zerar elementos não acelera o produto matriz-vetor. **O cutoff tem custo zero para remover:** use a soma completa por padrão.
- PBC com 1 imagem e cutoff é inconsistente. Para o artigo, use contornos livres (o objeto físico é finito) e soma completa. Ewald só seria necessário se o artigo precisar de bulk, o que provavelmente não é o caso.

### 2.3 Passo de tempo marginal (MÉDIO)
- O dt é definido a partir de ω0 calculado com B_ref (um único vizinho). O campo local real é cerca de 2× maior, então ω_loc·dt ≈ 0,47 com `dt_factor=0.035`. O monitor de estabilidade sinalizou **100 % dos passos**.
- Teste em b=0 e campo zero: com dt_factor 0,035 o erro de energia chega a 1,4 % de |E_dip|; com dt_factor 0,01 cai para 0,12 % (escala ∝ dt², ou seja, o Verlet funciona).
- **Correção:** definir dt a partir de max ω_loc (ou de B_ref multiplicado pelo fator de rede), com dt·ω_loc ≤ 0,1.

### 2.4 Q depende do protocolo (MÉDIO)
- Hoje Q = ω0 I/b com ω0 = √(m·max(B_ext, B_ref)/I). Com `B_max_factor 8`, o "Q" reportado é ≈ √8 ≈ 2,8× maior que o Q intrínseco da rede, e muda entre histerese, FORC e pulso.
- **Correção:** definir Q₀ = √(m B_loc I)/b, com B_loc igual ao campo dipolar local do estado fundamental (propriedade da rede e não do protocolo), e reportar o Q de campo separadamente.
- **Corrigido em out/2026:** o metadado `Q` agora é Q₀ = √(m B_lat I)/b, onde B_lat é a mediana de |B_i| nos estados uniformes (≈ 2,2 B_ref em redes grandes). Também são gravados `Q_ref` (com B_ref) e `Q_field` (a definição antiga).
- Achado adicional: o `damping_sweepV03.py` rotulava a grade com Q_ref, enquanto o motor gravava Q_field. O eixo de `figs/area_vs_Q.png` e `figs/tau_vs_Q.png` é, portanto, Q_ref. Para a mesma rede, Q₀ ≈ 1,2–1,5 × Q_ref, dependendo do tamanho.

### 2.5 Estatística de avalanches insuficiente e pipeline desatualizado (ALTO para o artigo)
- `figs/tau_vs_Q.png`: n_tail entre 30 e 40 eventos. Nenhum expoente com essa cauda é defensável, e as barras de erro não incluem a escolha de x_min.
- `damping_sweep_analysis.py` detecta avalanches por limiar MAD em dM/dB e **não usa** os contadores robustos (Schmitt + dwell) nem o `--event_log` da V2.1.
- As figuras `figs/area_vs_Q.png` e `figs/tau_vs_Q.png` foram geradas antes da consolidação de 13/07, quando `import compass` ainda apontava para a V2.0.2, ou seja, sem o contador corrigido. No repositório atual o import já aponta para o motor canônico, mas os dados precisam ser regenerados.
- **Recomendação:** retirar "expoentes/universalidade" do foco principal. Avalanches entram como observável secundário (distribuição de tamanhos e cutoff × Q), nunca como afirmação de classe de universalidade.

### 2.6 Física do modelo (para a robustez das conclusões)
- **Dipolo pontual:** L/a = 10/13 ≈ 0,77, então o campo de campo próximo difere bastante do dipolar. Um modelo de halteres (cargas ±q nas pontas) é necessário ao menos como teste de robustez, porque pode mudar até o estado fundamental da rede quadrada.
- **T = 0 e determinístico:** a única aleatoriedade são os ângulos iniciais e 5 % no amortecimento. Não há limiar de chaveamento intrínseco, então não há "histerons" no sentido de Preisach.
- **Atrito de Coulomb no pivô** é físico no aparato real e cria limiares estáticos (pinning). Ele é a forma natural de introduzir desordem controlada e histerons verdadeiros. Recomendo implementar como opção (com regularização suave de sgn(ω)).
- **Balanço de energia:** falta o livro-caixa de trabalho do campo e de dissipação (W = −∫ m·dB, D = ∫ b ω² dt). A identidade "área do laço = energia dissipada por ciclo" é um teste de validação excelente e uma figura de artigo. **Corrigido em out/2026:** colunas `W_field_J` e `E_diss_J` no CSV e `energy_ledger` nos metadados. O balanço fecha em O(dt²), e a identidade área = trabalho do campo é testada.

### 2.7 Organização
Grande parte disso já foi resolvida no repositório em 13/07 (motor único, testes de regressão, documentação reescrita; ver `docs/AUDIT.md`). Pendências:
- **A pasta do Google Drive está defasada** em relação ao GitHub. O trabalho deve ser feito só no clone (fora do Drive); o Drive serve apenas como arquivo.
- A espessura padrão ainda é inconsistente: 0,40 mm em `physics.compass_engine` do `config.yaml`, contra 0,26 mm em `damping_sweep`, `calc_defaults_temp` e no relatório. É preciso decidir pela agulha real.
- O repositório é **público** e contém PDFs de editoras em `papers/` (incluindo um livro baixado do Anna's Archive), além de resultados não publicados. Recomenda-se torná-lo privado ou remover `papers/` do histórico.
- Os testes existentes são de regressão: congelam os números atuais, inclusive os efeitos do cutoff e do dt. Ao corrigir a física, os arquivos de referência devem ser regenerados de propósito (`python3 tests/generate_reference.py`).

### 2.8 Sobre as críticas em `Code Suggestions.md` e `Compass Project.md`
Parte delas está **errada** e não deve guiar o trabalho:
- "Trocar Verlet por RK4": não deve ser feito. O Verlet com amortecimento implícito usado é simplético em b=0 e tem erro ∝ dt² (verificado acima). RK4 não é simplético e deriva em energia em integrações longas. O problema real é o tamanho do dt (seção 2.3).
- "Precisa de Ewald": só vale se usar PBC. A solução certa aqui é contorno livre com soma completa.
- "Flip counter infla avalanches por ringing": isso já foi corrigido na V2.1.

Outra parte está **correta** e é decisiva:
- O enquadramento como "análogo de spin ice nanomagnético" é vulnerável (LLG ≠ Newton).
- Dipolo pontual com L/a ≈ 0,77 é problemático.
- 36×36 não sustenta expoentes.
- FORC de histerons idênticos é sensível à taxa.

---

## 3. Onde estão as lacunas na literatura (avaliação honesta)

O que já existe:
- Agulhas de bússola como análogo dipolar (Olive & Molho, PRB 1998; Novak & Sinnecker, JMMM 2004).
- Rotores magnéticos macroscópicos com inércia e dissipação explícitas: Mellado, Concha & Mahadevan, *Macroscopic magnetic frustration*, PRL 109, 257203 (2012). Esse trabalho já argumenta que spin ice microscópico não captura inércia, então **"inércia como eixo novo" sozinho não é novidade suficiente**.
- Histerese projetada com cadeias dipolares de rotores, incluindo a escala de tempo de rotor que limita a taxa (Concha et al., arXiv:1710.08391).
- FORC em sistemas correlacionados (Ruta et al., Sci. Rep. 2017).
- Leis de escala de histerese dinâmica A − A₀ ∝ R^β: teoria de campo médio dá β=2/3, e ruído reduz β (PRE 108, 024101, 2023; arXiv:2603.24007, 2026). Nesses trabalhos, porém, não há inércia nem interação dipolar de longo alcance.
- Memória e redes de histerons em metamateriais mecânicos (Keim, Paulsen, Nagel, van Hecke, Mungan e outros; Sirote-Katz et al., Nat. Commun. 2024).

Lacunas que o projeto pode ocupar, em ordem de potencial:

**A. Memória de ponto de retorno (RPM) e inércia (maior potencial; alvo PRL/PRE/PNAS).**
O teorema de Sethna garante RPM para dinâmica superamortecida com interações "ferromagnéticas" (propriedade de no-passing). Interações dipolares frustradas **e** inércia violam as hipóteses do teorema por razões diferentes: o overshoot inercial pode atravessar barreiras que a dinâmica adiabática não atravessaria.
Hipótese testável: a fidelidade de RPM cai à medida que Q cresce, e o limiar depende da geometria (quadrada × triangular × honeycomb).
Por que é bom: a pergunta é nítida, se conecta à literatura ativa de "memory formation in matter", aproveita o protocolo FORC já implementado (subloops são o teste natural de RPM) e é verificável na bancada.

**B. Histerese dinâmica com crossover inercial (sólido; alvo PRE/PRB/PR Research).**
Mapa A(r, Q) para as 3 redes, com extração de A₀(Q) e do expoente de taxa β(Q). Pergunta: o colapso de dados A−A₀ ∝ r^β sobrevive ao crossover superamortecido→subamortecido? E com interação de longo alcance?
É o artigo que transforma o "problema 2.1" em resultado. É de baixo risco, porque o código quase já faz isso.

**C. Validação experimental (multiplica o valor de A e B).**
O autor tem o aparato de 2004. Medir Q real (decaimento livre de uma agulha), atrito de Coulomb e 2 ou 3 laços em taxas diferentes, e comparar com a simulação calibrada. Simulação mais bancada é o que leva de JMMM para PRE/PR Applied/PNAS.

**D. Lacunas menores, que entram como seções e não como artigo:**
dipolo pontual × halteres no estado fundamental e no campo coercivo; FORC "de histerons idênticos" como impressão digital de correlação (complementa Ruta 2017); distribuição de avalanches como observável secundário.

**O que abandonar como afirmação central:** "classe de universalidade de avalanches dependente de Q" e "análogo de spin ice artificial nanomagnético".

---

## 4. Devo usar agentes e configurar skills antes de mexer no código?

**Sim, mas pouco, e só depois da higiene básica.** Em ordem:

1. **Antes de tudo (sem agente):** trabalhar só no clone git fora do Drive (o motor único `compass.py` já existe desde 13/07) e manter um `CLAUDE.md` com as convenções físicas: unidades SI, definição de Q₀, dt por ω_loc, contorno livre com soma completa, e "nunca afirmar universalidade". Esse arquivo é o que mais melhora o trabalho de qualquer agente, porque evita que cada sessão redescubra as convenções.
2. **Duas skills de projeto valem a pena:**
   - `validate-run`: roda a bateria de validação (torque tensor × par a par, conservação de energia em b=0, identidade área = dissipação, convergência em dt, tamanho e taxa) e emite um `quality_control.json`. Será usada dezenas de vezes.
   - `campaign`: padroniza como se lança uma varredura (grade em Q₀ e r, sementes, manifesto, retomada) e onde ficam os dados (fora do Google Drive).
3. **Agentes (subagentes) são úteis em dois pontos:**
   - revisão de literatura, via a skill `deep-research`, sobre RPM/no-passing com inércia, histerese dinâmica e rotores magnéticos, para confirmar que a lacuna A está aberta antes de investir nela;
   - um agente "revisor hostil" ao fim de cada fase, lendo o rascunho como referee de PRE.

   Para rodar simulações, use scripts com multiprocessing e não agentes: agentes não aceleram computação.
4. **Não recomendo** já uma arquitetura multiagente complexa. O gargalo hoje é de definição física, não de mão de obra.

---

## 5. Roteiro sugerido

**Fase 0: higiene (≈ 1 semana)**
- Clone fora do Drive, `CLAUDE.md`, dados de campanha fora do git e do Drive. (Git e motor único já existem.)
- `tests/` com pytest: torque tensor × direto; energia em b=0; simetria; reprodutibilidade por semente.

**Fase 1: correções numéricas (≈ 2 semanas)**
- Soma completa por padrão; cutoff só opcional.
- dt por ω_loc; Q₀ intrínseco; taxa adimensional r como parâmetro de entrada (em vez de `t_sim`).
- Livro-caixa de trabalho e dissipação; teste da identidade área = dissipação.
- Wrapper de varredura apontando para a V2.1 e análise usando os contadores e o `event_log`.
- Espessura padrão unificada (0,26 ou 0,40 mm: decidir pela agulha real).

**Fase 2: estudo piloto da histerese dinâmica (≈ 3 semanas), que dá o artigo B**
- Grade Q₀ ∈ [0,03, 30] (log, 10 pontos) × r ∈ 4 décadas × 3 geometrias × N ∈ {12, 20, 30, 40} × 5 sementes.
- Produtos: A(r, Q₀), A₀(Q₀), β(Q₀), H_c e M_r quase-estáticos.

**Fase 3: memória e RPM (≈ 4 semanas), que dá o artigo A**
- Protocolo de subloops aninhados e FORC na taxa quase-estática.
- Métrica: distância de estado no retorno ‖θ_retorno − θ_ref‖ e fração de sítios que violam RPM, em função de Q₀ e da geometria.
- Controle: mesma análise com dinâmica puramente superamortecida (Q → 0) e com Coulomb ligado.

**Fase 4: robustez (≈ 2 semanas)**
- Modelo de halteres (campo próximo exato até 2–3 vizinhos, dipolo além disso).
- Atrito de Coulomb como desordem; tamanho finito; forma da amostra.

**Fase 5: bancada (em paralelo às fases 2–4)**
- Medir I, m, b e o atrito estático de uma agulha (decaimento livre); 3×3 a 6×6 agulhas; laços em 2 ou 3 taxas; teste de RPM com subloop.

**Fase 6: redação**
- Artigo A (se a Fase 3 der resultado positivo): PRL ou PNAS, com PRE como alternativa segura.
- Artigo B: PRE ou PR Research. O resumo do LAW3M deve ser ajustado para não afirmar o efeito de Q na área antes da Fase 2.

---

## 6. Ponto de decisão

Antes da Fase 2, decidir o foco:
- **A (RPM/memória):** mais impacto e mais risco.
- **B (histerese dinâmica):** mais seguro; muito do código já existe.

Recomendação: executar B primeiro, porque ele é pré-requisito técnico de A (é preciso conhecer o limite quase-estático), e decidir sobre A com os dados de B em mãos.

Fontes consultadas:
- [Mellado, Concha & Mahadevan, Macroscopic magnetic frustration (arXiv:1207.2374)](https://arxiv.org/abs/1207.2374)
- [Designing hysteresis with dipolar chains (arXiv:1710.08391)](https://arxiv.org/pdf/1710.08391)
- [Dynamic hysteresis at a noisy saddle-node (PRE 108, 024101)](https://dx.doi.org/10.1103/physreve.108.024101)
- [Universal scaling laws for dynamical-thermal hysteresis (arXiv:2603.24007)](https://arxiv.org/pdf/2603.24007)
- [Emergent disorder and mechanical memory in periodic metamaterials (Nat. Commun. 2024)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11109184/)
- [Dynamics of a macroscopic artificial spin ice (APS March Meeting 2024)](https://meetings-archive.aps.org/mar/2024/k22/8)
