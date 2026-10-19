# Protocolo de bancada: I, m, b e atrito de uma agulha real

Data: 19/10/2026. Fase 5 do roteiro (`docs/AUDIT_PHYSICS_ROADMAP.md`). Análise: `bench_fit.py` (testes em `tests/test_bench_fit.py`).

## 0. Modelo e objetivo

Uma agulha num campo horizontal uniforme B obedece a

  I θ̈ = −m B sin(θ − θ₀) − b θ̇ − τ_c sign(θ̇),

e fica parada enquanto |m B sin(θ − θ₀)| ≤ τ_s.

| símbolo | grandeza | uso na simulação |
|---|---|---|
| I | momento de inércia (kg·m²) | `--inertia` ou geometria |
| m | momento magnético (A·m²) | `--moment` ou geometria |
| b | amortecimento viscoso (N·m·s) | `--damping`; dá Q₀ = √(m B_lat I)/b |
| τ_c | atrito dinâmico (seco, de Coulomb) (N·m) | ainda não está no modelo (Fase 4) |
| τ_s | atrito estático (N·m) | ainda não está no modelo (Fase 4) |

O simulador assume τ_c = τ_s = 0. A grandeza que decide se isso é aceitável é **τ_s/(m B_ref)**. Se ela passar de alguns por cento, o atrito seco funciona como limiar de chaveamento nas avalanches e precisa entrar no modelo.

**Ordens de grandeza** (agulha da V79: L = 10 mm, W = 3 mm, 0,26 mm, pivô de latão; I ≈ 1,6×10⁻¹⁰ kg·m², m ≈ 4,9×10⁻³ A·m²):

| campo | f₀ = √(mB/I)/2π |
|---|---|
| Terra, componente horizontal (~20 µT) | ≈ 4 Hz |
| 100 µT | ≈ 9 Hz |
| B_ref da rede com r_nn = 15 mm (≈ 0,3 mT) | ≈ 15 Hz |

## 1. Material

- **Bobinas de Helmholtz** com campo horizontal de 0 a ~1 mT, calibradas com gaussímetro ou sensor Hall no centro. A região uniforme (±1 %) precisa conter a agulha inteira.
- **Campo da Terra:** alinhe o eixo das bobinas com o norte magnético. A medida 3 dá a componente residual; não é preciso cancelá-la.
- **Registro do ângulo:** câmera vista de cima, a pelo menos 20 quadros por período (240 quadros/s bastam até ~12 Hz; acima disso, use 480–960). Iluminação contínua (LED de corrente contínua), porque luz que pisca a 100/120 Hz produz batimento no vídeo. Extraia θ(t) com o Tracker (física, gratuito) ou com OpenCV, marcando as duas pontas da agulha.
  - Não cole espelho nem marcador na agulha: isso muda I.
- **Balança** de 0,1 mg (ou 0,01 mg) e paquímetro/micrômetro.
- **Sensor Hall ou magnetômetro** de 3 eixos para a medida do dipolo (o do celular serve para ordem de grandeza).
- Opcional: campânula com bomba de vácuo simples, para separar o arrasto do ar do atrito do pivô.

**Formato dos arquivos para o `bench_fit.py`:** CSV com cabeçalho.

| medida | colunas |
|---|---|
| frequência × campo | `B_T,f_Hz` |
| dipolo | `r_m,B_T` |
| decaimento livre | `t,theta` (radianos, ou graus com `--deg`) |
| banda morta | `theta_final` |

## 2. Momento de inércia I pela geometria

1. Pese a lâmina sem o pivô (M) e, separadamente, o pivô (m_p). Meça o comprimento L de ponta a ponta, a largura máxima W e o raio do pivô r_p.
2. Calcule:

```bash
python3 bench_fit.py geometry --mass 3.1e-5 --length 0.010 --width 0.003 --pivot_mass 5.3e-5 --pivot_radius 0.001
```

O script devolve I = M(L² + W²)/24 + ½ m_p r_p² e a **espessura efetiva** M/(ρ·LW/2). Ela resolve a pendência de 0,26 ou 0,4 mm. Com os números do exemplo (aço, 0,26 mm), I ≈ 1,67×10⁻¹⁰ kg·m².

## 3. Razão m/I e campo residual: frequência em função do campo

1. Para 6 a 8 valores de corrente (B_bobina de ~20 µT a ~0,5 mT), solte a agulha de **no máximo 5°** e registre pelo menos 10 períodos.
2. Meça a frequência pelo tempo entre cruzamentos do zero.
   - Se o amortecimento for grande (menos de ~20 oscilações até parar), corrija: ω₀² = ω_d² + (b/2I)², com b da medida 6.
3. Ajuste:

```bash
python3 bench_fit.py freq freq.csv --I 1.63e-10
```

A relação ω₀² = (m/I)(B_bobina + B_resid) dá **m/I** (inclinação) e **B_resid**, o campo residual ao longo do eixo (intercepto). Se B_resid diferir muito da componente horizontal local da Terra, há outro campo ou ímã por perto. Os resíduos da reta (`rel_rms`) devem ficar abaixo de ~1 %.

## 4. Momento de inércia I por massa adicionada (independente da geometria)

1. Num campo fixo, meça f₀.
2. Prenda simetricamente uma inércia conhecida ΔI. Por exemplo, dois pedaços iguais de fio de cobre de massa μ cada, a uma distância r do eixo, dão ΔI = 2μr². Escolha ΔI ≈ 0,3–1 I, para que a frequência mude de 15–40 %.
3. Meça f₁ no mesmo campo e calcule:

```bash
python3 bench_fit.py addmass --f0 8.73 --f1 7.38 --dI 6.5e-11
```

A fórmula é I = ΔI · f₁²/(f₀² − f₁²). Ela não depende de m nem de B.

Compare com a medida 2. Se as duas concordarem em ~5 %, use a média e, a partir dela, **m = (m/I) × I**.

## 5. Momento magnético m de forma independente (checagem)

1. Retire a agulha do pivô e fixe-a. Meça com o sensor Hall o campo **sobre o eixo dela**, a distâncias r de 3 a 8 cm do centro (r ≥ 3L, para valer a aproximação de dipolo). Desligue as bobinas e subtraia o fundo, ou deixe o ajuste estimá-lo.
2. Ajuste:

```bash
python3 bench_fit.py dipole dipolo.csv
```

A relação usada é B = (μ₀/4π)·2m/r³ + B_fundo. O valor de m deve concordar com o da medida 4 em ~10 %. A incerteza aqui vem sobretudo da posição do sensor (o elemento sensível fica alguns milímetros dentro do encapsulamento): meça r até o elemento, não até a capa.

## 6. Amortecimento viscoso b e atrito dinâmico τ_c: decaimento livre

1. Em 3 ou 4 campos (por exemplo 50, 100, 200 e 400 µT), solte a agulha de 30–60° e grave até ela parar. Repita 3 vezes em cada campo.
2. Exporte `t,theta` e ajuste:

```bash
python3 bench_fit.py decay solta_100uT.csv --I 1.63e-10 --m 4.9e-3 --B 1.18e-4 --deg
```

Em `--B` vai o campo **total** na agulha (B_bobina + B_resid da medida 3).

O script faz duas coisas:
- **Ajuste da envoltória:** usa os picos sucessivos com o mapa A_{k+1} = e^(−bT/4I)·A_k − 2τ_c/(mB). Ele separa as duas assinaturas: o atrito viscoso tira uma **fração** fixa da amplitude por oscilação, e o atrito seco tira uma **quantidade** fixa (4τ_c/(mB) por período). É aproximado em amplitudes grandes e serve de chute inicial.
- **Ajuste da equação completa:** integra o movimento e ajusta θ(t) inteiro até a agulha parar. Parâmetros livres: b, τ_c, θ₀ e as condições iniciais.

Ele devolve b e τ_c com incertezas, o Q de uma agulha isolada, √(mBI)/b, e τ_c/(mB).

**Checagens:**
- b deve sair **igual** em todos os campos. Se crescer com B, há perda que o modelo não descreve (correntes parasitas no pivô, por exemplo).
- τ_c também deve ser independente de B. A perda angular por período, 4τ_c/(mB), cai como 1/B, e é isso que permite separar os dois termos.
- No vácuo, b deve cair (a parte do ar) e τ_c ficar igual (o pivô).

Nos testes com dados sintéticos de ruído de 0,2° e 240 quadros/s, o ajuste completo recupera b e τ_c em ~2 %.

## 7. Atrito estático τ_s

Dois métodos complementares.

**a) Banda morta.**
1. Num campo fixo, solte a agulha 20 a 30 vezes, de ângulos e lados diferentes, e anote o ângulo final de repouso (no vídeo ou com transferidor).
2. Os ângulos finais preenchem a faixa |θ − θ₀| ≤ arcsin(τ_s/(mB)). Ajuste:

```bash
python3 bench_fit.py deadband repouso.csv --m 4.9e-3 --B 1.18e-4 --deg
```

O resultado é uma **cota inferior** de τ_s, porque a agulha raramente para exatamente na borda da faixa. Use campos baixos (20–50 µT), em que a faixa é mais larga e fácil de medir.

**b) Rotação lenta do campo.**
1. Gire a direção do campo devagar (≲ 1°/s), com duas bobinas ortogonais ou girando a base da agulha num prato giratório.
2. Meça o ângulo entre a agulha e o campo no instante em que ela "descola". Então τ_s = mB·sin θ_atraso.
3. Repita nos dois sentidos de rotação e em dois campos. τ_s deve sair independente de B. Se não sair, o problema é alinhamento (bobinas não ortogonais ou pivô fora do centro), não atrito.

## 8. O que levar para a simulação

1. **Q₀ real** da rede montada: Q₀ = √(m B_lat I)/b, com B_lat calculado pelo `compass.py` para a geometria e o espaçamento da rede (campo `B_lat_T` nos metadados de qualquer corrida). Com isso dá para saber em que regime da Fig. QS2 do relatório o experimento está.
2. **τ_s/(m B_ref)** e **τ_c/(m B_ref)**: se forem ≲ 1 %, o modelo sem atrito seco é adequado para comparação. Se forem maiores, a Fase 4 (atrito de Coulomb no `compass.py`) vem antes da comparação.
3. **Variação entre agulhas:** repita as medidas 3, 6 e 7 em 5 a 10 agulhas. A dispersão de m, I, b e τ_s é a desordem real do sistema, e alimenta diretamente a "desordem fraca" prevista para a Fase 3.
4. **Espessura:** a espessura efetiva da medida 2 fixa o valor a usar nas simulações (0,26 ou 0,4 mm).

## 9. Registro

Para cada agulha, guarde:
- identificação, massa, L, W, pivô, espessura efetiva;
- os CSVs brutos (fora do git, por exemplo em `~/compass_bench/agulha_NN/`);
- as saídas JSON do `bench_fit.py`;
- campo das bobinas por ampère (calibração), data e temperatura.

Uma tabela-resumo por agulha (I, m, b, τ_c, τ_s, com incertezas) pode entrar no repositório em `docs/bench/`.
