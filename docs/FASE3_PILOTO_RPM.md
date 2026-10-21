# Fase 3: piloto de memória de ponto de retorno (RPM)

Data: 19/10/2026. Branch `fase3/qs-subloops` (motor 2.5.0). Status: script pronto (`fase3_pilot.py`: `run` e `analyze`), testado em 6×6; piloto rodado na mumax em 19–20/10 (60 corridas, 22,3 h, sem erros; `~/compass_runs/fase3_pilot`). **Resultado: nenhuma violação de RPM** (seção 6).

## 1. Pergunta

Quando o campo volta a um ponto de virada anterior depois de um desvio, a rede volta ao mesmo estado (memória de ponto de retorno, RPM)? E a inércia (Q₀ alto) quebra essa memória?

O teorema de Sethna garante RPM para dinâmica superamortecida com interações "ferromagnéticas" (propriedade de *no-passing*). A interação dipolar é frustrada, então o fluxo gradiente pode já violar a RPM. Por isso o efeito inercial é medido como a **diferença** entre conjuntos de realizações com e sem inércia, no mesmo protocolo.

Pelo Relatório V61, a inércia só muda o laço quase-estático a partir de N = 24. O piloto usa 24×24.

## 2. Programa de campo

`--qs_program "4, A, A+Δp, A, A+Δg, A, -4"` (unidades de B_ref, passo de campo 0,005):

1. saturar em +4;
2. descer até A, no meio da reversão do ramo descendente;
3. desvio pequeno até A+Δp e volta a A: teste **RPM-1** (estado em A igual ao da primeira passagem);
4. desvio grande até A+Δg e volta a A: teste **RPM-2**;
5. continuar até −4: teste **RPM-3**, comparando o ramo com uma corrida de **referência** `"4,-4"` com a mesma semente, nos mesmos campos (`rpm.py --ref ... --legs 7:2`).

Parâmetros escolhidos a partir dos laços de 24×24 da grade QS2 (fluxo gradiente, semente 0):

| geometria | M de +0,5 a −0,5 no ramo descendente | A | Δp (≈ ¼ do trecho) | Δg (≈ ¾ do trecho) |
|---|---|---|---|---|
| quadrada | +0,235 → −0,34 | −0,10 | 0,145 | 0,43 |
| triangular | +0,015 → −0,36 | −0,20 | 0,095 | 0,28 |
| honeycomb | +0,055 → −0,24 | −0,20 | 0,075 | 0,22 |

Todos os valores são múltiplos de 0,005, de modo que os pontos de virada e o ramo de referência caem na mesma grade.

## 3. Corridas

- 3 geometrias × 2 dinâmicas (fluxo gradiente, Q₀ = 10) × 5 sementes × 2 corridas (programa e referência) = **60 corridas**.
- `--qs_state_log` em todas; amortecimento ajustado para Q₀ exato em cada geometria (como na Fase 2); `--damping_noise 0`.
- A saturação apaga o estado inicial (QS2), então a semente sorteia, na prática, a sequência dos chutes de 10⁻⁴ rad da verificação de sela.

**Custo estimado:** cerca de 60 % de um laço completo da QS2 por corrida (≈ 40 min no fluxo gradiente, ≈ 75 min em Q₀ = 10, no Mac); ≈ 57 h de CPU, **4–6 h de relógio** com o Mac (10 processos) e a mumax (7). Registros de estado: ≈ 0,7 GB, fora do git.

## 4. Leitura

Para cada condição, `rpm.py` dá a fração de agulhas em estado diferente (e o ΔM) nos três testes.

| resultado | interpretação | próximo passo |
|---|---|---|
| fluxo gradiente ≈ 0 e Q₀ = 10 > 0 na maioria das sementes, em alguma geometria | a inércia quebra a memória | campanha completa (N = 32, Q₀ = 3, mais sementes) |
| fluxo gradiente também viola | a frustração já quebra a memória; o efeito inercial é a diferença | quantificar com mais sementes |
| ninguém viola | memória robusta nessas condições | desvios maiores ou N = 32 |

**Cuidado:** os chutes de 10⁻⁴ rad também decidem caminhos perto de bifurcações (`docs/QS_SELECAO_DE_CAMINHO.md`). Uma violação no fluxo gradiente não é efeito de inércia; o controle existe para isso.

## 5. Como rodar

```bash
python3 fase3_pilot.py run --out_dir ~/compass_runs/fase3_pilot --workers 7
python3 fase3_pilot.py analyze --run_dirs ~/compass_runs/fase3_pilot
```

A análise grava `rpm_pilot.csv` (uma linha por geometria, dinâmica e semente) e imprime quantas sementes violam cada teste. Teste em 6×6 (fluxo gradiente, semente 0): nenhuma violação nas três geometrias, como esperado em rede pequena.

## 6. Resultado (21/10/2026)

**Nenhuma das 30 combinações (3 geometrias × 2 dinâmicas × 5 sementes) viola a RPM em nenhum dos três testes.** Ao voltar a A, o estado é o mesmo da primeira passagem (diferença rms de 3×10⁻⁶ a 2×10⁻⁵ rad, no nível da tolerância de repouso, e nenhuma agulha em estado diferente), e o ramo externo depois dos desvios coincide com o da corrida de referência em todos os 760–780 campos comparados.

Os desvios não foram triviais (semente 0, fração de agulhas que viraram durante o desvio e voltaram):

| geometria | dinâmica | M em A | desvio pequeno | desvio grande |
|---|---|---|---|---|
| quadrada | fluxo gradiente | −0,085 | 37 % | 67 % |
| quadrada | Q₀ = 10 | −0,015 | 6 % | 44 % |
| triangular | fluxo gradiente | −0,46 | 0,3 % | 50 % |
| triangular | Q₀ = 10 | −0,46 | 0,3 % | 50 % |
| honeycomb | fluxo gradiente | −0,18 | 2 % | 19 % |
| honeycomb | Q₀ = 10 | −0,18 | 2 % | 19 % |

Na quadrada, as duas dinâmicas chegam a A em estados diferentes (M = −0,085 contra −0,015): a inércia muda o caminho, como no V61. Mas cada dinâmica tem memória perfeita do seu próprio caminho.

**Leitura:** é o terceiro cenário da seção 4. Nessas condições a memória de ponto de retorno é robusta, apesar da frustração dipolar (que viola a hipótese de *no-passing* do teorema de Sethna) e da inércia. O desvio pequeno da triangular quase não virou agulhas (2 de 576) e é um teste fraco.

**Limites do piloto:** um único ponto A por geometria, desvios dentro do primeiro ramo descendente, só N = 24, só Q₀ = 10 contra o fluxo gradiente, sem desordem congelada.

**Opções para continuar:**
1. Desvios que atravessem a região em que a inércia seleciona caminhos no V61 (por exemplo, desvios a partir de vários pontos A ao longo de todo o ramo, ou subloops aninhados).
2. N = 32, onde o efeito inercial no laço é maior (até −21 %).
3. Desordem fraca (`--damping_noise`) e Q₀ = 3.
4. Tratar a robustez da RPM como resultado: memória exata num sistema frustrado e inercial, com contraste à literatura de no-passing.
