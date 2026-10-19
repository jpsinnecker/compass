# Fase 3: piloto de memória de ponto de retorno (RPM)

Data: 19/10/2026. Branch `fase3/qs-subloops` (motor 2.5.0). Status: **desenho**; o script da campanha e a análise ainda não foram escritos, e nada foi rodado.

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

## 5. A fazer antes de lançar

1. Script da campanha (reaproveitando a divisão Mac/mumax e o ajuste de Q₀ do `fase2_campaign.py`).
2. Análise que gera a tabela da seção 4 a partir do `rpm.py`.
