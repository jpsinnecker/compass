# Fase 2: histerese dinâmica — resultados preliminares

Data: 17/10/2026. Branch `fase2/dynamic-hysteresis`. Scripts: `fase2_campaign.py`, `fase2_analysis.py`.
Dados (fora do git): `~/compass_runs/fase2` (Mac, metade 1), `~/compass_runs/fase2_mumax` (mumax, metade 2), `~/compass_runs/fase2_analysis` (tabelas e figuras).
A extensão para redes grandes (48×48 e 64×64, GPU da mumax) ainda está rodando e **não** entra aqui.

## 1. O que foi rodado

- **600 ciclos completos de histerese.** Geometrias: quadrada, triangular e honeycomb. Tamanhos: N = 6, 10, 16, 24, 32. Q₀: 0,1 / 0,3 / 1 / 3 / 10. Taxa adimensional r: 10⁻¹ a 10⁻³, e também 3×10⁻⁴ para N ≤ 16. Mais duas sementes extras em N = 10.
- **Parâmetros fixos:** B_max = 4 B_ref, espessura 0,4 mm, contorno livre, soma dipolar completa. O amortecimento foi escolhido para que o Q₀ intrínseco de cada geometria e tamanho fosse **exatamente** o alvo.
- **Máquinas:** Mac M2 Max (359 corridas, 75 min) e mumax i9-9900K (237 corridas, 1 h 50).
- **Controle de qualidade, em todas as 600 corridas:**
  - erro de Q₀ ≤ 4×10⁻¹⁶;
  - resíduo do balanço de energia ≤ 1,2×10⁻³;
  - zero passos sinalizados pelo monitor de estabilidade;
  - M ≥ 0,995 em B_max;
  - área pelo livro-caixa = área por ∮M dB, com diferença ≤ 5×10⁻³.

A área do laço `a` é dada por agulha, em unidades de m·B_ref.

## 2. Resultados

### 2.1 No regime superamortecido, Q₀ não importa quando a taxa é medida corretamente
Com r medido em unidades do tempo de relaxação (viscoso, para Q₀ < 1), os laços para Q₀ = 0,1 e 0,3 são praticamente idênticos. A diferença de área fica entre 0,05 % e 2,4 % em todas as geometrias, tamanhos e taxas. Q₀ = 1 já está no cruzamento entre os regimes.

Isso era esperado: no limite superamortecido, a equação se reduz a b θ̇ = τ, e o amortecimento apenas reescala o tempo.

**Consequência:** a dependência área × Q da figura antiga (`figs/area_vs_Q.png`, feita com `t_sim` fixo) era, para Q < 1, quase toda efeito de taxa (atraso viscoso), e não física coletiva. Esse resultado sustenta a auditoria e precisa constar no artigo como resultado de controle.

### 2.2 Efeito inercial genuíno para Q₀ ≳ 3: laços mais estreitos na mesma taxa
Na mesma taxa r, aumentar Q₀ de 0,1 para 10 **reduz** a área do laço:

| geometria, N | r = 10⁻¹ | r = 10⁻² | r = 10⁻³ |
|---|---|---|---|
| quadrada, 16 | −44 % | −37 % | −30 % |
| quadrada, 32 | −38 % | −54 % | −5 % |
| triangular, 16 | −23 % | −28 % | −22 % |
| triangular, 32 | −24 % | −49 % | −34 % |
| honeycomb, 16 | −35 % | −33 % | −7 % |
| honeycomb, 32 | −38 % | −40 % | −11 % |

- Os laços subamortecidos mostram oscilações (ringing) e saltos discretos (avalanches); ver `fig4_loops`.
- A redução diminui à medida que r → 0 em honeycomb e na quadrada 32, mas **persiste na triangular** (−22 % a −34 % em r = 10⁻³). É um candidato a efeito inercial que sobrevive ao limite quase-estático e depende da frustração.
- **Isso ainda não está demonstrado** (ver 2.3).

### 2.3 O limite quase-estático não foi alcançado nas redes maiores
- Para N ≥ 24, a área muda de 40 % a 98 % entre os dois menores r (3×10⁻³ e 10⁻³). Nesses casos, o ajuste a(r) = a₀ + c r^β não é confiável: a₀ frequentemente bate no limite zero.
- Só honeycomb N = 16 está convergida dentro de ~5 % em r = 3×10⁻⁴.
- O expoente β fica em 0,6–0,9, mas não deve ser citado enquanto a₀ for incerto.

**Por que isso acontece:** reduzir r em 10× custa 10× mais passos. Para Q₀ < 1, o custo é ainda 1/Q₀ vezes maior. Chegar a r ~ 10⁻⁵ em 32×32 por varredura contínua é inviável.

### 2.4 Sensibilidade ao estado inicial (honeycomb subamortecida)
- Com amortecimento uniforme, sementes diferentes deveriam dar o mesmo laço, porque a saturação apaga a memória. Isso vale em 40 de 69 condições.
- Na **honeycomb com Q₀ ≥ 1**, a área varia até **27 %** entre sementes. Em B_max = 4 B_ref resta uma textura residual (M ≈ 0,9985), e a dinâmica inercial a amplifica em cascatas diferentes.
- Para essa rede, é obrigatório fazer média sobre sementes ou usar um B_max maior. Fisicamente, é uma forma de memória que só existe com inércia: um bom gancho para a Fase 3.

### 2.5 Tamanho finito
A área não varia monotonicamente com N. Na triangular, por exemplo, em r = 10⁻³: 0,51 (N=6), 0,31 (10), 0,90 (16), 0,66 (24), 0,73 (32). Efeitos de borda e de forma da amostra são fortes em contorno livre, e a forma do recorte triangular muda com N. Ainda não há base para extrapolar para o limite de rede infinita.

## 3. Recomendação principal: protocolo quase-estático por eventos

Em vez de tentar r → 0 por varredura contínua, implementar o limite quase-estático **diretamente**:

1. aumentar o campo em pequenos degraus ΔB;
2. após cada degrau, integrar a dinâmica real (com o Q₀ verdadeiro) a campo fixo até o repouso, ou seja, até ω_max e o torque ficarem abaixo de um limiar;
3. registrar o estado e seguir para o próximo degrau.

Vantagens:
- **Define a₀(Q₀) exatamente:** é o limite r → 0 em que a inércia ainda escolhe o destino de cada avalanche (overshoot entre bacias), sem nenhum atraso de taxa.
- **É muito mais barato:** o tempo de simulação só é gasto onde há avalanche. O custo não depende de r e quase não depende de Q₀.
- **Separa as avalanches de forma natural:** cada degrau de campo produz um evento discreto (tamanho, energia dissipada, região). É a estatística limpa que a análise de avalanches precisa.
- **É o protocolo natural para a Fase 3** (memória de ponto de retorno com subloops).
- **Em Q₀ → 0** coincide com o protocolo atérmico quase-estático (AQS) clássico. Isso dá uma comparação com a literatura de histerons e RFIM, e um teste de validação do código.

Implementar em `compass.py` exige uma nova opção de protocolo (por exemplo `--field_mode hysteresis_qs --qs_dB ... --qs_rest_tol ...`) e testes. Depois, refazer a grade Q₀ × geometria × N com sementes nesse protocolo.

## 4. Próximos passos sugeridos

1. Implementar o protocolo quase-estático por eventos, com testes (incluindo o limite Q₀ → 0).
2. Rodar com ele a grade Q₀ × geometria × N (com sementes) e medir a₀(Q₀). Testar se a redução inercial da triangular (2.2) sobrevive.
3. Estudar o tamanho finito com o novo protocolo, incluindo as redes grandes na GPU, onde ele também fica muito mais barato.
4. Manter a campanha de taxa atual como resultado de controle (2.1) e como medida da histerese dinâmica (A − A₀ em função de r).

## 5. Figuras (em `~/compass_runs/fase2_analysis/figs/`)
- `fig1_area_vs_rate`: a(r) por geometria e Q₀, N = 16.
- `fig2_a0_vs_Q0`: a₀ ajustado × Q₀ por tamanho (não confiável para N ≥ 24; ver 2.3).
- `fig3_excess_area`: a − a₀ × r em log-log.
- `fig4_loops`: laços M(B) para Q₀ = 0,1 e 10 em três taxas.
