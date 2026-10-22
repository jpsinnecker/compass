# Artigo B: roteiro

Data: 21/10/2026. Base: Relatório V61 (resultados I–III), piloto de RPM da Fase 3 (`docs/FASE3_PILOTO_RPM.md`) e o teste de estresse em execução. O texto do artigo será em inglês; este roteiro é para o autor.

**Rascunho:** `docs/paper_B/paper_B.tex` (REVTeX, PRE, duas colunas; PDF ao lado). Seções I–V redigidas em 22/10 a partir do V61; VI (memória), VII (discussão) e o apêndice sobre selas estão marcados em vermelho, à espera do teste de estresse.

## 1. Afirmação central

> Em redes dipolares de rotores inerciais, a inércia controla a histerese **dinâmica** e, em redes grandes, **seleciona o ramo** do laço quase-estático; a memória de ponto de retorno é preservada [a confirmar pelo teste de estresse].

Títulos possíveis:
- *Inertia selects reversal paths but preserves memory in dipolar rotor arrays*
- *Dynamic and quasi-static hysteresis of inertial dipolar rotor lattices*
- *Separating rate, inertia and memory in the hysteresis of compass-needle arrays*

O primeiro só vale se o teste de estresse confirmar a RPM nas redes grandes. Se aparecerem violações, o título muda (seção 6).

**Revista:** Physical Review E (alternativas: PRB, Physical Review Research). Formato: artigo regular, ~10 páginas, 6–7 figuras.

## 2. Estrutura

| seção | conteúdo | material existente | falta |
|---|---|---|---|
| I. Introdução | arranjos dipolares macroscópicos; rotores com inércia (Mellado 2012); histerese dinâmica (Rao 1990); memória e no-passing (Sethna 1993, Keim 2019); lacuna: separar A₀, A − A₀ e o papel da inércia em cada um | V61, introdução | revisar as três referências marcadas no V61; incluir a literatura de RPM sem no-passing |
| II. Modelo e escalas | rotores planares, soma dipolar completa, Q₀ intrínseco, taxa adimensional r | V61, seção 2 | figura das redes e da agulha (já existe) |
| III. Métodos | Verlet com amortecimento implícito, passo pela cota do campo local, livro-caixa de energia; protocolo quase-estático por eventos **com verificação de sela**; programas de campo | V61, seções 3–4 | parágrafo curto sobre os programas de campo (motor 2.5.0) |
| IV. Histerese dinâmica | colapso superamortecido; estreitamento inercial de 20–50 % em taxa finita; o limite r → 0 inacessível por varredura | V61, resultados II; figs. 1 e 4 | — |
| V. Limite quase-estático | Q₀ = 1 igual ao fluxo gradiente; sem efeito inercial até 16×16; estreitamento de 9–21 % em 24 e 32 por seleção de caminho; dependência da realização na quadrada | V61, resultados III; figs. QS1–QS2 | mapa espacial de um evento de seleção (já existe em `docs/figs/qs_paths/`) |
| VI. Memória de ponto de retorno | piloto (exato em 24×24 para desvios que viram até 67 % das agulhas) e teste de estresse (vários pontos de virada, subloop aninhado, N = 6–32) | piloto | **resultado do teste de estresse**; figura nova (fração de retornos que violam, por N e dinâmica) |
| VII. Discussão | inércia seleciona caminhos mas não apaga memória; hipótese do acoplamento cooperativo cabeça-cauda; relação com no-passing; previsões para a bancada | — | testar a hipótese com os estados registrados |
| Apêndice A | artefato de pontos de sela: como aparece e como se corrige | `docs/QS_SELECAO_DE_CAMINHO.md` | condensar |
| Apêndice B | validação numérica (testes físicos) | V61, validação | — |

## 3. Figuras

1. Redes e agulha (existente: `figs/`, fig. de geometria do V61).
2. Área do laço × taxa r, N = 16, cinco Q₀ (existente: `fig1_area_vs_rate`).
3. Laços M(B) superamortecido × subamortecido em três taxas (existente: `fig4_loops`).
4. a₀ × N e efeito inercial quase-estático (existentes: `figQS1`, `figQS2`; talvez juntar em uma figura de dois painéis).
5. Mapa espacial de um evento de seleção de caminho (existente: `docs/figs/qs_paths/`).
6. **Nova:** memória de ponto de retorno, com um exemplo de programa de campo e o estado nos retornos, e a fração de retornos violados por N e dinâmica.
7. Opcional: distribuição de avalanches por passo (`figQS3`) — pode ir para material suplementar, sem expoentes.

## 4. O que não afirmar (regras do projeto)

- Nenhum expoente de avalanche nem "classe de universalidade".
- Não é modelo de spin ice artificial nanomagnético (LLG ≠ Newton); ASI só como motivação.
- Efeitos de Q₀ sempre no limite quase-estático ou em função explícita de r.

## 5. Pendências antes da submissão

1. Resultado do teste de estresse (seção VI e figura 6).
2. Verificar as referências `menniti2025`, `taniguchi2025`, `concha2017`.
3. Decidir a espessura da agulha (0,26 × 0,4 mm); os resultados adimensionais não mudam, mas os valores em SI do texto sim.
4. Robustez mínima, se um revisor pedir: modelo de halteres num caso (Fase 4).
5. Bancada (Fase 5): se houver medidas a tempo, uma comparação de Q₀ real e um teste de RPM numa rede pequena valorizam muito o artigo; se não, entram como previsão.

## 6. Se o teste de estresse mostrar violações

- **Violações só no fluxo gradiente e em redes pequenas** (como no teste em 6×6): a frustração quebra a memória em redes pequenas e simétricas, e o efeito desaparece com o tamanho. A seção VI vira "memória emergente com o tamanho".
- **Violações que dependem da inércia:** o resultado sobe de nível e justifica um artigo separado sobre memória (o artigo A do roteiro), com o artigo B ficando com as seções IV–V.
