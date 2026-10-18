# Auditoria do relatório científico (`docs/project_report.tex`)

> **Status (17/10/2026):** a reestruturação proposta foi executada no Report V60, com resultados da Fase 2 e do protocolo quase-estático. As três referências marcadas em vermelho no PDF (`menniti2025`, `taniguchi2025`, `concha2017`) ainda precisam ser verificadas.

Data: 17/10/2026. Versão auditada: Report V53, com as notas de outubro (24 páginas, compila sem erros a partir de `docs/`).
O manual de uso passou a ter fonte única em `USER_GUIDE.md` (PDF em `docs/USER_GUIDE.pdf`). Esta auditoria trata só do relatório.

---

## 1. Diagnóstico geral

O relatório mistura quatro gêneros num documento só:

1. **Caderno de laboratório:** cronologia, motivação de cada mudança.
2. **Manual de software:** comandos, tabela de flags, formato do CSV e dos metadados, arquivos de saída. Isso **duplica** o `USER_GUIDE.md` e já envelheceu: comandos e tabela estão desatualizados.
3. **Proposta de pesquisa:** objetivos, cinco estudos propostos, plano de validação.
4. **Resultados:** **ausentes**. Não há nenhuma figura de dados, nem os achados da Fase 1 (as medidas da auditoria) e da Fase 2 (600 laços, colapso superamortecido, efeito inercial, limite quase-estático).

Para servir de base a um artigo, o documento precisa se tornar **relatório de métodos e resultados**. O material de software sai e vira referência ao manual. As hipóteses que os dados já testaram devem ser substituídas pelos resultados.

---

## 2. Erros e informações desatualizadas (com localização)

| # | Onde (linha aprox.) | Problema | Correção |
|---|---|---|---|
| 1 | §3.6, l. 360–362 | Afirma que 1/r³ em 2D é "apenas condicionalmente convergente". **Errado**: a soma dipolar no plano é absolutamente convergente (∫ r dr / r³), como diz o próprio parágrafo de atualização logo acima. | Reescrever o parágrafo; manter só a ressalva sobre PBC e efeitos de forma. |
| 2 | §3.4, l. 280–288 | ω₀ com B_eff aparece como "frequência natural de referência". Desde outubro, ele não define mais nem dt nem Q; é só unidade de tempo para *dwell* e taxas. | Rebaixar ω₀ a unidade e explicar onde ainda é usado. |
| 3 | §3.4, l. 308 | Diz que "Q ~ 1 é o regime mais sensível". Os dados da Fase 2 refutam: para Q₀ ≤ 1 tudo colapsa quando a taxa é medida em unidades viscosas. | Trocar a hipótese pelo resultado. |
| 4 | §4, l. 382–390 | Título "Research branch compass.py V2.0.2" e a frase "compass.py foi criado... não pretende substituir o compass.py". Isso é resíduo da renomeação `compassV2.py → compass.py` e ficou sem sentido. | Remover a seção. O histórico vai para um apêndice curto. |
| 5 | Tab. 1, l. 414–415 | Descreve `flip_field`/`flip_angle` como contagens simples. Os contadores atuais são os "endurecidos" da V2.1 (Schmitt + dwell). | Atualizar, ou remover e remeter ao manual. |
| 6 | §4.3, l. 424–435 | A lista de protocolos não tem `hysteresis_qs`, e `--sweep_rate_dimless` não é citado no corpo. | Os dois são centrais no método atual. |
| 7 | §4.4, l. 452 | A lista de metadados não tem B_lat, Q_ref, Q_field, T_loc, `sweep_rate`, `energy_ledger`, `quasi_static`. | Remeter ao manual. |
| 8 | §7, l. 554–665 | **Todos os comandos** usam `--cutoff_shells 3.5` (desaconselhado), `--dt_factor 0.035/0.04` (semântica antiga) e histerese por `--t_sim` (taxa não controlada). Contrariam diretamente as regras do projeto. | Remover. Os exemplos corretos estão no manual. |
| 9 | §7.6, l. 667–669 | "Status do `damping_sweepV02.py`": obsoleto. O V02 não roda mais; existem o V03 e o `fase2_campaign.py`. | Remover. |
| 10 | §8, l. 678–688 | Plano de validação escrito como futuro. Quase tudo já virou teste automático (`tests/test_physics.py`); o item de cutoff ficou obsoleto; os valores de `dt_factor` estão na semântica antiga. Faltam a convergência em r e a validação do modo quase-estático. | Transformar em seção de **validação realizada**, com números (ver §3). |
| 11 | §9, l. 695–706 | Os nomes dos arquivos de análise (`summary_loops.csv` etc.) não correspondem ao pipeline real (`manifest.jsonl`, `runs.csv`, `fits.csv`). | Atualizar ou remover. |
| 12 | §11, l. 748–756 | Limitações incompletas: faltam dipolo pontual (L/a ≈ 0,77), atrito de Coulomb, efeitos de borda e tamanho finito, sensibilidade à semente na honeycomb subamortecida. "Wrappers must be updated to V2" está obsoleto. | Reescrever. |
| 13 | §12, l. 763–772 | Próximos passos obsoletos. Pede a figura "loop area vs Q", que hoje se sabe ser dominada pela taxa. | Substituir pelo roteiro atual. |
| 14 | Apêndice A | Tabela de flags duplicada do manual e incompleta (faltam `--sweep_rate_dimless` e `--qs_*`). A descrição do cutoff não diz que ≤ 0 = soma completa. Afirma que a espessura de 0,4 mm difere "de propósito" da de 0,26 mm, o que **contradiz** a auditoria: a diferença não está conciliada. | Remover o apêndice e remeter ao manual. A espessura vira pendência explícita. |
| 15 | Apêndice C | O histórico de versões para em julho (V53). Não tem as mudanças de outubro, que alteram os números. | Acrescentar, ou trocar por referência ao git. |
| 16 | Fig. 1 | Desenha a agulha com t = 0,26 mm, enquanto as simulações usam 0,4 mm. | Alinhar com a decisão sobre a espessura. |

---

## 3. Lacunas científicas

**3.1 Afirmação de novidade sem a referência mais próxima.**
A l. 110 diz que "a dinâmica dipolar inercial é um eixo de controle ausente". O relatório **não cita** Mellado, Concha & Mahadevan, *Macroscopic magnetic frustration*, PRL **109**, 257203 (2012), que já estudou rotores magnéticos macroscópicos com inércia e dissipação explícitas. Um revisor apontaria isso de imediato. Faltam também:
- histerese dinâmica e leis de escala da área com a taxa (Rao, Krishnamurthy & Pandit, PRB 1990; trabalhos recentes, p. ex. PRE **108**, 024101, 2023);
- cadeias dipolares de rotores (Concha et al.);
- memória de ponto de retorno (Sethna et al., PRL 1993) e formação de memória em matéria (Keim et al., RMP 2019), que fundamentam a Fase 3.

**3.2 Referências a verificar.**
- `taniguchi2025` (Sci. Rep. 15, 9073) e `menniti2025` (PRB 112, 134410): conferir volume, página e DOI.
- `dahmen2011` trata de avalanches em materiais granulares sob cisalhamento. Como suporte para "expoente de avalanche × Q", é uma escolha estranha; avaliar se cabe.

**3.3 Os resultados obtidos não aparecem.** Devem entrar, com figuras:
- **Auditoria numérica (Fase 1):** viés do cutoff (15–17 %), erro de energia com o dt antigo (1,4–1,8 %), Q dependente do protocolo (√8), área ∝ taxa^~0,6.
- **Fase 2:** colapso para Q₀ ≤ 1 com a taxa em unidades viscosas; estreitamento inercial para Q₀ ≳ 3 em taxa finita; limite quase-estático não atingido por varredura para N ≥ 24; sensibilidade à semente.
- **Protocolo quase-estático:** validação (fluxo gradiente = inercial com Q₀ baixo; concordância com as varreduras lentas) e o primeiro indício de que, no limite quase-estático, a inércia muda pouco o laço (−0,7 % a −5,5 %). A grade em andamento deve confirmar.

**3.4 A pergunta de pesquisa precisa ser reformulada à luz dos dados.** A pergunta original era "o Q controla a histerese?". Os dados sugerem uma versão mais precisa e mais forte:

> No limite quase-estático, a inércia muda pouco o laço; o grande efeito dela está na **histerese dinâmica** (taxa finita) e, possivelmente, na **seleção de caminhos e na memória** (sensibilidade ao estado inicial na honeycomb).

Isso conecta diretamente com as Fases 3 e 4 do roteiro.

---

## 4. Estrutura proposta para o relatório revisado

1. **Resumo:** pergunta, método, principais resultados com números.
2. **Introdução e literatura:** agulhas de bússola, rotores macroscópicos (Mellado 2012), histerese dinâmica, ASI só como motivação. Lacuna formulada com precisão.
3. **Modelo:** equação de movimento, interação dipolar, parâmetros da agulha (pendência da espessura), escalas B_ref, B_lat, Q₀, τ_relax e taxa r.
4. **Métodos numéricos:** soma completa, Verlet com amortecimento implícito, dt pela cota do campo local, livro-caixa de energia, protocolos (varredura com r controlado; quase-estático por eventos; fluxo gradiente).
5. **Validação:** testes automáticos e números (energia de 2ª ordem, balanço, cota de dt, Q₀ invariante, convergência em r, validação do quase-estático).
6. **Resultados I, auditoria:** o que estava errado e quanto isso afetava as conclusões antigas.
7. **Resultados II, histerese dinâmica (Fase 2):** a(r, Q₀) por geometria e tamanho; colapso superamortecido; estreitamento inercial.
8. **Resultados III, limite quase-estático:** a₀(Q₀) por geometria e tamanho (grade em andamento); avalanches por degrau.
9. **Discussão:** o que é efeito de taxa e o que é inercial; comparação com histerese dinâmica e com rotores macroscópicos.
10. **Limitações:** dipolo pontual, atrito, T = 0, bordas e tamanho finito, sementes.
11. **Próximos passos:** memória e RPM, modelo de halteres, atrito de Coulomb, experimento de bancada.
12. **Apêndices:** histórico condensado com ponteiro para o git; tabela de reprodutibilidade (commit, comandos das campanhas, onde estão os dados).

O que **sai**: comandos, tabela de flags, formato do CSV e dos metadados, arquivos de saída. Fica uma frase remetendo ao `USER_GUIDE.md`.

## 5. Como executar a revisão

- Fazer a reescrita **depois que a grade quase-estática terminar** (em andamento), para que a seção de resultados III já entre com dados.
- Figuras de resultados geradas por script versionado (`fase2_analysis.py`) e salvas como PDF em `docs/figs/` (PDF não é bloqueado pelo `.gitignore`).
- Corrigir agora os erros factuais simples (itens 1, 4 e 14 da tabela) não traz ganho: a reestruturação vai substituir essas seções de qualquer forma.
