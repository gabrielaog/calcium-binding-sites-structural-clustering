# Pipeline — Classificação Estrutural de Sítios de Ligação ao Cálcio

Documento de referência para retomada do projeto. Descreve as etapas do pipeline, os scripts utilizados em cada uma, e as justificativas para as escolhas metodológicas.

---

## Visão Geral

O objetivo é classificar sítios de ligação ao cálcio (Minimal Functional Sites, MFS) com base em similaridade estrutural, seguindo a metodologia de Andreini et al. (2011) para zinco.

**Dado de entrada:** estruturas PDB de MFS extraídas do MetalPDB, organizadas por superfamília (CATH e SCOPe).
**Dado de saída:** 11 clusters globais de MFS representativos, com threshold d < 0.3 (TM-score ≥ 0.7).

---

## Etapa 1 — Alinhamento estrutural intra-superfamília

**Ferramenta:** MM-align
**O que faz:** Compara pares de MFS dentro de cada superfamília e gera arquivos `.txt` com os dois TM-scores do alinhamento.
**Saída:** Um arquivo `.txt` por par de sítios, nomeado `site1_site2.txt`.

**Métrica de similaridade usada:**
```
TM-score_médio = (TM-score_1 + TM-score_2) / 2
distância = 1 − TM-score_médio
```
A média dos dois TM-scores é usada porque MM-align normaliza cada score pelo comprimento de uma das estruturas; a média garante simetria e robustez.

---

## Etapa 2 — Clustering intra-superfamília (Stage 1)

**Script:** `similarity_matrix_and_clustering.py`
**Algoritmo:** Single-linkage hierárquico (scipy), threshold = 0.5
**O que faz:**
1. Para cada superfamília, lê os arquivos `.txt` do MM-align e monta uma matriz de similaridade.
2. Converte para matriz de distância (d = 1 − TM-score).
3. Aplica single-linkage com threshold 0.5.
4. Seleciona um representante por cluster.

**Saídas:**
- `similarity_matrices_CATH.xlsx` / `similarity_matrices_SCOPe.xlsx` — matrizes de similaridade por superfamília.
- `CATH_all_clusters05.txt` / `SCOPe_all_clusters05.txt` — atribuição de clusters por sítio.

**Erro conhecido no código (não afeta resultado):**
O comentário na linha ~78 diz `# Missing alignment assumed as 0 similarity` mas está colado ao bloco que trata `i==j` (diagonal = 1.0), não ao bloco de pares ausentes. É apenas um erro de posicionamento de comentário.

**Por que single-linkage aqui?**
Dentro de uma superfamília os sítios já compartilham fold; o critério mais permissivo (menor distância entre quaisquer dois membros) é adequado para agrupar variantes estruturais próximas antes de selecionar representantes.

**Por que threshold 0.5?**
Threshold intermediário que captura similaridade estrutural significativa sem colapsar toda a superfamília num único cluster.

---

## Etapa 3 — Seleção de representantes e resolução PDB

**Scripts:** `APIresolutionsCATH.py`, `APIresolutionsSCOPe.py`
**O que fazem:** Consultam as APIs do CATH e SCOPe para obter a resolução cristalográfica de cada estrutura PDB.
**Critério de seleção:** Dentro de cada cluster intra-superfamília, o sítio com a **menor resolução** (melhor qualidade cristalográfica) é escolhido como representante.

**Erros conhecidos no código (não afetam resultado):**
- `APIresolutionsCATH.py`: contém `print(f"match_superfamily {match_superfamily}")` — print de debug deixado em produção.
- `APIresolutionsSCOPe.py`: contém `print(f'nfo pdb {info_pdb}')` — print de debug com typo ("nfo" em vez de "info").
- Os arquivos de saída (`best_resolutionsCATH.txt`, `best_resolutionsSCOPe.txt`) têm nomes diferentes dos que `efhand_check.py` lê (`CATH_Representatives_SelectedByBestResolution_TMscore_0.5.txt`, `SCOPe_Representatives_SelectedByBestResolution_TMscore_0.5.txt`). Houve uma etapa manual de conversão/renomeação não documentada.

---

## Etapa 4 — Alinhamento estrutural global

**Ferramenta:** MM-align (rodada novamente)
**O que faz:** Compara todos os pares de representantes selecionados na Etapa 3, entre todas as superfamílias.
**Saída:** `DistanceMatrix.csv` — matriz de distâncias N×N entre todos os representantes.

**Valores especiais na matriz:**
- `d = 0.0` na diagonal (auto-distância).
- `d = 1.0` para pares sem alinhamento válido (MM-align não convergiu); esses valores são tratados como ausentes, não como distância máxima real.

---

## Etapa 5 — Clustering global (Stage 2)

**Script:** `global_clustering_tm07.py`
**Algoritmo:** Greedy nearest-neighbour (NÃO é single-linkage)
**Threshold:** d < 0.3 (equivalente a TM-score ≥ 0.7)

**Como funciona (passo a passo):**
1. Para cada sítio (linha da matriz), encontra o vizinho mais próximo (coluna com menor distância, excluindo d=0).
2. Se ambos já estão em clusters → ignora (sem fusão de clusters).
3. Se a distância for d=1.0 ou d≥0.3 → ignora.
4. Se apenas um já tem cluster → adiciona o outro ao mesmo cluster.
5. Se nenhum tem cluster → cria um novo cluster com os dois.
6. Após a varredura, descarta clusters com menos de 10 membros.

**Saída:** `resultadoMaiorQue07.txt` — lista de sítios com seu cluster, vizinho mais próximo e distância.

**Erro conhecido no código (não afeta resultado):**
O docstring diz "similar to single-linkage" — está errado. O algoritmo é greedy NN. A lógica do código em si está correta.

**Por que greedy NN e não single-linkage?**
Single-linkage sofre do *chaining effect*: um único par com distância baixa pode encadear clusters muito heterogêneos. O greedy NN com a regra "se ambos já têm cluster → ignora" impede fusões indiretas, garantindo que cada sítio entre no cluster apenas pelo seu vizinho mais próximo imediato.

**Por que descartar clusters com < 10 membros?**
Para focar em arquiteturas estruturalmente recorrentes. Clusters pequenos representam configurações raras ou isoladas sem relevância estatística.

---

## Etapa 6 — Análise de threshold (validação da escolha d=0.3)

**Script:** `threshold_analysis.py`
**O que faz:** Roda o mesmo algoritmo greedy NN da Etapa 5 para 9 thresholds (0.20 a 0.60) e computa 4 métricas de qualidade.

**Métricas e justificativas:**

| Métrica | O que mede | Por que usar |
|---------|-----------|--------------|
| Nº de clusters válidos (≥10 membros) | Resolução da partição | Thresholds muito baixos colapsam tudo; muito altos fragmentam demais |
| Tamanho médio dos clusters | Compacidade | Crescimento rápido indica chaining effect |
| % distâncias intra-cluster < 0.5 | Coerência estrutural interna | Clusters com muitas distâncias ≥ 0.5 são heterogêneos |
| Silhouette score | Separação entre clusters | Score positivo = clusters bem separados; negativo = sobreposição |

**Cálculo do `pct_below_0.5`:**
É uma **média ponderada** pelo número de pares de cada cluster:
```
pct = (Σ_k  #{pares em C_k com d < 0.5}) / (Σ_k #{pares em C_k})  × 100
```
Isso dá peso igual a cada par (não a cada cluster). Uma média simples por cluster daria valores inflados pelos clusters pequenos — esse foi o erro da Table 4 original (que usava média simples e resultava em 76.08% ao invés de 60.46%).

**Saídas:**
- `threshold_analysis.csv` — tabela com as 4 métricas para cada threshold.
- `threshold_analysis.png` — figura 4-painéis (usada como Fig. 3 no paper).

**Conclusão:** d=0.3 é o threshold mais alto com silhouette claramente positivo (0.0375), coerência interna acima de 60% e número de clusters estável e interpretável (11 clusters, 224 membros).

---

## Etapa 7 — Análise dos clusters finais

**Arquivo principal:** `cluster_statistics_final_11_clusters.xlsx`
**Conteúdo:** Para cada um dos 11 clusters: sítio, PDB_ID, superfamília CATH, código SCOPe, geometria de coordenação, número de coordenação.

**Verificação de integridade do xlsx (realizada em 2026-04-13):**

| Verificação | Resultado |
|---|---|
| Total de sítios | 224 — correto |
| Sites por cluster | 45, 34, 12, 28, 19, 19, 15, 15, 16, 11, 10 |
| Sites duplicados | Nenhum |
| Sites em múltiplos clusters | Nenhum |
| Formato `pdbid_N_chain` | Todos os 224 válidos |
| Formato CATH (`N.N.N.N`) | Todos válidos |
| CN fora do intervalo [1,12] | Nenhum |
| CN > 1 sem Geometry | Nenhum |
| CN = 1 com Geometry | Nenhum (45 sítios CN=1, todos sem geometry — esperado) |
| Cross-check com Stage2 | 100% dos 224 sítios presentes; cluster = grupo correto |

**Observações sobre os dados:**
- 45 sítios têm CN=1 e geometry ausente — comportamento esperado do MetalPDB.
- SCOPe ausente em 184/224 sítios (82%) — limitação de cobertura do SCOPe, não erro.
- CATH ausente em 2/224 sítios (0.9%).
- 1 sítio com geometria incomum: `5b66_74_AEFV` (Cluster 1, CN=7, hexagonal bipyramid with a vacancy) — válido pelo esquema MetalPDB.

**Script auxiliar:** `conta_dist=1.py`
**Propósito:** Investigação da discrepância entre Table 4 (valores incorretos) e threshold_analysis.py.
**Conclusão:** Apenas 1 par dos 2767 tem d=1.0 nos clusters finais — descartada a hipótese de que pares sem alinhamento distorciam o cálculo. A causa real era média simples vs. ponderada (ver Etapa 6).

---

## Etapa 8 — Correspondência com o motivo EF-hand

**Script:** `efhand_check.py`
**Propósito:** Verificar se algum dos 11 clusters corresponde ao motivo EF-hand, em resposta ao comentário do Revisor 2.

**Estratégia em 3 camadas (em ordem de custo):**

| Camada | Fonte | Critério |
|--------|-------|---------|
| 1ª | SCOPe (já no xlsx) | código começa com `a.7.` (EF-hand fold) |
| 2ª | CATH (já no xlsx) | contém `1.10.238` (EF-hand superfamily) |
| 3ª | API PDBe/Pfam | presença de PF00036, PF13499, PF13405 ou PF13833 |

**Saídas:**
- `efhand_per_site.csv` — flag EF-hand e fonte para cada um dos 224 sítios.
- `efhand_summary.csv` — contagem e % de membros EF-hand por cluster.

**Resultado:**
Apenas **1 proteína EF-hand** identificada no dataset inteiro: **2bko** (*Pyrococcus horikoshii*, proteína hipotética PH0236), com dois sítios de ligação em clusters distintos:

| Sítio | Cluster | CN | Geometria | SCOPe | CATH |
|-------|---------|-----|-----------|-------|------|
| `2bko_5_A` | Cluster 2 | 7 | pentagonal bipyramid (regular) | a.7.12;d.286.1 | 1.20.58.220;3.30.70.1450 |
| `2bko_1_A` | Cluster 3 | 4 | irregular | a.7.12;d.286.1 | 1.20.58.220;3.30.70.1450 |

**CORREÇÃO IMPORTANTE — SCOPe a.7.12:**
SCOPe a.7.12 = **PhoU-like** (superfamília 109755), NÃO calmodulin-like. A 2bko é uma proteína hipotética de canal de potássio putativo, não calmodulina. Isso foi verificado diretamente nos dados SCOPe pelo usuário. O texto anterior do PIPELINE e do paper estava errado nesse ponto.

**Parágrafo final para o paper (Seção 4.5, 3º parágrafo) — com \cite{}:**

```latex
No cluster was found to be exclusively composed of EF-hand motifs. The only
protein classified in the EF-hand fold (SCOPe a.7.12, PhoU-like superfamily)
\cite{Chandonia2022SCOPe} present in the dataset is \textit{Pyrococcus
horikoshii} hypothetical protein 2bko, which contributes two structurally
distinct calcium-binding sites to different clusters: site 2bko\_5\_A
(Cluster~2, coordination number~7, pentagonal bipyramidal geometry) displays
the canonical EF-hand calcium coordination arrangement \cite{Kirberger2008},
while site 2bko\_1\_A (Cluster~3, coordination number~4, irregular geometry)
represents an atypical secondary site. This demonstrates that the clustering
captures local coordination geometry rather than global protein function, and
that the structural diversity of calcium-binding MFSs extends beyond classical
EF-hand architectures.
```

---

## Etapa 9 — Composição de folds por cluster (nova)

**Script:** `fold_composition.py` (criado nesta sessão)
**Propósito:** Detalhar quais superfamílias CATH e SCOPe estão presentes em cada cluster, em resposta ao Revisor.

**O que faz:**
1. Carrega o xlsx com o mesmo loader de `efhand_check.py`.
2. Para cada cluster, conta quantos sítios têm cada superfamília (trata códigos múltiplos separados por `;`).
3. Gera `fold_composition.csv` (uma linha por cluster/fonte/superfamília) e `fold_composition.txt` (relatório legível).

**Saídas:**
- `fold_composition.csv`
- `fold_composition.txt`

**Resultados principais:**
- 230 superfamílias CATH distintas em 224 sítios (mais superfamílias do que sítios — muitos têm múltiplas anotações).
- Nenhuma superfamília contribui com mais de 3 sítios para qualquer cluster.
- Superfamílias que aparecem em mais clusters:

| CATH | Nome verificado | Clusters |
|------|----------------|---------|
| 2.60.40.10 | Immunoglobulins | 5 clusters |
| 3.20.20.140 | Metal-dependent hydrolases | 4 clusters |
| 3.40.50.720 | NAD(P)-binding Rossmann-like Domain | 4 clusters |

**CORREÇÕES de nomenclatura CATH verificadas no site CATH (cath.info):**
- `2.60.40.10` = **Immunoglobulins** — NÃO "jelly-roll" (erro corrigido)
- `3.20.20.80` = **Glycosidases** — NÃO "TIM barrel" (erro corrigido; removido do parágrafo)
- `3.40.50.720` = **NAD(P)-binding Rossmann-like Domain** — correto

**Nova Seção 4.6 para o paper — com \cite{}:**

```latex
\subsection{Fold Composition of the Structural Clusters}

The fold composition of the 11 clusters reveals high structural heterogeneity:
across the 224 clustered sites, 230 distinct CATH superfamilies
\cite{Orengo1997CATH} were identified, and no single superfamily contributes
more than three sites to any cluster. Several superfamilies recur across
multiple clusters — notably immunoglobulin-like domains (CATH 2.60.40.10,
present in 5~clusters), metal-dependent hydrolases (CATH 3.20.20.140,
4~clusters), and NAD(P)-binding Rossmann-like domains (CATH 3.40.50.720,
4~clusters) \cite{Orengo1997CATH} — indicating that the same global fold can
support calcium-binding sites with distinct local coordination geometries. This
pattern is consistent with the view that structurally similar metal-binding
microenvironments arise independently in proteins with different structural
frameworks, as previously demonstrated for zinc-binding sites
\cite{Andreini2011}.
```

**Frase adicional para Seção 5.2 (Discussion):**

```latex
The fold composition analysis further showed that 230 distinct CATH
superfamilies are represented across 224 sites, with no superfamily dominating
any single cluster, reinforcing the conclusion that convergent evolution
produces similar calcium-binding microenvironments in unrelated protein folds
\cite{Andreini2011}.
```

---

## Estado atual do paper (2026-04-13)

**Arquivo:** `iwbbio2026_sabrina.pdf` (versão mais recente)

**Referências BibTeX relevantes (`references.bib`):**

| Chave BibTeX | Ref no paper | O que é | Onde é citada |
|---|---|---|---|
| `Andreini2011` | [1] | Metodologia original para zinco | Methods (justificativa pipeline), Seções 4.6 e 5.2 |
| `Chandonia2022SCOPe` | [4] | SCOPe database | Seção 4.5 (EF-hand, classificação a.7.12) |
| `Kirberger2008` | [9] | Coordenação canônica EF-hand, CN=7, geometria pentagonal bipyramidal | Introdução p.2, Seção 4.5 |
| `Orengo1997CATH` | [12] | CATH database | Seção 4.6 (fold composition) |
| `Mukherjee2009MMalign` | [11] | MM-align — alinhamento de complexos multi-chain | Seção 3 (Methods), justificativa da ferramenta |
| `rousseeuw1987silhouettes` | [17] | Silhouette score | Seção 3 / Table 3 |
| `pedregosa2011scikit` | — | Scikit-learn — usado para calcular silhouette | Methods / implementação |
| `Putignano2017` | — | MetalPDB — banco de dados de origem das estruturas MFS | Introdução / Methods |
| `Xu2010TMscoreSignificance` | — | TM-score = 0.5 tem significância estrutural — justifica o threshold da Stage 1 | Methods (threshold 0.5 Stage 1) |
| `Valasatava2015` | — | Relações ocultas entre metaloproteínas via comparação de metal sites — trabalho relacionado | Related Work |

### Já corrigido:
- Seção 3.8: algoritmo descrito corretamente como greedy nearest-neighbour (não single-linkage).
- Table 3: valores corretos de threshold_analysis.py para os 9 thresholds.
- Fig. 3: figura de 4-painéis do threshold_analysis.py.
- Silhouette com citação Rousseeuw [17].

### Pendente:
- [ ] **Caption da Table 3** — ainda como "Placeholder Caption" no LaTeX — substituir.
- [ ] **Caption da Fig. 3** — ainda como "Enter Caption" no LaTeX — substituir.
- [ ] **Referência [?] na p.9** — deve ser [11] (Mukherjee & Zhang, 2009, MM-align).
- [ ] **Seção 4.5** — inserir o parágrafo EF-hand (Etapa 8 acima) como 3º parágrafo da seção existente.
- [ ] **Nova Seção 4.6** — inserir a seção de fold composition (Etapa 9 acima) após a 4.5.
- [ ] **Seção 5.2** — inserir a frase adicional sobre fold composition + EF-hand ao final do parágrafo relevante.
- [ ] **SCOPe c.37.1** — verificar nome manualmente em scop.berkeley.edu (WebFetch bloqueado por SSL durante a sessão).
- [ ] **Table 4 removida** — confirmar se foi removida no LaTeX; se não, adicionar frase explicativa sobre pct_below_0.5 ser média ponderada.
- [ ] **Análise de sensibilidade para raio de 5Å** — avaliado como inviável a curto prazo (exigiria reprocessar o pipeline do zero).

### Texto sugerido para substituir Table 4 (se ainda existir):
> The percentage of intra-cluster pairwise distances below 0.5 — computed as a weighted average over all valid clusters, where each cluster contributes proportionally to its number of pairwise comparisons — decreased from 100% at d = 0.2 to 60.46% at d = 0.3...

---

## Diferenças metodológicas em relação a Andreini et al. (2011)

| Aspecto | Andreini 2011 (zinco) | Este trabalho (cálcio) | Justificativa |
|---|---|---|---|
| Ferramenta de alinhamento | FAST | MM-align | MM-align suporta multi-chain; mais adequado para MFS com múltiplas cadeias |
| Algoritmo Stage 2 | Single-linkage | Greedy NN | Evita chaining effect |
| Threshold Stage 2 | Variável | d < 0.3 fixo | Escolhido por silhouette + coerência interna |
| Pseudo-clusters | Sim (7 clusters, 16% dos sítios) | Não | Clusters < 10 membros descartados |

---

## Referências-chave

- **Andreini et al. (2011):** metodologia original para zinco que este trabalho adapta para cálcio.
- **Mukherjee & Zhang (2009):** TM-score ≥ 0.7 → >98% dos complexos proteicos compartilham função biológica idêntica. Citado para justificar d=0.3.
- **Rousseeuw (1987):** definição do silhouette score.
- **MetalPDB:** banco de dados de origem das estruturas MFS.
- **Kirberger (2008):** coordenação canônica EF-hand — CN=7, geometria pentagonal bipyramidal.
- **Orengo et al. (1997):** CATH database — nomes e classificação de superfamílias.
- **Chandonia et al. (2022):** SCOPe database — classificação a.7.12 = PhoU-like.
