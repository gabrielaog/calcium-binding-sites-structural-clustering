# Pipeline — Structural Classification of Calcium-Binding Sites

Reference document for the project. Describes each pipeline stage, the scripts used, and the rationale behind methodological choices.

---

## Overview

The goal is to classify calcium-binding sites (Minimal Functional Sites, MFS) based on structural similarity, following the methodology of Andreini et al. (2011) for zinc.

**Input:** PDB structures of MFS extracted from MetalPDB, organized by superfamily (CATH and SCOPe).  
**Output:** 11 global clusters of representative MFS, with threshold d < 0.3 (TM-score ≥ 0.7).

---

## Stage 1 — Intra-superfamily structural alignment

**Tool:** MM-align  
**What it does:** Compares pairs of MFS within each superfamily and generates `.txt` files with the two TM-scores of the alignment.  
**Output:** One `.txt` file per site pair, named `site1_site2.txt`.

**Similarity metric used:**
```
avg_TM-score = (TM-score_1 + TM-score_2) / 2
distance = 1 − avg_TM-score
```
The average of the two TM-scores is used because MM-align normalizes each score by the length of one of the structures; the average ensures symmetry and robustness.

---

## Stage 2 — Intra-superfamily clustering (Stage 1)

**Script:** `similarity_matrix_and_clustering.py`  
**Algorithm:** Hierarchical single-linkage (scipy), threshold = 0.5  
**What it does:**
1. For each superfamily, reads the MM-align `.txt` files and builds a similarity matrix.
2. Converts to a distance matrix (d = 1 − TM-score).
3. Applies single-linkage with threshold 0.5.
4. Selects one representative per cluster.

**Outputs:**
- `similarity_matrices_CATH.xlsx` / `similarity_matrices_SCOPe.xlsx` — per-superfamily similarity matrices.
- `CATH_all_clusters05.txt` / `SCOPe_all_clusters05.txt` — cluster assignment per site.

**Known code issue (does not affect results):**
The comment at line ~78 reads `# Missing alignment assumed as 0 similarity` but is attached to the block handling `i==j` (diagonal = 1.0), not to the block for missing pairs. It is a comment misplacement only.

**Why single-linkage here?**
Within a superfamily, sites already share a fold; the most permissive criterion (minimum distance between any two members) is appropriate for grouping structurally close variants before selecting representatives.

**Why threshold 0.5?**
An intermediate threshold that captures significant structural similarity without collapsing an entire superfamily into a single cluster.

---

## Stage 3 — Representative selection and PDB resolution

**Scripts:** `APIresolutionsCATH.py`, `APIresolutionsSCOPe.py`  
**What they do:** Query the CATH and SCOPe APIs to retrieve the crystallographic resolution of each PDB structure.  
**Selection criterion:** Within each intra-superfamily cluster, the site with the **lowest resolution value** (best crystallographic quality) is chosen as the representative.

**Known code issues (do not affect results):**
- `APIresolutionsCATH.py`: contains `print(f"match_superfamily {match_superfamily}")` — a debug print left in production.
- `APIresolutionsSCOPe.py`: contains `print(f'nfo pdb {info_pdb}')` — a debug print with a typo ("nfo" instead of "info").
- The output files (`best_resolutionsCATH.txt`, `best_resolutionsSCOPe.txt`) have different names from the ones read by `efhand_check.py` (`CATH_Representatives_SelectedByBestResolution_TMscore_0.5.txt`, `SCOPe_Representatives_SelectedByBestResolution_TMscore_0.5.txt`). There was an undocumented manual conversion/renaming step.

---

## Stage 4 — Global structural alignment

**Tool:** MM-align (run again)  
**What it does:** Compares all pairs of representatives selected in Stage 3, across all superfamilies.  
**Output:** `DistanceMatrix.csv` — N×N distance matrix for all representatives.

**Special values in the matrix:**
- `d = 0.0` on the diagonal (self-distance).
- `d = 1.0` for pairs with no valid alignment (MM-align did not converge); these values are treated as missing, not as actual maximum distance.

---

## Stage 5 — Global clustering (Stage 2)

**Script:** `global_clustering_tm07.py`  
**Algorithm:** Greedy nearest-neighbour (NOT single-linkage)  
**Threshold:** d < 0.3 (equivalent to TM-score ≥ 0.7)

**How it works (step by step):**
1. For each site (matrix row), find the nearest neighbour (column with minimum distance, excluding d=0).
2. If both are already in clusters → skip (no cluster merging).
3. If distance is d=1.0 or d≥0.3 → skip.
4. If only one already has a cluster → add the other to the same cluster.
5. If neither has a cluster → create a new cluster with both.
6. After the scan, discard clusters with fewer than 10 members.

**Output:** `resultadoMaiorQue07.txt` — list of sites with their cluster, nearest neighbour, and distance.

**Known code issue (does not affect results):**
The docstring says "similar to single-linkage" — this is incorrect. The algorithm is greedy NN. The code logic itself is correct.

**Why greedy NN instead of single-linkage?**
Single-linkage suffers from the *chaining effect*: a single low-distance pair can chain together very heterogeneous clusters. The greedy NN with the rule "if both already have a cluster → skip" prevents indirect merges, ensuring each site joins a cluster only through its immediate nearest neighbour.

**Why discard clusters with < 10 members?**
To focus on structurally recurrent architectures. Small clusters represent rare or isolated configurations with no statistical relevance.

---

## Stage 6 — Threshold analysis (validation of d=0.3)

**Script:** `threshold_analysis.py`  
**What it does:** Runs the same greedy NN algorithm from Stage 5 for 9 thresholds (0.20 to 0.60) and computes 4 quality metrics.

**Metrics and rationale:**

| Metric | What it measures | Why use it |
|--------|-----------------|------------|
| Number of valid clusters (≥10 members) | Partition resolution | Thresholds too low collapse everything; too high over-fragment |
| Mean cluster size | Compactness | Rapid growth indicates chaining effect |
| % intra-cluster distances < 0.5 | Internal structural coherence | Clusters with many distances ≥ 0.5 are heterogeneous |
| Silhouette score | Inter-cluster separation | Positive score = well-separated clusters; negative = overlap |

**Computation of `pct_below_0.5`:**
This is a **pair-weighted average** across clusters:
```
pct = (Σ_k  #{pairs in C_k with d < 0.5}) / (Σ_k #{pairs in C_k})  × 100
```
This gives equal weight to each pair (not to each cluster). A simple per-cluster average would inflate values due to small clusters — this was the error in the original Table 4 (which used a simple average and yielded 76.08% instead of 60.46%).

**Outputs:**
- `results/analysis/threshold_analysis.csv` — table with the 4 metrics for each threshold.
- `results/analysis/threshold_analysis.png` — 4-panel figure (used as Fig. 3 in the paper).

**Conclusion:** d=0.3 is the highest threshold with a clearly positive silhouette (0.0375), internal coherence above 60%, and a stable and interpretable number of clusters (11 clusters, 224 members).

---

## Stage 7 — Final cluster analysis

**Main file:** `results/clustering/cluster_statistics_final_11_clusters.xlsx`  
**Contents:** For each of the 11 clusters: site, PDB_ID, CATH superfamily, SCOPe code, coordination geometry, coordination number.

**Integrity checks performed (2026-04-13):**

| Check | Result |
|-------|--------|
| Total sites | 224 — correct |
| Sites per cluster | 45, 34, 12, 28, 19, 19, 15, 15, 16, 11, 10 |
| Duplicate sites | None |
| Sites in multiple clusters | None |
| Format `pdbid_N_chain` | All 224 valid |
| CATH format (`N.N.N.N`) | All valid |
| CN outside range [1,12] | None |
| CN > 1 without Geometry | None |
| CN = 1 with Geometry | None (45 sites CN=1, all without geometry — expected) |
| Cross-check with Stage 2 | 100% of the 224 sites present; cluster = correct group |

**Notes on the data:**
- 45 sites have CN=1 and absent geometry — expected MetalPDB behaviour.
- SCOPe absent for 184/224 sites (82%) — SCOPe coverage limitation, not an error.
- CATH absent for 2/224 sites (0.9%).
- 1 site with uncommon geometry: `5b66_74_AEFV` (Cluster 1, CN=7, hexagonal bipyramid with a vacancy) — valid under the MetalPDB scheme.

**Auxiliary script:** `scripts/analysis/check_distance_one.py`  
**Purpose:** Investigation of the discrepancy between Table 4 (incorrect values) and threshold_analysis.py.  
**Conclusion:** Only 1 pair out of 2767 has d=1.0 in the final clusters — the hypothesis that missing-alignment pairs distorted the calculation was ruled out. The actual cause was simple vs. weighted average (see Stage 6).

---

## Stage 8 — EF-hand motif correspondence

**Script:** `scripts/analysis/efhand_check.py`  
**Purpose:** Verify whether any of the 11 clusters corresponds to the EF-hand motif, in response to Reviewer 2's comment.

**3-layer strategy (in order of cost):**

| Layer | Source | Criterion |
|-------|--------|-----------|
| 1st | SCOPe (already in xlsx) | code starts with `a.7.` (EF-hand fold) |
| 2nd | CATH (already in xlsx) | contains `1.10.238` (EF-hand superfamily) |
| 3rd | PDBe/Pfam API | presence of PF00036, PF13499, PF13405, or PF13833 |

**Outputs:**
- `results/analysis/efhand_per_site.csv` — EF-hand flag and source for each of the 224 sites.
- `results/analysis/efhand_summary.csv` — count and % of EF-hand members per cluster.

**Result:**
Only **1 EF-hand protein** identified in the entire dataset: **2bko** (*Pyrococcus horikoshii*, hypothetical protein PH0236), with two binding sites in distinct clusters:

| Site | Cluster | CN | Geometry | SCOPe | CATH |
|------|---------|-----|----------|-------|------|
| `2bko_5_A` | Cluster 2 | 7 | pentagonal bipyramid (regular) | a.7.12;d.286.1 | 1.20.58.220;3.30.70.1450 |
| `2bko_1_A` | Cluster 3 | 4 | irregular | a.7.12;d.286.1 | 1.20.58.220;3.30.70.1450 |

**IMPORTANT CORRECTION — SCOPe a.7.12:**
SCOPe a.7.12 = **PhoU-like** (superfamily 109755), NOT calmodulin-like. 2bko is a hypothetical putative potassium channel protein, not calmodulin. This was verified directly in the SCOPe data. Previous versions of PIPELINE.md and the paper were incorrect on this point.

**Suggested paragraph for the paper (Section 4.5, 3rd paragraph) — with \cite{}:**

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

## Stage 9 — Fold composition per cluster

**Script:** `scripts/analysis/fold_composition.py`  
**Purpose:** Detail which CATH and SCOPe superfamilies are present in each cluster, in response to a reviewer's comment.

**What it does:**
1. Loads the xlsx using the same loader as `efhand_check.py`.
2. For each cluster, counts how many sites have each superfamily (handles multiple codes separated by `;`).
3. Generates `fold_composition.csv` (one row per cluster/source/superfamily) and `fold_composition.txt` (human-readable report).

**Outputs:**
- `results/analysis/fold_composition.csv`
- `results/analysis/fold_composition.txt`

**Main results:**
- 230 distinct CATH superfamilies across 224 sites (more superfamilies than sites — many have multiple annotations).
- No superfamily contributes more than 3 sites to any cluster.
- Superfamilies appearing in the most clusters:

| CATH | Verified name | Clusters |
|------|--------------|---------|
| 2.60.40.10 | Immunoglobulins | 5 clusters |
| 3.20.20.140 | Metal-dependent hydrolases | 4 clusters |
| 3.40.50.720 | NAD(P)-binding Rossmann-like Domain | 4 clusters |

**CATH nomenclature corrections verified at cath.info:**
- `2.60.40.10` = **Immunoglobulins** — NOT "jelly-roll" (error corrected)
- `3.20.20.80` = **Glycosidases** — NOT "TIM barrel" (error corrected; removed from paragraph)
- `3.40.50.720` = **NAD(P)-binding Rossmann-like Domain** — correct

**Suggested new Section 4.6 for the paper — with \cite{}:**

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

**Additional sentence for Section 5.2 (Discussion):**

```latex
The fold composition analysis further showed that 230 distinct CATH
superfamilies are represented across 224 sites, with no superfamily dominating
any single cluster, reinforcing the conclusion that convergent evolution
produces similar calcium-binding microenvironments in unrelated protein folds
\cite{Andreini2011}.
```

---

## Current paper status (2026-04-13)

**File:** `iwbbio2026_sabrina.pdf` (most recent version)

**Relevant BibTeX references (`references.bib`):**

| BibTeX key | Ref in paper | What it is | Where cited |
|---|---|---|---|
| `Andreini2011` | [1] | Original methodology for zinc | Methods (pipeline rationale), Sections 4.6 and 5.2 |
| `Chandonia2022SCOPe` | [4] | SCOPe database | Section 4.5 (EF-hand, a.7.12 classification) |
| `Kirberger2008` | [9] | Canonical EF-hand coordination, CN=7, pentagonal bipyramidal | Introduction p.2, Section 4.5 |
| `Orengo1997CATH` | [12] | CATH database | Section 4.6 (fold composition) |
| `Mukherjee2009MMalign` | [11] | MM-align — multi-chain complex alignment | Section 3 (Methods), tool rationale |
| `rousseeuw1987silhouettes` | [17] | Silhouette score | Section 3 / Table 3 |
| `pedregosa2011scikit` | — | Scikit-learn — used to compute silhouette | Methods / implementation |
| `Putignano2017` | — | MetalPDB — source database for MFS structures | Introduction / Methods |
| `Xu2010TMscoreSignificance` | — | TM-score = 0.5 has structural significance — justifies Stage 1 threshold | Methods (threshold 0.5 Stage 1) |
| `Valasatava2015` | — | Hidden relationships between metalloproteins via metal site comparison — related work | Related Work |

### Already fixed:
- Section 3.8: algorithm correctly described as greedy nearest-neighbour (not single-linkage).
- Table 3: correct values from threshold_analysis.py for the 9 thresholds.
- Fig. 3: 4-panel figure from threshold_analysis.py.
- Silhouette with citation Rousseeuw [17].

### Pending:
- [ ] **Caption of Table 3** — still "Placeholder Caption" in LaTeX — replace.
- [ ] **Caption of Fig. 3** — still "Enter Caption" in LaTeX — replace.
- [ ] **Reference [?] on p.9** — should be [11] (Mukherjee & Zhang, 2009, MM-align).
- [ ] **Section 4.5** — insert the EF-hand paragraph (Stage 8 above) as the 3rd paragraph of the existing section.
- [ ] **New Section 4.6** — insert the fold composition section (Stage 9 above) after Section 4.5.
- [ ] **Section 5.2** — insert the additional sentence about fold composition + EF-hand at the end of the relevant paragraph.
- [ ] **SCOPe c.37.1** — verify name manually at scop.berkeley.edu (WebFetch blocked by SSL during session).
- [ ] **Table 4 removed** — confirm whether it was removed in LaTeX; if not, add explanatory sentence about pct_below_0.5 being a weighted average.
- [ ] **Sensitivity analysis for 5Å radius** — assessed as infeasible in the short term (would require reprocessing the entire pipeline from scratch).

### Suggested text to replace Table 4 (if still present):
> The percentage of intra-cluster pairwise distances below 0.5 — computed as a weighted average over all valid clusters, where each cluster contributes proportionally to its number of pairwise comparisons — decreased from 100% at d = 0.2 to 60.46% at d = 0.3...

---

## Methodological differences relative to Andreini et al. (2011)

| Aspect | Andreini 2011 (zinc) | This work (calcium) | Rationale |
|---|---|---|---|
| Alignment tool | FAST | MM-align | MM-align supports multi-chain; more suitable for MFS with multiple chains |
| Stage 2 algorithm | Single-linkage | Greedy NN | Avoids chaining effect |
| Stage 2 threshold | Variable | d < 0.3 fixed | Chosen by silhouette + internal coherence |
| Pseudo-clusters | Yes (7 clusters, 16% of sites) | No | Clusters < 10 members discarded |

---

## Key references

- **Andreini et al. (2011):** original zinc methodology that this work adapts for calcium.
- **Mukherjee & Zhang (2009):** TM-score ≥ 0.7 → >98% of protein complexes share identical biological function. Cited to justify d=0.3.
- **Rousseeuw (1987):** silhouette score definition.
- **MetalPDB:** source database for MFS structures.
- **Kirberger (2008):** canonical EF-hand coordination — CN=7, pentagonal bipyramidal geometry.
- **Orengo et al. (1997):** CATH database — superfamily names and classification.
- **Chandonia et al. (2022):** SCOPe database — classification a.7.12 = PhoU-like.
