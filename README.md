# Structural Classification and Clustering of Calcium-Binding Sites in Proteins Using Minimal Functional Sites

Paper accepted at IWBBIO 2026.

## Abstract

We present a two-stage structural clustering pipeline for calcium-binding sites represented as Minimal Functional Sites (MFS). In Stage 1, sites are grouped within CATH and SCOPe superfamilies using TM-score-based similarity matrices, reducing redundancy and selecting high-resolution representatives. In Stage 2, a global nearest-neighbor clustering of 2,960 representatives produces 11 structurally recurrent calcium-binding site architectures, revealing conserved coordination geometries across diverse protein folds.

## Pipeline

See [PIPELINE_ENGLISH.md](PIPELINE_ENGLISH.md) for full details.

## Requirements

```bash
pip install -r requirements.txt
```

## How to Reproduce

All scripts must be run from the **project root** directory.

**Stage 0 — Data collection from MetalPDB:**

```bash
python scripts/stage1/APIMetalPDB.py
```

**Stage 1 — Intra-superfamily clustering:**

```bash
python scripts/stage1/similarity_matrix_and_clustering.py
python scripts/stage1/APIresolutionsCATH.py
python scripts/stage1/APIresolutionsSCOPe.py
```

**Stage 2 — Global clustering:**

```bash
python scripts/stage2/global_clustering_tm07.py
python scripts/stage2/threshold_analysis.py
```

**Complementary analyses:**

```bash
python scripts/analysis/efhand_check.py
python scripts/analysis/fold_composition.py
python scripts/analysis/efhand_stage1_analysis.py
```

## Results

Main results are in `results/`:

- [`results/clustering/cluster_statistics_final_11_clusters.xlsx`](results/clustering/cluster_statistics_final_11_clusters.xlsx) — membership and statistics for the 11 final clusters
- [`results/analysis/threshold_analysis.png`](results/analysis/threshold_analysis.png) — effect of distance threshold on clustering
- [`results/analysis/fold_composition.csv`](results/analysis/fold_composition.csv) — CATH/SCOPe superfamily breakdown per cluster
- [`results/analysis/efhand_summary.csv`](results/analysis/efhand_summary.csv) — EF-hand presence per cluster

## Supplementary Data

Large files are available in [Release v1.0](https://github.com/gabrielaog/calcium-binding-site-structural-clustering/releases/tag/v1.0):

- `DistanceMatrix.csv` — pairwise distance matrix for 2,960 Stage-1 representatives
- `representativeSites.zip` — PDB structures of the representative sites
- `similarity_matrices_CATH.xlsx` / `similarity_matrices_SCOPe.xlsx` — Stage 1 similarity matrices

See [`data/README_data.md`](data/README_data.md) for download instructions.

## Authors

Gabriela Dias, Vinícius Paiva, Sandro Izidoro, Cláudia Andreini, Sabrina Silveira  
UFV / UNIFEI / University of Florence

## Citation

BibTeX will be added upon publication.
