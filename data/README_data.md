# Supplementary Data

The large data files are not stored in this repository due to GitHub's size limits.
They are available as assets in **[Release v1.0](https://github.com/gabrielaog/calcium-binding-sites-structural-clustering/releases/tag/v1.0)** — download each file directly from there.

## Files

| File | Size | Description |
|------|------|-------------|
| `DistanceMatrix.csv` | ~70 MB | Pairwise distance matrix (d = 1 − avg TM-score) for all 2,960 Stage-1 representatives |
| `representativeSites.zip` | ~11 MB | PDB structures of the Stage-1 representative sites |
| `similarity_matrices_CATH.xlsx` | ~128 MB | TM-score similarity matrices for each CATH superfamily (Stage 1) |
| `similarity_matrices_SCOPe.xlsx` | ~105 MB | TM-score similarity matrices for each SCOPe superfamily (Stage 1) |

## External Data Sources

The following files are from public databases and must be downloaded directly from their sources:

| File                           | Source    | URL                                                                |
|--------------------------------|-----------|--------------------------------------------------------------------|
| `13jan-cath-domain-list.txt`   | CATH v4.3 | [cathdb.info](https://www.cathdb.info/wiki/doku/?id=release_notes) |
| `dir.des.scope.2.08-stable.txt`| SCOPe 2.08| [scop.berkeley.edu](https://scop.berkeley.edu/downloads/)          |

## Usage

Place `DistanceMatrix.csv` in the project root before running any Stage 2 scripts.
Place the similarity matrix `.xlsx` files in the project root before running Stage 1 scripts.
