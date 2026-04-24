"""
EF-hand identification — FULL PIPELINE (Stage 1 → Stage 2 → Final clusters).

Checks for EF-hand at three levels:

  A. Stage 1 (intra-superfamily clustering)
       Reads CATH_StructuralClusters_TMscore_0.5.txt and
       SCOPe_StructuralClusters_TMscore_0.5.txt to count how many
       EF-hand sites were present before the Stage-1 compression and
       which representative was selected for each EF-hand cluster.

  B. Stage 2 (global clustering)
       Reads SecondStage_Clustering_MembershipList_0.3.txt to report
       whether each EF-hand representative ended up clustered or
       unclustered in the final 11-cluster solution.

  C. Final 11 clusters (3-layer check)
       Reads cluster_statistics_final_11_clusters.xlsx and checks every
       one of the 224 final sites via:
         1. SCOPe code starts with 'a.7.'   (Calmodulin-like / EF-hand)
         2. CATH  code starts with '1.10.238' (EF-hand superfamily)
         3. PDBe Pfam API: PF00036 / PF13499 / PF13405 / PF13833

       Note: the CATH codes stored in the xlsx are domain-level annotations
       from MetalPDB and may differ from the superfamily code 1.10.238.10,
       so the SCOPe and Pfam layers are the primary detectors for section C.

Output files:
  efhand_stage1_report.txt  — Stage 1 EF-hand counts, cluster sizes,
                              representatives
  efhand_stage2_fate.csv    — EF-hand representative → Stage-2 outcome
                              (clustered in Group N, or unclustered)
  efhand_per_site.csv       — one row per site in the 224 final-cluster
                              members, with EF-hand flag and source
  efhand_summary.csv        — per-cluster count and % EF-hand in the
                              final 11 clusters
"""

import re
import ast
import time

import pandas as pd
import requests
from collections import defaultdict

# ── file paths ─────────────────────────────────────────────────────────────
XLSX              = "results/clustering/cluster_statistics_final_11_clusters.xlsx"
CATH_CLUSTERS     = "results/clustering/CATH_StructuralClusters_TMscore_0.5.txt"
SCOPE_CLUSTERS    = "results/clustering/SCOPe_StructuralClusters_TMscore_0.5.txt"
CATH_REPS         = "results/clustering/CATH_Representatives_SelectedByBestResolution_TMscore_0.5.txt"
SCOPE_REPS        = "results/clustering/SCOPe_Representatives_SelectedByBestResolution_TMscore_0.5.txt"
STAGE2_MEMBERSHIP = "results/clustering/SecondStage_Clustering_MembershipList_0.3.txt"

EF_CATH_PREFIX  = "1.10.238"
EF_SCOPE_PREFIX = "a.7."
EF_PFAM         = {"PF00036", "PF13499", "PF13405", "PF13833"}


# =============================================================================
# A. STAGE 1 — EF-hand sites before compression
# =============================================================================

def parse_stage1_clusters(filepath):
    """
    Read a Stage-1 cluster file and return:
        { superfamily: { cluster_id (int): [site, ...] } }
    File format:
        Clusters {superfamily} threshold_0.5:
        {site}: {cluster_id}
        ...
    """
    result = {}
    current_sf = None
    with open(filepath, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("Clusters "):
                current_sf = line.split()[1]
                result[current_sf] = defaultdict(list)
            elif current_sf and ":" in line:
                site, cid_str = line.split(":", 1)
                result[current_sf][int(cid_str.strip())].append(site.strip())
    return result


def parse_stage1_representatives(filepath):
    """
    Read a Stage-1 representatives file and return:
        { superfamily: { cluster_id (int): representative_site } }
    Handles two line formats:
        Cluster representative: {site}
        Direct cluster representative: {site} - ({resolution} A)
    """
    reps = {}
    current_sf = None
    current_cid = None
    with open(filepath, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("Superfamily:"):
                current_sf = line.split(":", 1)[1].strip()
                reps.setdefault(current_sf, {})
            elif line.startswith("Cluster:"):
                current_cid = int(line.split(":", 1)[1].strip())
            elif "Cluster representative:" in line and current_sf and current_cid is not None:
                rep = line.split("Cluster representative:", 1)[1].strip()
                reps[current_sf][current_cid] = rep
            elif "Direct cluster representative:" in line and current_sf and current_cid is not None:
                rep = line.split("Direct cluster representative:", 1)[1].strip()
                rep = rep.split(" - ")[0].strip()
                reps[current_sf][current_cid] = rep
    return reps


def collect_stage1_efhand(cath_clusters, scope_clusters, cath_reps, scope_reps):
    """
    Return a list of dicts, one per EF-hand cluster found in Stage 1.
    """
    records = []
    for source, clusters, reps, prefix in [
        ("CATH",  cath_clusters,  cath_reps,  EF_CATH_PREFIX),
        ("SCOPe", scope_clusters, scope_reps, EF_SCOPE_PREFIX),
    ]:
        for sf, by_cid in clusters.items():
            if not sf.startswith(prefix):
                continue
            sf_reps = reps.get(sf, {})
            for cid, members in sorted(by_cid.items()):
                records.append({
                    "source":         source,
                    "superfamily":    sf,
                    "stage1_cluster": cid,
                    "n_members":      len(members),
                    "representative": sf_reps.get(cid, "—"),
                })
    return records


def write_stage1_report(records, filepath):
    total_sites    = sum(r["n_members"] for r in records)
    total_clusters = len(records)

    lines = [
        "=" * 70,
        "EF-HAND — STAGE 1  (intra-superfamily clustering, threshold 0.5)",
        "=" * 70,
        "",
        f"  EF-hand sites across all Stage-1 clusters : {total_sites}",
        f"  EF-hand Stage-1 clusters (CATH + SCOPe)   : {total_clusters}",
        "",
        f"  {'Source':<6}  {'Superfamily':<20}  {'Cluster':>7}  "
        f"{'Members':>8}  Representative",
        "  " + "-" * 66,
    ]
    for r in records:
        lines.append(
            f"  {r['source']:<6}  {r['superfamily']:<20}  "
            f"{r['stage1_cluster']:>7}  {r['n_members']:>8}  {r['representative']}"
        )
    lines += ["", "=" * 70]

    with open(filepath, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"Saved: {filepath}")
    return total_sites


# =============================================================================
# B. STAGE 2 — fate of EF-hand representatives
# =============================================================================

def parse_stage2_membership(filepath):
    """
    Read Stage-2 membership file and return { site: group_id (int) }.
    File format:
        Group {N}: ['{site}', ...]
    """
    membership = {}
    with open(filepath, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            m = re.match(r"Group\s+(\d+):\s*(\[.*\])", line)
            if not m:
                continue
            group_id = int(m.group(1))
            for site in ast.literal_eval(m.group(2)):
                membership[site] = group_id
    return membership


def build_stage2_fate(stage1_records, stage2_membership):
    """
    For each EF-hand Stage-1 representative, report its Stage-2 outcome.
    """
    rows = []
    for r in stage1_records:
        rep   = r["representative"]
        group = stage2_membership.get(rep)
        rows.append({
            "source":             r["source"],
            "superfamily":        r["superfamily"],
            "stage1_cluster":     r["stage1_cluster"],
            "n_stage1_members":   r["n_members"],
            "representative":     rep,
            "stage2_outcome":     f"Group {group}" if group is not None else "unclustered",
        })
    return rows


# =============================================================================
# C. FINAL 11 CLUSTERS — per-site EF-hand check
# =============================================================================

def load_final_clusters(xlsx_path):
    raw = pd.read_excel(xlsx_path, header=None)
    cluster_header_rows = []
    for i, row in raw.iterrows():
        val = str(row[0])
        if val.startswith("Cluster") and row[1] != "PDB_ID":
            cluster_header_rows.append((i, val.strip()))

    records = []
    for idx, (header_row, cluster_name) in enumerate(cluster_header_rows):
        col_row  = header_row + 1
        data_end = (cluster_header_rows[idx + 1][0]
                    if idx + 1 < len(cluster_header_rows) else len(raw))
        block = raw.iloc[col_row + 1 : data_end].copy()
        block.columns = raw.iloc[col_row].tolist()
        block = block.dropna(subset=["Site"])
        block = block[block["Site"] != "Site"]
        block["Cluster"] = cluster_name
        records.append(block)

    df = pd.concat(records, ignore_index=True)
    df["CATH_superfamily"]       = df["CATH_superfamily"].fillna("NA").astype(str).str.strip()
    df["SCOPe_superfamily_code"] = df["SCOPe_superfamily_code"].fillna("NA").astype(str).str.strip()
    df["PDB_ID"]                 = df["PDB_ID"].astype(str).str.strip().str.lower()
    return df


def scop_is_efhand(scop_str):
    if scop_str in ("NA", "nan", ""):
        return False
    return any(c.strip().startswith(EF_SCOPE_PREFIX) for c in scop_str.split(";"))


def cath_is_efhand(cath_str):
    if cath_str in ("NA", "nan", ""):
        return False
    return any(c.strip().startswith(EF_CATH_PREFIX) for c in cath_str.split(";"))


def check_pfam(pdb_id, cache):
    if pdb_id in cache:
        return cache[pdb_id]
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{pdb_id}"
    try:
        r      = requests.get(url, timeout=10)
        result = (r.status_code == 200 and
                  bool(EF_PFAM & set(r.json().get(pdb_id, {}).get("Pfam", {}).keys())))
    except Exception:
        result = False
    cache[pdb_id] = result
    return result


def classify_site(row):
    if row["efhand_scop"]:
        return True, "SCOPe"
    if row["efhand_cath"]:
        return True, "CATH"
    if row["efhand_pfam"]:
        return True, "Pfam"
    return False, "-"


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":

    # ── A. Stage 1 ────────────────────────────────────────────────────────────
    print("=" * 60)
    print("A. STAGE 1 — EF-hand in intra-superfamily clusters")
    print("=" * 60)

    cath_clusters  = parse_stage1_clusters(CATH_CLUSTERS)
    scope_clusters = parse_stage1_clusters(SCOPE_CLUSTERS)
    cath_reps      = parse_stage1_representatives(CATH_REPS)
    scope_reps     = parse_stage1_representatives(SCOPE_REPS)

    stage1_records  = collect_stage1_efhand(
        cath_clusters, scope_clusters, cath_reps, scope_reps
    )
    total_ef_s1 = write_stage1_report(stage1_records, "results/analysis/efhand_stage1_report.txt")

    total_s1 = (
        sum(len(m) for sf in cath_clusters.values()  for m in sf.values()) +
        sum(len(m) for sf in scope_clusters.values() for m in sf.values())
    )
    print(f"  EF-hand sites in Stage 1 : {total_ef_s1} "
          f"(out of {total_s1} total site-superfamily pairs)")
    print(f"  EF-hand Stage-1 clusters : {len(stage1_records)}")

    # ── B. Stage 2 fate ───────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("B. STAGE 2 — fate of EF-hand representatives")
    print("=" * 60)

    stage2_membership = parse_stage2_membership(STAGE2_MEMBERSHIP)
    fate_rows = build_stage2_fate(stage1_records, stage2_membership)
    fate_df   = pd.DataFrame(fate_rows)
    fate_df.to_csv("results/analysis/efhand_stage2_fate.csv", index=False)

    print(fate_df[["source", "superfamily", "n_stage1_members",
                   "representative", "stage2_outcome"]].to_string(index=False))
    n_unclustered = (fate_df["stage2_outcome"] == "unclustered").sum()
    print(f"\n  Unclustered in Stage 2 : {n_unclustered} / {len(fate_df)} "
          f"EF-hand representatives")
    print("Saved: efhand_stage2_fate.csv")

    # ── C. Final clusters ─────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("C. FINAL 11 CLUSTERS — per-site EF-hand check")
    print("=" * 60)

    df = load_final_clusters(XLSX)
    print(f"  Total sites loaded: {len(df)}")

    df["efhand_scop"] = df["SCOPe_superfamily_code"].apply(scop_is_efhand)
    df["efhand_cath"] = df["CATH_superfamily"].apply(cath_is_efhand)

    needs_api = df[~(df["efhand_scop"] | df["efhand_cath"])]["PDB_ID"].unique()
    print(f"  Sites needing Pfam API check: {len(needs_api)} unique PDB IDs")

    pfam_cache = {}
    for i, pdb_id in enumerate(needs_api):
        check_pfam(pdb_id, pfam_cache)
        if (i + 1) % 50 == 0 or (i + 1) == len(needs_api):
            print(f"    Pfam API: {i+1}/{len(needs_api)} done")
        time.sleep(0.05)

    already_detected = df["efhand_scop"] | df["efhand_cath"]
    df["efhand_pfam"] = df.apply(
        lambda row: (False if already_detected[row.name]
                     else pfam_cache.get(row["PDB_ID"], False)),
        axis=1,
    )

    df[["is_efhand", "efhand_source"]] = df.apply(
        classify_site, axis=1, result_type="expand"
    )

    site_out = df[[
        "Cluster", "Site", "PDB_ID",
        "CATH_superfamily", "SCOPe_superfamily_code",
        "Geometry", "Coordination_number",
        "is_efhand", "efhand_source",
    ]]
    site_out.to_csv("results/analysis/efhand_per_site.csv", index=False)
    print("Saved: results/analysis/efhand_per_site.csv")

    summary_rows = []
    for cluster, group in df.groupby("Cluster", sort=False):
        n_total  = len(group)
        n_ef     = int(group["is_efhand"].sum())
        sources  = (group.loc[group["is_efhand"], "efhand_source"]
                        .value_counts().to_dict())
        summary_rows.append({
            "Cluster":    cluster,
            "N_members":  n_total,
            "N_efhand":   n_ef,
            "pct_efhand": round(n_ef / n_total * 100, 1),
            "sources":    str(sources),
        })

    summary = pd.DataFrame(summary_rows)
    summary.to_csv("results/analysis/efhand_summary.csv", index=False)

    print("\n=== EF-hand summary — final 11 clusters ===")
    print(summary[["Cluster", "N_members", "N_efhand", "pct_efhand"]].to_string(index=False))
    print("Saved: results/analysis/efhand_summary.csv")

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("KEY FINDINGS")
    print("=" * 60)
    print(f"  EF-hand sites in full dataset (Stage 1)   : {total_ef_s1}")
    print(f"  Stage-1 EF-hand clusters -> representatives: {len(stage1_records)}")
    print(f"  Representatives unclustered in Stage 2    : {n_unclustered}")
    print(f"  EF-hand sites in final 11 clusters (C)    : {int(df['is_efhand'].sum())}")
