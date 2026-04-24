"""
fold_composition.py

CATH and SCOPe fold breakdown for each of the 11 final clusters.

For each cluster, reports:
  - How many distinct CATH superfamilies are present and their member counts
  - How many distinct SCOPe superfamilies are present and their member counts

Input:  cluster_statistics_final_11_clusters.xlsx
Output: fold_composition.csv   — one row per (cluster, source, superfamily)
        fold_composition.txt   — human-readable summary
"""

import re
import pandas as pd
from collections import defaultdict

XLSX = "results/clustering/cluster_statistics_final_11_clusters.xlsx"


# ── reuse the same loader as efhand_check.py ────────────────────────────────

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


# ── helpers ──────────────────────────────────────────────────────────────────

def split_codes(code_str):
    """Split a semicolon-separated code string into a list of non-empty codes."""
    if code_str in ("NA", "nan", ""):
        return []
    return [c.strip() for c in code_str.split(";") if c.strip() not in ("NA", "nan", "")]


def build_composition(df, col, source_label):
    """
    For each (Cluster, superfamily_code) pair, count how many sites have
    that code.  Returns a list of dicts.
    """
    rows = []
    for cluster, group in df.groupby("Cluster", sort=False):
        counter = defaultdict(int)
        n_annotated = 0
        for _, site_row in group.iterrows():
            codes = split_codes(site_row[col])
            if codes:
                n_annotated += 1
                for code in codes:
                    counter[code] += 1
        n_total = len(group)
        n_unannotated = n_total - n_annotated
        for sf, count in sorted(counter.items(), key=lambda x: -x[1]):
            rows.append({
                "Cluster":      cluster,
                "Source":       source_label,
                "Superfamily":  sf,
                "N_sites":      count,
                "N_total":      n_total,
                "pct":          round(count / n_total * 100, 1),
            })
        if n_unannotated > 0:
            rows.append({
                "Cluster":      cluster,
                "Source":       source_label,
                "Superfamily":  "NA",
                "N_sites":      n_unannotated,
                "N_total":      n_total,
                "pct":          round(n_unannotated / n_total * 100, 1),
            })
    return rows


# ── write human-readable report ──────────────────────────────────────────────

def write_report(df, cath_rows, scope_rows, filepath):
    lines = [
        "=" * 72,
        "FOLD COMPOSITION — FINAL 11 CLUSTERS",
        "(CATH superfamily and SCOPe superfamily breakdown per cluster)",
        "=" * 72,
        "",
    ]

    clusters = df["Cluster"].unique()

    # Index rows for quick lookup
    cath_by_cluster  = defaultdict(list)
    scope_by_cluster = defaultdict(list)
    for r in cath_rows:
        cath_by_cluster[r["Cluster"]].append(r)
    for r in scope_rows:
        scope_by_cluster[r["Cluster"]].append(r)

    for cluster in clusters:
        n_total = len(df[df["Cluster"] == cluster])
        lines.append(f"{'─'*72}")
        lines.append(f"  {cluster}  (n = {n_total} sites)")
        lines.append(f"{'─'*72}")

        # CATH
        lines.append("  CATH superfamilies:")
        cath_entries = [r for r in cath_by_cluster[cluster] if r["Superfamily"] != "NA"]
        na_cath      = next((r for r in cath_by_cluster[cluster] if r["Superfamily"] == "NA"), None)
        if cath_entries:
            for r in cath_entries:
                lines.append(f"    {r['Superfamily']:<30}  {r['N_sites']:>3} sites  ({r['pct']:>5.1f}%)")
        else:
            lines.append("    (none annotated)")
        if na_cath and na_cath["N_sites"] > 0:
            lines.append(f"    {'NA':<30}  {na_cath['N_sites']:>3} sites  ({na_cath['pct']:>5.1f}%)")

        # SCOPe
        lines.append("  SCOPe superfamilies:")
        scope_entries = [r for r in scope_by_cluster[cluster] if r["Superfamily"] != "NA"]
        na_scope      = next((r for r in scope_by_cluster[cluster] if r["Superfamily"] == "NA"), None)
        if scope_entries:
            for r in scope_entries:
                lines.append(f"    {r['Superfamily']:<30}  {r['N_sites']:>3} sites  ({r['pct']:>5.1f}%)")
        else:
            lines.append("    (none annotated)")
        if na_scope and na_scope["N_sites"] > 0:
            lines.append(f"    {'NA':<30}  {na_scope['N_sites']:>3} sites  ({na_scope['pct']:>5.1f}%)")

        lines.append("")

    lines += ["=" * 72]

    with open(filepath, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"Saved: {filepath}")


# ── main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("Loading final clusters ...")
    df = load_final_clusters(XLSX)
    print(f"  {len(df)} sites in {df['Cluster'].nunique()} clusters")

    cath_rows  = build_composition(df, "CATH_superfamily",       "CATH")
    scope_rows = build_composition(df, "SCOPe_superfamily_code", "SCOPe")

    all_rows = cath_rows + scope_rows
    out_df = pd.DataFrame(all_rows)
    out_df.to_csv("results/analysis/fold_composition.csv", index=False)
    print("Saved: results/analysis/fold_composition.csv")

    write_report(df, cath_rows, scope_rows, "results/analysis/fold_composition.txt")

    # Quick console summary
    print("\n=== Distinct superfamilies per cluster ===")
    print(f"  {'Cluster':<12}  {'CATH sfs':>8}  {'SCOPe sfs':>9}  {'N sites':>7}")
    print("  " + "-" * 40)
    for cluster in df["Cluster"].unique():
        n = len(df[df["Cluster"] == cluster])
        n_cath  = len([r for r in cath_rows  if r["Cluster"] == cluster and r["Superfamily"] != "NA"])
        n_scope = len([r for r in scope_rows if r["Cluster"] == cluster and r["Superfamily"] != "NA"])
        print(f"  {cluster:<12}  {n_cath:>8}  {n_scope:>9}  {n:>7}")
