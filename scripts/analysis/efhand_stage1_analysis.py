"""
efhand_stage1_analysis.py

Identifica sítios EF-hand nos arquivos de clustering da Stage 1
(CATH_StructuralClusters_TMscore_0.5.txt e SCOPe_StructuralClusters_TMscore_0.5.txt)
e salva os resultados em efhand_stage1_results.txt.

Critérios EF-hand:
  CATH  → superfamília começa com "1.10.238"
  SCOPe → superfamília começa com "a.7."
"""

import os
from collections import defaultdict

# --- configuração -----------------------------------------------------------

CATH_FILE  = "results/clustering/CATH_StructuralClusters_TMscore_0.5.txt"
SCOPE_FILE = "results/clustering/SCOPe_StructuralClusters_TMscore_0.5.txt"
OUTPUT     = "efhand_stage1_results.txt"

CATH_EF_PREFIX  = "1.10.238"
SCOPE_EF_PREFIX = "a.7."

# ----------------------------------------------------------------------------

def parse_clusters(filepath):
    """
    Lê um arquivo de clustering da Stage 1 e retorna um dicionário:
        { superfamilia: [ (site, cluster_id), ... ] }
    """
    result = {}
    current_sf = None

    with open(filepath, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith("Clusters "):
                # ex: "Clusters 1.10.238.10 threshold_0.5:"
                parts = line.split()
                current_sf = parts[1]
                result[current_sf] = []
            elif current_sf and ":" in line:
                site, cluster_str = line.split(":", 1)
                result[current_sf].append((site.strip(), int(cluster_str.strip())))

    return result


def is_efhand(superfamily, source):
    if source == "CATH":
        return superfamily.startswith(CATH_EF_PREFIX)
    if source == "SCOPe":
        return superfamily.startswith(SCOPE_EF_PREFIX)
    return False


def summarise(sf_dict, source):
    """
    Filtra apenas as superfamílias EF-hand e devolve lista de dicts com:
        superfamily, source, n_sites, n_clusters,
        clusters: { cluster_id: [site, ...] }
    """
    records = []
    for sf, sites in sf_dict.items():
        if not is_efhand(sf, source):
            continue
        by_cluster = defaultdict(list)
        for site, cid in sites:
            by_cluster[cid].append(site)
        records.append({
            "superfamily": sf,
            "source": source,
            "n_sites": len(sites),
            "n_clusters": len(by_cluster),
            "clusters": dict(by_cluster),
        })
    records.sort(key=lambda r: r["superfamily"])
    return records


def write_results(records_cath, records_scope,
                  total_cath, total_scope, output_path):

    total_ef_cath  = sum(r["n_sites"] for r in records_cath)
    total_ef_scope = sum(r["n_sites"] for r in records_scope)
    total_all      = total_cath + total_scope
    total_ef_all   = total_ef_cath + total_ef_scope

    lines = []
    lines.append("=" * 70)
    lines.append("ANÁLISE EF-HAND — STAGE 1 (clustering intra-superfamília)")
    lines.append("=" * 70)
    lines.append("")

    # --- resumo geral -------------------------------------------------------
    lines.append("RESUMO GERAL")
    lines.append("-" * 40)
    lines.append(f"  CATH : {total_cath:>6} sítios totais  |  {total_ef_cath:>4} EF-hand "
                 f"({100*total_ef_cath/total_cath:.2f}%)")
    lines.append(f"  SCOPe: {total_scope:>6} sítios totais  |  {total_ef_scope:>4} EF-hand "
                 f"({100*total_ef_scope/total_scope:.2f}%)")
    lines.append(f"  Total: {total_all:>6} sítios totais  |  {total_ef_all:>4} EF-hand "
                 f"({100*total_ef_all/total_all:.2f}%)")
    lines.append("")

    # --- CATH ----------------------------------------------------------------
    lines.append("SUPERFAMÍLIAS EF-HAND — CATH  (prefixo '1.10.238')")
    lines.append("-" * 70)
    for r in records_cath:
        lines.append(f"\n  Superfamília : {r['superfamily']}")
        lines.append(f"  Total sítios : {r['n_sites']}")
        lines.append(f"  Nº clusters  : {r['n_clusters']}")
        for cid in sorted(r["clusters"]):
            members = r["clusters"][cid]
            lines.append(f"    Cluster {cid} ({len(members)} sítios): "
                         + ", ".join(members))
    lines.append("")

    # --- SCOPe ---------------------------------------------------------------
    lines.append("SUPERFAMÍLIAS EF-HAND — SCOPe  (prefixo 'a.7.')")
    lines.append("-" * 70)
    for r in records_scope:
        lines.append(f"\n  Superfamília : {r['superfamily']}")
        lines.append(f"  Total sítios : {r['n_sites']}")
        lines.append(f"  Nº clusters  : {r['n_clusters']}")
        for cid in sorted(r["clusters"]):
            members = r["clusters"][cid]
            lines.append(f"    Cluster {cid} ({len(members)} sítios): "
                         + ", ".join(members))
    lines.append("")

    # --- lista plana de todos os sítios EF-hand ------------------------------
    lines.append("LISTA COMPLETA DE SÍTIOS EF-HAND (Stage 1)")
    lines.append("-" * 70)
    lines.append(f"{'Site':<25} {'Source':<6} {'Superfamily':<20} {'Cluster':>7}")
    lines.append("-" * 70)
    for r in records_cath + records_scope:
        for cid in sorted(r["clusters"]):
            for site in r["clusters"][cid]:
                lines.append(f"{site:<25} {r['source']:<6} "
                              f"{r['superfamily']:<20} {cid:>7}")
    lines.append("")
    lines.append("=" * 70)

    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"Resultado salvo em: {output_path}")
    return total_ef_all, total_all


# ----------------------------------------------------------------------------

if __name__ == "__main__":
    script_dir = os.path.dirname(os.path.abspath(__file__))

    cath_path  = os.path.join(script_dir, CATH_FILE)
    scope_path = os.path.join(script_dir, SCOPE_FILE)
    out_path   = os.path.join(script_dir, OUTPUT)

    print("Lendo CATH ...")
    cath_dict  = parse_clusters(cath_path)
    print("Lendo SCOPe ...")
    scope_dict = parse_clusters(scope_path)

    total_cath  = sum(len(v) for v in cath_dict.values())
    total_scope = sum(len(v) for v in scope_dict.values())

    records_cath  = summarise(cath_dict,  "CATH")
    records_scope = summarise(scope_dict, "SCOPe")

    n_ef, n_total = write_results(
        records_cath, records_scope,
        total_cath, total_scope, out_path
    )

    print(f"EF-hand encontrados: {n_ef} / {n_total} sítios totais "
          f"({100*n_ef/n_total:.2f}%)")
