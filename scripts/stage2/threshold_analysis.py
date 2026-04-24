import pandas as pd
import numpy as np
from collections import defaultdict
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

"""
Threshold Sweep Analysis — Global Clustering of Ca2+-Binding Sites

Runs the same clustering algorithm as global_clustering_tm07.py for
multiple distance thresholds and computes quality metrics:
  - Number of valid clusters (>= 10 members)
  - Total members clustered
  - Mean and max cluster size
  - % intra-cluster pairwise distances < 0.5
  - Silhouette score

Input:  DistanceMatrix.csv  (same file used in global_clustering_tm07.py)
Output: threshold_analysis.csv   — metrics table
        threshold_analysis.png   — 4-panel figure
"""

# ── Clustering (exact same logic as global_clustering_tm07.py) ────────────────

def run_clustering(df, threshold):
    """Runs the nearest-neighbour greedy clustering for a given threshold."""
    groups = {}
    group_id = 1

    for row_name, row in df.iterrows():
        filtered_row = row[row != 0.0]
        if filtered_row.empty:
            continue

        column_name = filtered_row.idxmin()
        row_group    = groups.get(row_name)
        column_group = groups.get(column_name)

        if row_group and column_group:
            continue
        if row[column_name] == 1.0 or row[column_name] >= threshold:
            continue

        if row_group:
            groups[column_name] = {
                'group': row_group['group'],
                'closest_pdb': row_name,
                'distance_value': row[column_name]
            }
        elif column_group:
            groups[row_name] = {
                'group': column_group['group'],
                'closest_pdb': column_name,
                'distance_value': row[column_name]
            }
        else:
            groups[row_name]    = {'group': group_id, 'closest_pdb': column_name, 'distance_value': row[column_name]}
            groups[column_name] = {'group': group_id, 'closest_pdb': row_name,    'distance_value': row[column_name]}
            group_id += 1

    # Group members
    members_by_group = defaultdict(list)
    for key, val in groups.items():
        members_by_group[val['group']].append(key)

    # Keep only clusters with >= 10 members
    valid = {gid: mbs for gid, mbs in members_by_group.items() if len(mbs) >= 10}
    return valid


# ── Metrics ───────────────────────────────────────────────────────────────────

def compute_metrics(df, valid_clusters, threshold):
    if not valid_clusters:
        return None

    sizes       = [len(m) for m in valid_clusters.values()]
    all_members = [m for mbs in valid_clusters.values() for m in mbs]

    # % intra-cluster pairwise distances < 0.5
    intra = []
    for mbs in valid_clusters.values():
        mbs_valid = [m for m in mbs if m in df.index]
        for i in range(len(mbs_valid)):
            for j in range(i + 1, len(mbs_valid)):
                intra.append(df.loc[mbs_valid[i], mbs_valid[j]])
    pct_below_05 = np.mean([d < 0.5 for d in intra]) * 100 if intra else 0.0

    # Silhouette score (requires >= 2 clusters and >= 4 members total)
    silhouette = None
    valid_sites, labels = [], []
    for gid, mbs in valid_clusters.items():
        for m in mbs:
            if m in df.index:
                valid_sites.append(m)
                labels.append(gid)

    if len(set(labels)) >= 2 and len(valid_sites) >= 4:
        dist_sub = df.loc[valid_sites, valid_sites].values.copy()
        np.fill_diagonal(dist_sub, 0)
        dist_sub = np.clip(dist_sub, 0, None)
        try:
            silhouette = silhouette_score(dist_sub, labels, metric='precomputed')
        except Exception as e:
            print(f"  [silhouette error @ threshold={threshold}] {e}")

    return {
        'threshold':     threshold,
        'n_clusters':    len(valid_clusters),
        'total_members': len(all_members),
        'mean_size':     round(float(np.mean(sizes)), 2),
        'max_size':      int(max(sizes)),
        'pct_below_05':  round(pct_below_05, 2),
        'silhouette':    round(silhouette, 4) if silhouette is not None else None,
    }


# ── Main ──────────────────────────────────────────────────────────────────────

df = pd.read_csv("DistanceMatrix.csv", index_col=0)
print(f"Distance matrix loaded: {df.shape[0]} x {df.shape[1]}\n")

thresholds = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60]
records = []

for t in thresholds:
    valid   = run_clustering(df, t)
    metrics = compute_metrics(df, valid, t)
    if metrics:
        records.append(metrics)
        sil_str = f"{metrics['silhouette']:.4f}" if metrics['silhouette'] is not None else "N/A"
        print(
            f"threshold={t:.2f} | clusters={metrics['n_clusters']:2d} | "
            f"members={metrics['total_members']:4d} | mean_size={metrics['mean_size']:5.1f} | "
            f"pct<0.5={metrics['pct_below_05']:5.1f}% | silhouette={sil_str}"
        )
    else:
        print(f"threshold={t:.2f} | no valid clusters (all < 10 members)")

results_df = pd.DataFrame(records)
results_df.to_csv("results/analysis/threshold_analysis.csv", index=False)
print(f"\nMetrics saved → results/analysis/threshold_analysis.csv")

# ── Figure ────────────────────────────────────────────────────────────────────

fig, axes = plt.subplots(2, 2, figsize=(11, 8))
fig.suptitle(
    "Effect of Distance Threshold on Global Clustering of Ca²⁺-Binding Sites",
    fontsize=13, fontweight='bold'
)

x         = results_df['threshold']
highlight = 0.30

def mark(ax):
    ax.axvline(highlight, color='red', linestyle='--', linewidth=1.2, label='d = 0.3 (chosen)')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.grid(True, alpha=0.35)
    ax.legend(fontsize=8)

# (A) Number of clusters
ax = axes[0, 0]
ax.plot(x, results_df['n_clusters'], 'bo-', linewidth=1.8, markersize=7)
mark(ax)
ax.set_xlabel("Distance threshold  (d = 1 − TM-score)")
ax.set_ylabel("Number of clusters  (≥ 10 members)")
ax.set_title("(A) Number of clusters")

# (B) Mean cluster size
ax = axes[0, 1]
ax.plot(x, results_df['mean_size'], 'ro-', linewidth=1.8, markersize=7)
mark(ax)
ax.set_xlabel("Distance threshold  (d = 1 − TM-score)")
ax.set_ylabel("Mean cluster size")
ax.set_title("(B) Mean cluster size")

# (C) Intra-cluster coherence
ax = axes[1, 0]
ax.plot(x, results_df['pct_below_05'], 'go-', linewidth=1.8, markersize=7)
mark(ax)
ax.set_xlabel("Distance threshold  (d = 1 − TM-score)")
ax.set_ylabel("% intra-cluster distances < 0.5")
ax.set_title("(C) Intra-cluster structural coherence")

# (D) Silhouette score
ax = axes[1, 1]
sil_df = results_df.dropna(subset=['silhouette'])
ax.plot(sil_df['threshold'], sil_df['silhouette'], 'mo-', linewidth=1.8, markersize=7)
mark(ax)
ax.set_xlabel("Distance threshold  (d = 1 − TM-score)")
ax.set_ylabel("Silhouette score")
ax.set_title("(D) Silhouette score")

plt.tight_layout()
plt.savefig("results/analysis/threshold_analysis.png", dpi=300, bbox_inches='tight')
plt.show()
print("Figure saved → results/analysis/threshold_analysis.png")
