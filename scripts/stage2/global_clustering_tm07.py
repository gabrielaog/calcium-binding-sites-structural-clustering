import pandas as pd
from collections import defaultdict

"""
Global Clustering of Representative Calcium-Binding Sites

Description
-----------
This script performs the second-stage global clustering of representative
Minimal Functional Sites (MFSs) using a distance matrix derived from
average TM-scores.

The input is a pairwise distance matrix (DistanceMatrix.csv), where the
distance between two sites is defined as:

    d = 1 - average(TM-score)

The script iteratively scans the matrix to identify the nearest neighbor
for each site and assigns sites to clusters when their minimum pairwise
distance satisfies a strict structural similarity cutoff.

Clustering Strategy
-------------------
The clustering procedure follows a nearest-neighbor strategy similar to
single-linkage clustering, but with an additional strict distance constraint.
Two sites (or clusters) are merged only when their minimum pairwise distance
is strictly below the predefined cutoff:

    d < 0.3   (equivalent to TM-score ≥ 0.7)

This rule prevents the chaining effect typically observed in classical
single-linkage clustering and ensures that only highly similar structural
environments are grouped together.

Filtering of Small Clusters
---------------------------
After cluster construction, clusters containing fewer than ten members
are discarded. This filtering step focuses the analysis on structurally
recurrent calcium-binding site architectures while removing rare or
isolated structural configurations.

Usage
-----
The script is intended to be executed after the selection of representative
sites from intra-superfamily clustering, using the distance matrix produced
from their structural alignments.
"""
count=0
# Read the distance matrix (first column used as index)
df = pd.read_csv("DistanceMatrix.csv", index_col=0)
groups = {}
group_id = 1


def find_group(key):
    return groups.get(key, None)

for row_name, row in df.iterrows():   
    # Filter out zero values (self-distance)
    filtered_row = row[row != 0.0]

    if filtered_row.empty:
        continue

    # Identify the column with the minimum distance
    column_name = filtered_row.idxmin()

    # Check whether the row or column already belongs to a cluster
    row_group = find_group(row_name)
    column_group = find_group(column_name)

    if row_group and column_group:
        continue
    elif  (row[column_name]==1.0):
        continue
    elif  (row[column_name]>=0.3):
        count+=1
        continue
    elif row_group:
        groups[column_name] = {'group': row_group['group'], 'closest_pdb': row_name, 'distance_value': row[column_name]}
    elif column_group:
        groups[row_name] = {'group': column_group['group'], 'closest_pdb': column_name, 'distance_value': row[column_name]}
    else:
        # Neither belongs to a cluster — create a new cluster
        groups[row_name] = {'group': group_id, 'closest_pdb': column_name, 'distance_value': row[column_name]}
        groups[column_name] = {'group': group_id, 'closest_pdb': row_name, 'distance_value': row[column_name]}
        group_id += 1

members_by_group = defaultdict(list)
for key, valor in groups.items():
    members_by_group[valor['group']].append(key)

# Identify clusters with fewer than 10 members
small_groups = {gid for gid, membros in members_by_group.items() if len(membros) < 10}

filtered_groups = {key: valor for key, valor in groups.items() if valor['group'] not in small_groups}

with open('resultadoMaiorQue07.txt', 'w') as f:
    for key, valor in filtered_groups.items():
        f.write(f"'{key}': {valor},\n")



