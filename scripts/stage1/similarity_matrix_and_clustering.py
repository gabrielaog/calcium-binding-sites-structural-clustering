import os
import re
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster

"""
Structural Clustering of Calcium-Binding Minimal Functional Sites (MFSs)

Description
-----------
This script performs intra-superfamily structural clustering based on
pairwise TM-score similarity obtained from MM-align output files.

Methodology
-----------
- Similarity metric: average of the two TM-scores reported by MM-align
- Distance metric: 1 - TM-score
- Clustering algorithm: Single Linkage hierarchical clustering
- Default clustering threshold: 0.5 (distance metric)
"""

def extract_tm_scores(file_path):
    """Extract both TM-scores from an MM-align output file and return their average."""
    with open(file_path, 'r') as file:
        content = file.read()
    
    scores = re.findall(r'TM-score= ([0-9\.]+)', content)
    if len(scores) >= 2:
        media= (float(scores[0]) + float(scores[1])) / 2
        return media
    return None

def process_superfamily(superfamily_path):
    """Build a similarity matrix for a given superfamily."""
    print('-------------------------------------------------------------------------')
    print(f'Processing: {superfamily_path} \n')
    
    # Recursively collect all PDB site names
    all_files = []
    for dirpath, _, filenames in os.walk(superfamily_path):
        for file_name in filenames:
            if file_name.endswith('.pdb'):
              all_files.append(file_name.replace('.pdb', ''))

    site_names = set(all_files)
    alignments = {}

    # Read alignment files (.txt)
    for dirpath, dirnames, filenames in os.walk(superfamily_path):
        for file_name in filenames:
            if file_name.endswith('.txt'):  
                file_path = os.path.join(dirpath, file_name)
                parts = file_name.replace('.txt', '').split('_')
                site1 = '_'.join(parts[:3])
                site2 = '_'.join(parts[3:])
                
                score = extract_tm_scores(file_path)
                
                if score is not None:
                    alignments[(site1, site2)] = score
                    alignments[(site2, site1)] = score  # Ensure symmetry
    
    site_names = sorted(site_names)
    
    n= len(site_names)
    similarity_matrix = np.zeros((n, n))
    
    # Fill similarity matrix
    for i, site1 in enumerate(site_names):
        for j, site2 in enumerate(site_names):
            if i!=j and (site1, site2) in alignments:
                similarity_matrix[i, j] = alignments[(site1, site2)]
            elif i==j:
                similarity_matrix[i, j] = 1.0  # Missing alignment assumed as 0 similarity
    return site_names, similarity_matrix

def single_linkage_clustering(similarity_matrix, site_names, superfamily, all_clusters):
    """Perform Single Linkage hierarchical clustering."""
    if len(site_names) < 2: 
    # Sítio único é o cluster 1 
        all_clusters.append(f"Clusters {superfamily} threshold_0.5:") 
        all_clusters.append(f"{site_names[0]}: 1") 
        all_clusters.append("\n") 
        return
    
    distance_matrix = 1 - similarity_matrix  
    condensed_dist = squareform(distance_matrix, checks=False) 
    
    linkage_matrix = sch.linkage(condensed_dist, method='single')  

    
    threshold = 0.5
    clusters = fcluster(linkage_matrix, threshold, criterion='distance')
    
    cluster_dict = {site_names[i]: clusters[i] for i in range(len(site_names))}

    all_clusters.append(f"Clusters {superfamily} threshold_{threshold}:")
    for site, cluster in cluster_dict.items():
        all_clusters.append(f"{site}: {cluster}")
    all_clusters.append("\n")  

def main(root_dir, output_dir, db_to_run):
    """Process all superfamilies from a selected database (CATH or SCOPe)."""
    results = {}
    all_clusters = []
    
    db_path = os.path.join(root_dir, db_to_run)
    output_file = os.path.join(output_dir, f"{db_to_run}_all_clusters05.txt")
    
    if os.path.exists(db_path):
        for superfamily in os.listdir(db_path):
            superfamily_path = os.path.join(db_path, superfamily)
            if os.path.isdir(superfamily_path):
                sites, matrix = process_superfamily(superfamily_path)
                if sites is not None:
                    results[superfamily] = pd.DataFrame(matrix, index=sites, columns=sites)
                    single_linkage_clustering(matrix, sites, superfamily, all_clusters)
    
    output_filename = os.path.join(output_dir, f"similarity_matrices_{db_to_run}.xlsx")
    
    with pd.ExcelWriter(output_filename) as writer:
        for superfamily, df in results.items():
            df.to_excel(writer, sheet_name=superfamily[:31]) 
    print(f'Similarity matrices saved to {output_filename}')

    with open(output_file, "w") as f:
        f.write("\n".join(all_clusters))

    print(f"Cluster assignments saved to {output_file}")


if __name__ == "__main__":
    root_directory = "./"  
    output_directory = "./"  
    db_to_run = 'SCOPe'  # CATH or 'SCOPe'
    
    
    main(root_directory, output_directory, db_to_run)
    

