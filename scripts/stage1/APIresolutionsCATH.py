import requests
import os
import re
from collections import defaultdict

"""
Cluster Representative Selection Based on Experimental Resolution

Description
-----------
This script identifies representative structures for each structural cluster
by selecting the PDB entry (or entries) with the best experimental resolution.

The input consists of previously generated cluster assignment files
(e.g., intra-superfamily clustering results). For each superfamily and cluster,
the script queries the RCSB PDB REST API to retrieve experimental resolution
data and selects all structures sharing the minimum resolution value.
-----
If multiple PDB entries share the same best resolution within a cluster,
all are retained as co-representatives.
"""

pdb_directory = "./"

def parse_cluster_file(file_path):
    """Parse the cluster file and extract PDB identifiers, chain labels, superfamily IDs, and cluster assignments"""
    clusters = {}
    current_superfamily = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip() 

            match_superfamily = re.match(r"Clusters (\d+\.\d+\.\d+\.\d+) threshold_\d+\.\d+:", line) #CATH
            print(f"match_superfamily {match_superfamily}")
            if match_superfamily:
                current_superfamily = match_superfamily.group(1)
                
                clusters[current_superfamily] = defaultdict(list)
                continue

            match_pdb = re.match(r"(\w+_\d+_[A-Za-z0-9]+): (\d+)", line)
            if match_pdb and current_superfamily:
                pdb_entry = match_pdb.group(1)
                cluster_id = int(match_pdb.group(2))

                pdb_id, pdb_number, chain = pdb_entry.rsplit("_", 2)

                clusters[current_superfamily][cluster_id].append((pdb_id,pdb_number, chain,  pdb_entry))
                

    return clusters

def get_pdb_resolution(pdb_id):
    print(f'Querying the RCSB PDB API to retrieve the resolution for {pdb_id}')
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        resolution = data.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
        return resolution
    return None

def get_best_resolutions(clusters, pdb_directory):
    """For each cluster, select all PDB entries with the best experimental resolution."""
    best_resolutions = {}
    ignored_superfamilies = set()

    for superfamily, cluster_dict in clusters.items():
        best_resolutions[superfamily] = {}

        for cluster_id, pdb_list in cluster_dict.items():
            resolutions = {}  # Dicionário para armazenar resoluções dos PDBs
            min_resolution = float("inf")

            for info_pdb in pdb_list:
                resolution = get_pdb_resolution(info_pdb[0])
                print(f'resolution: {resolution}\n')
                

                if resolution is not None:
                    resolutions[(info_pdb[0], info_pdb[1], info_pdb[2])] = resolution
                    min_resolution = min(min_resolution, resolution)
            
            best_pdbs = [
                (pdb_key[0], pdb_key[1], pdb_key[2], resolution)
                for pdb_key, resolution in resolutions.items() if resolution == min_resolution      
            ]
            
            if best_pdbs:
                best_resolutions[superfamily][cluster_id] = best_pdbs
            else:
                print(f"No PDB entry with available resolution found for cluster {cluster_id} in superfamily {superfamily}")
                all_empty = True
                for c in cluster_dict:
                    if c in best_resolutions[superfamily]:
                        all_empty = False
                        break
                if all_empty:
                    ignored_superfamilies.add(superfamily)
    
    # print("DEBUG - Best Resolutions:", best_resolutions)
    # exit()
    return best_resolutions, ignored_superfamilies 

def save_best_resolutions(best_resolutions,  ignored_superfamilies,output_file="best_resolutionsCATH.txt"):
    if ignored_superfamilies:
        print("\n----------------------------------------------- Superfamilies excluded due to absence of resolution data:")
        for sf in sorted(ignored_superfamilies):
            print(f" - {sf}")
    else:
        print("\n All superfamilies contained at least one PDB entry with available resolution.")

    with open(output_file, "w") as file:
        for superfamily, cluster_dict in best_resolutions.items():
            file.write(f"Superfamily {superfamily}:\n")
            for cluster_id, pdb_list in cluster_dict.items():
                file.write(f"  Cluster {cluster_id}:\n")
                for info_pdb in pdb_list:
                    print(f'nfo pdb {info_pdb}')
                    file.write(f"PDB {info_pdb[0]}_{info_pdb[1]}_{info_pdb[2]} - Resolution {info_pdb[3]} A\n")
            file.write("\n")

file_path = "results/clustering/CATH_StructuralClusters_TMscore_0.5.txt"

clusters= parse_cluster_file(file_path)
print(clusters)
best_resolutions, ignored_superfamilies = get_best_resolutions(clusters, pdb_directory)
save_best_resolutions(best_resolutions, ignored_superfamilies)

with open("best_resolutionsCATH.txt", "r") as file:
    print(file.read())

