import requests
import ast
import csv
import os
import time

# ==========================================
# CONFIGURAÇÕES
# ==========================================

sites = [
"2j45_1_A", "4m0k_3_D", "3iqt_1_A", "4b5w_2_C", "4eqb_4_A",
"4pih_3_B", "5a8k_5_A", "5bqt_1_BD", "5nsf_1_A", "6uff_12_D"
]

geometry_file = "geometrias07Cluster100.txt"
scope_file = "dir.des.scope.2.08-stable.txt"
output_file = "supplementary_table11.tsv"

# ==========================================
# BAIXAR SCOPe SE NECESSÁRIO
# ==========================================

if not os.path.exists(scope_file):
    print("Baixando arquivo SCOPe...")
    os.system("wget https://scop.berkeley.edu/downloads/parse/dir.des.scope.2.08-stable.txt")

# ==========================================
# PARSE GEOMETRIA
# ==========================================

def parse_geometry(filepath):
    geometry_dict = {}

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("{"):
                continue

            data = ast.literal_eval(line)

            site = data["site"].strip()
            geometry = data.get("geometria", "NA")
            coord = data.get("coord_number", "NA")

            geometry_dict[site] = (geometry, coord)

    print("Geometrias carregadas:", len(geometry_dict))
    return geometry_dict

# ==========================================
# CARREGAR SCOPe MAPA NOME → CÓDIGO
# ==========================================

def load_scope_mapping(scope_file):
    mapping = {}

    with open(scope_file) as f:
        lines = f.readlines()[4:]

    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) < 5:
            continue

        if parts[1] == 'sf':
            code = parts[2]
            name = parts[4]
            mapping[name] = code

    print("SCOPe superfamilies carregadas:", len(mapping))
    return mapping

# ==========================================
# CONSULTAR API UMA VEZ POR PDB
# ==========================================

def get_superfamilies(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/search/pdb/select?q=pdb_id:{pdb_id}&wt=json"

    try:
        r = requests.get(url)
        data = r.json()
        docs = data['response']['docs'][0]

        cath_codes = list(set(docs.get('cath_code', [])))
        scope_names = list(set(docs.get('scop_superfamily', [])))

        return cath_codes, scope_names

    except:
        return [], []

# ==========================================
# EXECUÇÃO
# ==========================================

geometry_data = parse_geometry(geometry_file)
scope_mapping = load_scope_mapping(scope_file)

pdb_cache = {}

with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")

    writer.writerow([
        "Site",
        "PDB_ID",
        "CATH_superfamily",
        "SCOPe_superfamily_code",
        "Geometry",
        "Coordination_number"
    ])

    for site in sites:

        pdb_id = site.split("_")[0]
        short_site = "_".join(site.split("_")[:2])

        # evitar repetição API
        if pdb_id not in pdb_cache:
            cath_codes, scope_names = get_superfamilies(pdb_id)
            pdb_cache[pdb_id] = (cath_codes, scope_names)
            time.sleep(0.2)
        else:
            cath_codes, scope_names = pdb_cache[pdb_id]

        # converter SCOPe nome → código
        scope_codes = []
        for name in scope_names:
            if name in scope_mapping:
                scope_codes.append(scope_mapping[name])

        scope_codes = list(set(scope_codes))

        # pegar geometria
        geometry, coord = geometry_data.get(short_site, ("NA", "NA"))

        writer.writerow([
            site,
            pdb_id,
            ";".join(cath_codes) if cath_codes else "NA",
            ";".join(scope_codes) if scope_codes else "NA",
            geometry,
            coord
        ])

        print("Processado:", site)

print("\nTabela final gerada:", output_file)