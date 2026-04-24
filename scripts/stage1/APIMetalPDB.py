import requests
import re

def removeNaoCalcio(dados, chave):
    """Remove todos os registros cujo valor da chave não seja 'Calcium'."""
    return [d for d in dados if d.get(chave) != "Calcium"]


def remove_duplicados(dados, chave):
    """Remove duplicados de uma lista de dicionários com base no valor de uma chave."""
    vistos, unicos = set(), []
    for d in dados:
        valor = d.get(chave)
        if valor not in vistos:
            unicos.append(d)
            vistos.add(valor)
    return unicos


def consulta_metalpdb(pdb_id):
    """
    Consulta a API do MetalPDB para obter informações sobre geometrias de metais em um PDB.
    Retorna lista de dicionários com metal, cadeia, geometria e número de coordenação.
    """
    url = f"http://metalpdb.cerm.unifi.it/api?query=pdb:{pdb_id}&columns=metal,chain,geometry,coordination"
    resp = requests.get(url)

    if resp.status_code != 200:
        print(f"Erro ao acessar {pdb_id}: {resp.status_code}")
        return []

    data = resp.json()
    if not data:
        print(f"Nenhum sítio metálico encontrado em {pdb_id}")
        return []

    resultados = []
    for site in data:
        site_id = site.get("site")
        if site_id in sitios:
            for metal in site.get("metals", []):
                if metal["name"] != "Calcium":
                  pdbNaoCalcio.append({
                      "pdb_id": pdb_id,
                      "site": site_id,
                      "metal": metal["symbol"],
                      "nome": metal["name"],
                      "geometria": metal["geometry"],
                      "coord_number": metal["coordination"]
                  })
                else:
                  resultados.append({
                      "pdb_id": pdb_id,
                      "site": site_id,
                      "metal": metal["symbol"],
                      "nome": metal["name"],
                      "geometria": metal["geometry"],
                      "coord_number": metal["coordination"]
                  })
    return resultados


# 🔹 Exemplo de uso
pdbs = [
    '2j45_1_A', '4m0k_3_D', '3iqt_1_A', '4b5w_2_C', '4eqb_4_A', '4pih_3_B', '5a8k_5_A', '5bqt_1_BD', '5nsf_1_A', '6uff_12_D'
]
pdbNaoCalcio = []
# Regex: captura até o segundo "_"
pattern = r'^([^_]+_[^_]+)'
sitios = [re.match(pattern, pdb).group(1) for pdb in pdbs if re.match(pattern, pdb)]

# 🔹 Processamento e salvamento
for pdb in pdbs:
    pdb_id = pdb[:4]

    dados = consulta_metalpdb(pdb_id)
    dados = remove_duplicados(dados, "site")

    with open("geometrias07Cluster100.txt", "a", encoding="utf-8") as f:
        f.writelines(str(d) + "\n" for d in dados)

with open("naoCalcio.txt", "w", encoding="utf-8") as f:
    f.writelines(str(d) + "\n" for d in pdbNaoCalcio)
