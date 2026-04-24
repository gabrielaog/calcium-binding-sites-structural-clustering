import pandas as pd
import matplotlib.pyplot as plt

# ==== Dados ====
data = {
    "Cluster": [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,16,16,17,17,17,18,18,18,19,19,19,20,20,20,21,21,21,22,22,22,23,23,23,24,24,24,25,25,25],
    "Geometria": [
        "irregular (n/a)", "pentagonal bipyramid (distorted)", "n/a",
        "square antiprism (distorted)", "n/a", "square antiprism with a vacancy (distorted)",
        "irregular (n/a)", "pentagonal bipyramid with a vacancy (axial) (distorted)", "n/a",
        "irregular (n/a)", "n/a", "pentagonal bipyramid (distorted)",
        "irregular (n/a)", "pentagonal bipyramid (distorted)", "octahedron (regular)",
        "irregular (n/a)", "pentagonal bipyramid (distorted)", "pentagonal bipyramid with a vacancy (axial) (distorted)",
        "irregular (n/a)", "square pyramid with a vacancy (equatorial) (distorted)", "square pyramid (regular)",
        "irregular (n/a)", "pentagonal bipyramid (distorted)", "octahedron (regular)",
        "irregular (n/a)", "octahedron (regular)", "pentagonal bipyramid (distorted)",
        "irregular (n/a)", "octahedron, face monocapped with a vacancy (non-capped face) (distorted)", "square pyramid with a vacancy (equatorial) (distorted)",
        "irregular (n/a)", "n/a", "square pyramid with a vacancy (equatorial) (distorted)",
        "irregular (n/a)", "pentagonal bipyramid (distorted)", "square pyramid with a vacancy (equatorial) (distorted)",
        "irregular (n/a)", "pentagonal bipyramid (distorted)", "n/a",
        "irregular (n/a)", "pentagonal bipyramid with a vacancy (equatorial) (distorted)", "n/a",
        "irregular (n/a)", "n/a", "octahedron (distorted)",
        "irregular (n/a)", "pentagonal bipyramid with a vacancy (equatorial) (distorted)", "pentagonal bipyramid with a vacancy (equatorial) (regular)",
        "n/a", "irregular (n/a)", "pentagonal bipyramid (distorted)",
        "irregular (n/a)", "n/a", "square pyramid with a vacancy (equatorial) (distorted)",
        "irregular (n/a)", "square pyramid (regular)", "n/a",
        "irregular (n/a)", "octahedron (regular)", "trigonal plane with a vacancy (regular)",
        "irregular (n/a)", "square pyramid with a vacancy (equatorial) (distorted)", "pentagonal bipyramid with a vacancy (equatorial) (distorted)",
        "irregular (n/a)", "square pyramid with a vacancy (equatorial) (distorted)", "n/a",
        "irregular (n/a)", "square pyramid with a vacancy (equatorial) (distorted)", "pentagonal bipyramid (distorted)",
        "trigonal plane with a vacancy (regular)", "irregular (n/a)", "pentagonal bipyramid (distorted)",
        "irregular (n/a)", "square pyramid with a vacancy (equatorial) (distorted)", "trigonal prism, square-face monocapped (distorted)"
    ],
    "Porcentagem": [
        33.28, 10.87, 7.36, 21.43, 14.29, 14.29, 42.02, 13.83, 6.38, 31.86, 10, 6.78, 27.5, 15, 10, 26.47, 11.76, 8.82,
        32.14, 14.29, 7.14, 33.33, 9.09, 9.09, 33.33, 7.62, 7.62, 59.09, 9.09, 4.55, 33, 9.36, 6.4, 36.36, 13.64, 9.09,
        38.1, 14.29, 4.76, 50, 30, 10, 35.76, 8.61, 7.28, 31.58, 15.79, 10.53, 27.78, 22.22, 11.11, 43.33, 26.67, 6.67,
        33.33, 16.67, 11.11, 53.85, 15.38, 7.69, 27.78, 11.11, 8.33, 38.46, 23.08, 7.69, 31.82, 13.64, 9.09, 20, 20, 20,
        23.08, 11.54, 11.54
    ]
}

df = pd.DataFrame(data)

# ==== Pivot para formato largo ====
pivot_df = df.pivot_table(index="Cluster", columns="Geometria", values="Porcentagem", fill_value=0)

# ==== Gráfico de barras empilhadas ====
pivot_df.plot(kind='bar', stacked=True, figsize=(14, 8), width=0.8)

plt.title("Distribuição das Geometrias por Cluster (%)", fontsize=14, pad=15)
plt.xlabel("Cluster", fontsize=12)
plt.ylabel("Porcentagem (%)", fontsize=12)
plt.legend(title="Geometria", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
plt.tight_layout()
plt.show()
