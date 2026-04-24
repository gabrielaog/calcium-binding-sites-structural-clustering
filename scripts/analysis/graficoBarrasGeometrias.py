import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# ==== Dados (os seus) ====
data = {
    "Cluster": [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,16,16,17,17,17,18,18,18,19,19,19,20,20,20,21,21,21,22,22,22,23,23,23,24,24,24,25,25,25],
    "Geometry": [
        "Irregular (n/a)", "Pentagonal bipyramid (distorted)", "N/A",
        "Square antiprism (distorted)", "N/A", "Square antiprism with vacancy \n(distorted)",
        "Irregular (n/a)", "Pentagonal bipyramid with vacancy \n(axial) (distorted)", "N/A",
        "Irregular (n/a)", "N/A", "Pentagonal bipyramid (distorted)",
        "Irregular (n/a)", "Pentagonal bipyramid (distorted)", "Octahedron (regular)",
        "Irregular (n/a)", "Pentagonal bipyramid (distorted)", "Pentagonal bipyramid with vacancy \n(axial) (distorted)",
        "Irregular (n/a)", "Square pyramid with vacancy \n(equatorial) (distorted)", "Square pyramid (regular)",
        "Irregular (n/a)", "Pentagonal bipyramid (distorted)", "Octahedron (regular)",
        "Irregular (n/a)", "Octahedron (regular)", "Pentagonal bipyramid (distorted)",
        "Irregular (n/a)", "Octahedron, mono-capped face with \nvacancy (uncapped face) (distorted)", "Square pyramid with vacancy \n(equatorial) (distorted)",
        "Irregular (n/a)", "N/A", "Square pyramid with vacancy \n(equatorial) (distorted)",
        "Irregular (n/a)", "Pentagonal bipyramid (distorted)", "Square pyramid with vacancy \n(equatorial) (distorted)",
        "Irregular (n/a)", "Pentagonal bipyramid (distorted)", "N/A",
        "Irregular (n/a)", "Pentagonal bipyramid with vacancy \n(equatorial) (distorted)", "N/A",
        "Irregular (n/a)", "N/A", "Octahedron (distorted)",
        "Irregular (n/a)", "Pentagonal bipyramid with vacancy \n(equatorial) (distorted)", "Pentagonal bipyramid with vacancy \n(equatorial) (regular)",
        "N/A", "Irregular (n/a)", "Pentagonal bipyramid (distorted)",
        "Irregular (n/a)", "N/A", "Square pyramid with vacancy \n(equatorial) (distorted)",
        "Irregular (n/a)", "Square pyramid (regular)", "N/A",
        "Irregular (n/a)", "Octahedron (regular)", "Trigonal planar with vacancy (regular)",
        "Irregular (n/a)", "Square pyramid with vacancy \n(equatorial) (distorted)", "Pentagonal bipyramid with vacancy \n(equatorial) (distorted)",
        "Irregular (n/a)", "Square pyramid with vacancy \n(equatorial) (distorted)", "N/A",
        "Irregular (n/a)", "Square pyramid with vacancy \n(equatorial) (distorted)", "Pentagonal bipyramid (distorted)",
        "Trigonal planar with vacancy (regular)", "Irregular (n/a)", "Pentagonal bipyramid (distorted)",
        "Irregular (n/a)", "Square pyramid with vacancy \n(equatorial) (distorted)", "Trigonal prism with square face \nmono-capped (distorted)"

    ],
    "Percentage": [
        33.28,10.87,7.36,21.43,14.29,14.29,42.02,13.83,6.38,31.86,10,6.78,27.5,15,10,26.47,11.76,8.82,
        32.14,14.29,7.14,33.33,9.09,9.09,33.33,7.62,7.62,59.09,9.09,4.55,33,9.36,6.4,36.36,13.64,9.09,
        38.1,14.29,4.76,50,30,10,35.76,8.61,7.28,31.58,15.79,10.53,27.78,22.22,11.11,43.33,26.67,6.67,
        33.33,16.67,11.11,53.85,15.38,7.69,27.78,11.11,8.33,38.46,23.08,7.69,31.82,13.64,9.09,20,20,20,
        23.08,11.54,11.54
    ]
}

df = pd.DataFrame(data)

# ==== Pivot para formato largo ====
pivot_df = df.pivot_table(index="Cluster", columns="Geometry", values="Percentage", fill_value=0)

# (Maior porcentagem média embaixo)
mean_order = pivot_df.mean().sort_values(ascending=False).index
pivot_df = pivot_df[mean_order]

# ==== Paleta sem repetição ====
colors = sns.color_palette("tab20", len(pivot_df.columns))  # paleta categórica sem repetição

# ==== Gráfico de barras empilhadas ====
ax = pivot_df.plot(kind='bar', width= 0.9, stacked=True, color=colors, figsize=(14,8))
for c in ax.containers:
    ax.bar_label(c, fmt='%.0f', label_type='center', fontsize=11)

plt.title("Distribution of Geometries by Group (Distance ≤ 0.5)", fontsize=16, pad=15)
plt.xlabel("Group", fontsize=16)
plt.ylabel("Percentage (%)", fontsize=16)
plt.legend(title="Geometry", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=11)
plt.tight_layout()
plt.show()
