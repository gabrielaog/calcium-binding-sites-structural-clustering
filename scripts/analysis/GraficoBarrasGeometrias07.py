import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# ==== Dados (os seus) ====
data = {
    "Cluster": [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11],
    "Geometry": [
        "Irregular (n/a)",
        "N/A",
        "Square pyramid with one vacancy \n(equatorial) (distorted)",
        "N/A",
        "Irregular (n/a)",
        "Trigonal planar with one vacancy (regular)",
        "Irregular (n/a)",
        "Octahedron, mono-capped face with one vacancy \n(uncapped face) (distorted)",
        "Square pyramid with one vacancy \n(equatorial) (distorted)",
        "Irregular (n/a)",
        "N/A",
        "Square pyramid with one vacancy \n(equatorial) (distorted)",
        "Irregular (n/a)",
        "N/A",
        "Square planar with one vacancy (regular)",
        "Irregular (n/a)",
        "Pentagonal bipyramid with one vacancy \n(axial) (distorted)",
        "Pentagonal bipyramid with one vacancy \n(equatorial) (distorted)",
        "N/A",
        "Irregular (n/a)",
        "Square pyramid with one vacancy \n(equatorial) (distorted)",
        "Irregular (n/a)",
        "N/A",
        "Pentagonal bipyramid (distorted)",
        "Irregular (n/a)",
        "Trigonal planar with one vacancy (regular)",
        "N/A",
        "Square pyramid with one vacancy \n(equatorial) (distorted)",
        "Irregular (n/a)",
        "Trigonal prism, mono-capped square face (distorted)",
        "Irregular (n/a)",
        "Square antiprism with one vacancy (distorted)",
        "Square pyramid (distorted)"
    ],

    "Percentage": [
        42.55,
        21.28,
        6.38,
        41.67,
        27.78,
        13.89,
        58.33,
        16.67,
        8.33,
        46.43,
        14.29,
        10.71,
        26.32,
        15.79,
        10.53,
        57.89,
        10.53,
        10.53,
        40.00,
        33.33,
        13.33,
        46.67,
        26.67,
        6.67,
        37.50,
        18.75,
        12.50,
        27.27,
        18.18,
        18.18,
        40.00,
        20.00,
        10.00,
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

plt.title("Distribution of Geometries by Group (Distance ≤ 0.3)", fontsize=16, pad=15)
plt.xlabel("Group", fontsize=16)
plt.ylabel("Percentage (%)", fontsize=16)
plt.legend(title="Geometry", bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)
plt.tight_layout()
plt.show()
