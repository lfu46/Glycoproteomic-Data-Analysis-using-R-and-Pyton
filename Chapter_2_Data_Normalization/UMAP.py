# import packages
import numpy as np
import pandas as pd
import umap
import matplotlib.pyplot as plt
import seaborn as sns

# load formatted normalized data
protein_quantification_raw_sl_tmm_data = pd.read_csv(
  "training_data/protein_quantification_raw_sl_tmm_data.csv"
)

# traing reducer using data
reducer = umap.UMAP()

normlized_data = protein_quantification_raw_sl_tmm_data.iloc[:, 2:].values

embedding = reducer.fit_transform(normlized_data)
normlized_data.shape
embedding.shape

# Map categorical values to numbers
color_mapping = protein_quantification_raw_sl_tmm_data.Treatment.map({
    "Iron": 0,
    "Copper": 1,
    "Ctrl": 2
})

# Generate colors from Seaborn palette
colors = [sns.color_palette()[x] for x in color_mapping]

# UMAP projection
plt.scatter(
    embedding[:, 0],
    embedding[:, 1],
    c=colors)
plt.gca().set_aspect('equal', 'datalim')
plt.title('UMAP projection of the normalized dataset', fontsize=16)
handles = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=sns.color_palette()[0], markersize=10, label="Iron"),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=sns.color_palette()[1], markersize=10, label="Copper"),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=sns.color_palette()[2], markersize=10, label="Ctrl")
]
plt.legend(handles=handles, title="Treatment")

plt.savefig("training_data/umap_plot.png", dpi=300)
