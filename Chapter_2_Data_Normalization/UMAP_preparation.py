# import packages
import numpy as np
import pandas as pd
import umap
import matplotlib.pyplot as plt
import seaborn as sns

# load normalized data
protein_quantification_raw_sl_tmm = pd.read_csv(
  "training_data/protein_quantification_raw_sl_tmm.csv"
)


