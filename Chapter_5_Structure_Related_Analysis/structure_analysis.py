# import packages
import pandas as pd
import numpy as np
import tqdm as tqdm
import warnings
import os

# Suppress FutureWarnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# localCIDER-related packages
from localcider.sequenceParameters import SequenceParameters

# StructureMap-related packages
import structuremap.utils
structuremap.utils.set_logger()
from structuremap.processing import download_alphafold_cif, download_alphafold_pae, format_alphafold_data, annotate_accessibility, get_smooth_score, annotate_proteins_with_idr_pattern, get_extended_flexible_pattern, get_proximity_pvals, perform_enrichment_analysis, perform_enrichment_analysis_per_protein, evaluate_ptm_colocalization, extract_motifs_in_proteome

### Configuration: Set your data directories here
# Users should update these paths to point to their own data directories
# For the tutorial, place your AlphaFold CIF and PAE files in the training_data folder
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CIF_DIRECTORY = os.path.join(SCRIPT_DIR, "training_data", "alphafold_cif")
PAE_DIRECTORY = os.path.join(SCRIPT_DIR, "training_data", "alphafold_pae")

### Peptide Sequence Physicochemical Property Analysis
# import your sequences here
# here I used Nterm project result
HEK_Nterm_Kd_half_life_sequence = pd.read_csv(
  os.path.join(SCRIPT_DIR, 'training_data/HEK_Nterm_Kd_half_life_sequence.csv'),
  header = 0,
  index_col = None
)

# get sequence parameters for Nterm_13mer
HEK_Nterm_Kd_half_life_sequence.loc[:, "Nterm_13mers_sequence_parameter"]= HEK_Nterm_Kd_half_life_sequence["Nterm_13mer"].astype(pd.StringDtype()).apply(SequenceParameters)

# hydropathy
HEK_Nterm_Kd_half_life_sequence.loc[:, "hydropathy"] = HEK_Nterm_Kd_half_life_sequence["Nterm_13mers_sequence_parameter"].apply(
  lambda obj: obj.get_mean_hydropathy()
)

# isoelectric point
HEK_Nterm_Kd_half_life_sequence.loc[:, "isoelectric_point"] = HEK_Nterm_Kd_half_life_sequence["Nterm_13mers_sequence_parameter"].apply(
  lambda obj: obj.get_isoelectric_point()
)

## calculate hydropathy, FCR, NCPR, isoelectric point and kappa values
# methods = {
#     "hydropathy": "get_mean_hydropathy",
#     "FCR": "get_FCR",
#     "NCPR": "get_NCPR",
#     "isoelectric_point": "get_isoelectric_point",
#     "kappa": "get_kappa"
# }
# 
# for col, method in tqdm(methods.items()):
#     HEK_Nterm_Kd_half_life_sequence.loc[:, col] = HEK_Nterm_Kd_half_life_sequence["Nterm_13mers_sequence_parameter"].apply(
#         lambda obj: getattr(obj, method)() if obj is not None else None
#     )

### Protein Modification Structure Analysis
# import your modification
# here I used Nterm project result
common_Nterm_protein = pd.read_csv(
  os.path.join(SCRIPT_DIR, 'training_data/common_Nterm_protein.csv'),
  header = 0,
  index_col = None
)

common_Nterm_protein_list = common_Nterm_protein["UniProt_Accession"].unique().tolist()

## download AlphaFold data
## I already have the files here so I didn't run the code.
# crystallographic information file
# valid_proteins_cif, invalid_proteins_cif, existing_proteins_cif = download_alphafold_cif(
#     proteins = total_Nterm_protein_list,
#     out_folder = "data_source/Nterm_structuremap/Nterm_degradation_cif"
# )

# predicted aligned error
# valid_proteins_pae, invalid_proteins_pae, existing_proteins_pae = download_alphafold_pae(
#     proteins = total_Nterm_protein_list,
#     out_folder = "data_source/Nterm_structuremap/Nterm_degradation_pae"
# )

# format AlphaFold data input
# Update CIF_DIRECTORY at the top of this script to point to your AlphaFold CIF files
common_Nterm_alphafold_annotation = format_alphafold_data(
  directory = CIF_DIRECTORY,
  protein_ids = common_Nterm_protein_list
)

# annotate prediction-aware part-sphere exposure (pPSE) values
# Update PAE_DIRECTORY at the top of this script to point to your AlphaFold PAE files
common_Nterm_full_sphere_exposure = annotate_accessibility(
    df = common_Nterm_alphafold_annotation,
    max_dist = 24,
    max_angle = 180,
    error_dir = PAE_DIRECTORY
)

common_Nterm_alphafold_accessibility = common_Nterm_alphafold_annotation.merge(
  common_Nterm_full_sphere_exposure,
  how = "left",
  on = ["protein_id", "AA", "position"]
)

common_Nterm_part_sphere_exposure = annotate_accessibility(
    df = common_Nterm_alphafold_annotation,
    max_dist = 12,
    max_angle = 70,
    error_dir = PAE_DIRECTORY
)

common_Nterm_alphafold_accessibility = common_Nterm_alphafold_accessibility.merge(
  common_Nterm_part_sphere_exposure,
  how = "left",
  on = ["protein_id", "AA", "position"]
)

common_Nterm_alphafold_accessibility["high_acc_5"] = np.where(
  common_Nterm_alphafold_accessibility.nAA_12_70_pae <= 5, 1, 0
)

common_Nterm_alphafold_accessibility["low_acc_5"] = np.where(
  common_Nterm_alphafold_accessibility.nAA_12_70_pae > 5, 1, 0
)

# annotate intrinsically disordered region (IDR)
common_Nterm_alphafold_accessibility_smooth = get_smooth_score(
  common_Nterm_alphafold_accessibility,
  np.array(['nAA_24_180_pae']),
  [10]
).reset_index(drop=True)

common_Nterm_alphafold_accessibility_smooth["IDR"] = np.where(
  common_Nterm_alphafold_accessibility_smooth.nAA_24_180_pae_smooth10 <= 34.27, 1, 0
)

# annotate common Nterm data
common_Nterm_alphafold_N_terminus = common_Nterm_alphafold_accessibility_smooth.merge(
  common_Nterm_protein,
  how = "left",
  left_on = ["protein_id", "position"],
  right_on = ["UniProt_Accession", "start.position"]
)

common_Nterm_alphafold_N_terminus = common_Nterm_alphafold_N_terminus.fillna(0)

common_Nterm_alphafold_N_terminus["common_Nterm"] = common_Nterm_alphafold_N_terminus["common_Nterm"].apply(lambda x: 1 if x != 0 else x)

Nterm_site_dict = {
  "common_Nterm" : [
    "A",
    "R",
    "N",
    "D",
    "C",
    "E",
    "Q",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V"
  ]
}

# perform N-terminus enrichment analysis
enrichment_N_terminus = perform_enrichment_analysis(
  df = common_Nterm_alphafold_N_terminus,
  ptm_types = ["common_Nterm"],
  rois = ["BEND", "HELX", "STRN", "TURN", "IDR", "high_acc_5", "low_acc_5"],
  ptm_site_dict = Nterm_site_dict,
  quality_cutoffs = [0]
)

# N-terminus 3D cluster
common_Nterm_proximity = get_proximity_pvals(
    df = common_Nterm_alphafold_N_terminus,
    ptm_types = ['common_Nterm'],
    ptm_site_dict = Nterm_site_dict,
    error_dir = PAE_DIRECTORY,
    per_site_metric = 'mean',
    error_operation = 'plus',
    n_random = 10000,
    random_seed = 44
)
