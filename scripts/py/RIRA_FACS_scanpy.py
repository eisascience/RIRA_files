# Import required libraries
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import pandas as pd

# Ensure DEV_MODE is correctly defined
DEV_MODE = True

def save_fig(fig, filename):
    # Only save figures if not in development mode
    if not DEV_MODE:
        fig.savefig(filename)

# Set global plot settings for publication quality
mpl.rcParams["savefig.dpi"] = 300  # High resolution
mpl.rcParams["font.size"] = 12     # Adjust font size
mpl.rcParams["figure.figsize"] = (8, 6)  # Set default figure size

# Define output directory for figures
output_dir = "/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/Figs/scanpy"
os.makedirs(output_dir, exist_ok=True)

# Load the saved AnnData object
adata = sc.read_h5ad("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/FACS/RhesusFACS_TNK_Nov3023.seurat_554361.h5ad")

print(adata.obs_names[:5])  # Cell names
print(adata.var_names[:5])  # Gene names

# --- Quality Control ---
# Annotate mitochondrial, ribosomal, and hemoglobin genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")  # Mitochondrial genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))  # Ribosomal genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")  # Hemoglobin genes

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

if not DEV_MODE:
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], 
                 jitter=0.4, multi_panel=True, save="qc_metrics.pdf")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", save="qc_scatter.pdf")

# --- Filtering ---
# Filter cells with minimum gene count and genes expressed in at least 3 cells
# sc.pp.filter_cells(adata, min_genes=100)
# sc.pp.filter_genes(adata, min_cells=3)

# --- Doublet Detection ---
# Run Scrublet to detect doublets
# sc.pp.scrublet(adata)

# Inspect doublet results
# if not DEV_MODE:
#     sc.pl.umap(adata, color=["doublet_score", "predicted_doublet"], save="doublet_detection.pdf")

# Optional: Filter out doublets
# adata = adata[adata.obs["predicted_doublet"] == False]

# --- Normalization ---
# Save raw counts layer
adata.layers["counts"] = adata.X.copy()

# Normalize counts and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# --- Feature Selection ---
# Identify highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
if not DEV_MODE:
    sc.pl.highly_variable_genes(adata, save="hvg.pdf")

print(adata.var["highly_variable"].head())

# Keep only highly variable genes
adata = adata[:, adata.var["highly_variable"]]

# --- Dimensionality Reduction ---
# PCA
sc.tl.pca(adata)
if not DEV_MODE:
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save="pca_variance.pdf")
    sc.pl.pca(adata, color=["pct_counts_mt", "n_genes_by_counts"], dimensions=[(0, 1), (2, 3)], save="pca_dimensions.pdf")

# --- Nearest Neighbor Graph ---
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# --- UMAP Visualization ---
sc.tl.umap(adata)
if not DEV_MODE:
    sc.pl.umap(adata, color="sample", size=2, save="umap_sample.pdf")

# --- Clustering ---
# Run Leiden clustering
sc.tl.leiden(adata, resolution=1.2, flavor="igraph")
if not DEV_MODE:
    sc.pl.umap(adata, color=["leiden"], save="umap_leiden.pdf")

print(adata.obs.columns)

if not DEV_MODE:
    sc.pl.umap(adata, color="Population", size=3, save="umap_population.pdf")

# Explore multiple clustering resolutions
for res in [0.02, 0.5, 1.2, 2.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_res_{res:4.2f}", flavor="igraph")

if not DEV_MODE:
    sc.pl.umap(adata, color=["leiden_res_0.02", "leiden_res_0.50"], 
               legend_loc="on data", save="umap_multiple_resolutions.pdf")
    sc.pl.umap(adata, color=["leiden_res_1.20", "leiden_res_2.00"], 
               legend_loc="on data", save="umap_multiple_resolutions2.pdf")

# --- Save Clustering Results for R ---
# Create a DataFrame with rownames (cell IDs) and clustering assignments
clustering_results = adata.obs[["leiden"] + [f"leiden_res_{res:4.2f}" for res in [0.02, 0.5, 1.2, 2.0]]]
clustering_results.index.name = "barcode"
clustering_results.to_csv("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/FACS/RhesusFACS_TNK_Nov3023_scanpy_processed_clustering_results.csv")

print("Clustering results saved to CSV.")

# Save highly variable gene names for R
hvg_genes = adata.var[adata.var["highly_variable"]].index.tolist()
hvg_genes_df = pd.DataFrame(hvg_genes, columns=["Gene"])
hvg_genes_df.to_csv("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/FACS/RhesusFACS_TNK_Nov3023_hvg_genes.csv", index=False)

# Extract UMAP coordinates
umap_coords = adata.obsm["X_umap"]
# Create a DataFrame with barcodes as rownames and UMAP coordinates
umap_df = pd.DataFrame(umap_coords, index=adata.obs_names, columns=["UMAP1", "UMAP2"])
# Save UMAP coordinates to CSV
umap_df.to_csv("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/FACS/RhesusFACS_TNK_Nov3023_umap_coordinates.csv")
print("UMAP coordinates saved to CSV.")

# Ensure cell names are stored correctly
adata.obs_names.name = "barcode"
# Ensure gene names are stored correctly
adata.var_names.name = "gene"
# Verify that the names are correctly set
print("Cell names (obs_names):", adata.obs_names.name)
print("Gene names (var_names):", adata.var_names.name)

# Save the updated AnnData
adata.write("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/FACS/RhesusFACS_TNK_Nov3023.seurat_554361.h5ad")

print("Pipeline completed. Processed data saved.")
