
#these load the libraries and files needed 
source("./rbin/RIRA4v3/LoadLibraries.R")
source("./rbin/FACSdata/FACS_Seurat_final_load.R")

#more specific libraries
library(reticulate)
library(zellkonverter)
library(SeuratDisk)
library(S4Vectors)
library(SingleCellExperiment)

ComboSerObj$barcode = colnames(ComboSerObj)

#thse files are made in python... see the python code
sce_Meta = data.table::fread("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/FACS/RhesusFACS_TNK_Nov3023_scanpy_processed_clustering_results.csv") %>% as.data.frame()
sce_top2K = read.csv("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/FACS/RhesusFACS_TNK_Nov3023_hvg_genes.csv")[,1]
sce_UMAP = read.csv("/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/data/FACS/RhesusFACS_TNK_Nov3023_umap_coordinates.csv")
rownames(sce_Meta) = sce_Meta$barcode
sce_Meta = sce_Meta[,-1]

rownames(sce_UMAP) = sce_UMAP$barcode
sce_UMAP = sce_UMAP[,-1]

# add the UMAP from scanpy in our seurat
ComboSerObj[["umap_sc"]] <- CreateDimReducObject(
  embeddings = as.matrix(sce_UMAP),
  key = "UMAP_",
  assay = DefaultAssay(ComboSerObj)
)

# add the metadata 
ComboSerObj = AddMetaData(ComboSerObj, sce_Meta)



# Process ComboSerObj with default 2000 top variable genes
ComboSerObj_top2K = FindVariableFeatures(ComboSerObj, selection.method = "vst", nfeatures = 2000)

# Venn Diagram comparison of genes
# library(VennDiagram)
genes_supervised = VariableFeatures(ComboSerObj)
genes_scanpy = sce_top2K
genes_combo = VariableFeatures(ComboSerObj_top2K)


## some venn diagrams 
ggvenn_plot <- ggvenn::ggvenn(
  list(#RedSupervised = genes_supervised, 
       Top2K_ScanPy = genes_scanpy, 
       Top2K_Seurat = genes_combo),
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"),
  stroke_size = 0.5,
  text_size = 6.5,
  set_name_size = 6.5
); ggvenn_plot


ggplot2::ggsave(
  filename = "/Volumes/Maggie/Work/OHSU/Bimber/Expts/RIRA_manuscript/Figs/scanpy/venn_comparison2.pdf",
  plot = ggvenn_plot,
  width = 9, height = 7, units = "in", dpi = 300
)



# Process ComboSerObj with default 2000 top variable genes

ComboSerObj_top2K <- ScaleData(ComboSerObj_top2K, verbose = FALSE)
ComboSerObj_top2K <- RunPCA(ComboSerObj_top2K, npcs = 50, verbose = FALSE)
ElbowPlot(ComboSerObj_top2K, ndims = 50)
ComboSerObj_top2K <- FindNeighbors(ComboSerObj_top2K, dims = 1:15, verbose = FALSE)
ComboSerObj_top2K = RunUMAP(ComboSerObj_top2K, reduction.name = "umap_2k",reduction = "pca", dims = 1:15)



# Clustering comparisons using ComboSerObj_top2K
## 2.1 Run Scanpy-like Leiden clustering (through Seurat)
ComboSerObj_top2K <- FindClusters(ComboSerObj_top2K, resolution = 1.2, method = "igraph", 
                                  algorithm = 4, verbose = FALSE,
                                  cluster.name = "RNA_Leiden_res.1.2")
Idents(ComboSerObj_top2K) = "RNA_Leiden_res.1.2"
# RedFACSgenes_Leid_clusters <- Idents(ComboSerObj_top2K)

## 2.2 Louvain algorithm with multilevel refinement
ComboSerObj_top2K <- FindClusters(ComboSerObj_top2K, resolution = 1.2, method = "igraph", 
                                  algorithm = 2, verbose = FALSE,
                                  cluster.name = "RNA_Louv_refined_res.1.2")
Idents(ComboSerObj_top2K) = "RNA_Louv_refined_res.1.2"
# RedFACSgenes_LouvRef_clusters <- Idents(ComboSerObj_top2K)

## 2.3 SLM algorithm
ComboSerObj_top2K <- FindClusters(ComboSerObj_top2K, resolution = 1.2, method = "igraph", 
                                  algorithm = 3, verbose = FALSE,
                                  cluster.name = "RNA_SLM_res.1.2")
Idents(ComboSerObj_top2K) = "RNA_SLM_res.1.2"
# RedFACSgenes_SLM_clusters <- Idents(ComboSerObj_top2K)


## 2.4 defualt leiden algorithm
ComboSerObj_top2K <- FindClusters(ComboSerObj_top2K, resolution = 1.2, #method = "igraph", 
                                  algorithm = 1, verbose = FALSE,
                                  cluster.name = "RNA_res.1.2")
Idents(ComboSerObj_top2K) = "RNA_res.1.2"
# RedFACSgenes_SLM_clusters <- Idents(ComboSerObj_top2K)



## 3.0 Sankey Plot Visualization with highcharter
library(highcharter)
library(RColorBrewer)


# Prepare data for Sankey plot (clusters only)
sankey_clusters <- data.frame(
  Louvain_2KSeurat_Clusters = paste0("Louv_ser2k_", ComboSerObj_top2K$RNA_res.1.2),
  FACS_2KSeurat_Population = paste0("FACS_", ComboSerObj_top2K$Population),
  Leiden_2KScanPy_Clusters = paste0("Leid_sc2k_", ComboSerObj_top2K$leiden_res_1.20)
)

hchart(data_to_sankey(sankey_clusters), "sankey", name = "Comparison of Cluster Results")


# Prepare data for Sankey plot (clusters only)
sankey_clusters <- data.frame(
  FACS_2KSeurat_Population = paste0("FACS_", ComboSerObj_top2K$Population),
  Louvain_2KSeurat_Clusters = paste0("Louv_ser2k_", ComboSerObj_top2K$RNA_res.1.2)
)

hchart(data_to_sankey(sankey_clusters), "sankey", name = "Comparison of Cluster Results")



# Prepare data for Sankey plot (clusters only)
sankey_clusters <- data.frame(
  FACS_2KSeurat_Population = paste0("FACS_", ComboSerObj_top2K$Population),
  Leiden_2KScanPy_Clusters = paste0("Leid_sc2k_", ComboSerObj_top2K$leiden_res_1.20)
)

hchart(data_to_sankey(sankey_clusters), "sankey", name = "Comparison of Cluster Results")




mytbl = round(table(paste0("Louv_ser2k_", ComboSerObj_top2K$RNA_res.1.2),
                    paste0("Leid_sc2k_", ComboSerObj_top2K$leiden_res_1.20)), 2)

breaksList = seq(0, max(mytbl), by = 10)
pheatmap::pheatmap(mytbl,
                   color = colorRampPalette((brewer.pal(n = 7, name = "BrBG")[4:7]))(length(breaksList)),
                   # breaks = breaksList, 
                   fontsize = 14, display_numbers = T,
                   number_format = "%.0f", fontsize_number = 12, 
                   clustering_method = "ward.D2", 
                   clustering_distance_rows = "correlation", 
                   clustering_distance_cols = "correlation")



 # Prepare data for Sankey plot (clusters only)

sankey_clusters <- data.frame(
  Louvain_2KSeurat_Clusters = paste0("Louv_ser2k_", ComboSerObj_top2K$RNA_res.1.2),
  Louvain_175Seurat_Clusters = paste0("Louv_ser175_", ComboSerObj$RNA_snn_res.1.2),
  LouvRefined_2KSeurat_Clusters = paste0("LouvRef_ser2k_", ComboSerObj_top2K$RNA_Louv_refined_res.1.2),
  SLM_2KSeurat_Clusters = paste0("SLM_ser2k_", ComboSerObj_top2K$RNA_SLM_res.1.2),
  Leiden_2KScanPy_Clusters = paste0("Leid_sc2k_", ComboSerObj_top2K$leiden_res_1.20)
)

 
 hchart(data_to_sankey(sankey_clusters), "sankey", name = "Comparison of Cluster Results")
 
 

#supervised 175 gene clusters
DimPlot(ComboSerObj, group.by = "RNA_snn_res.1.2", label = T, label.size = 8) + NoLegend()
DimPlot(ComboSerObj_top2K, group.by = "leiden_res_1.20", label = T, label.size = 8) + NoLegend()
DimPlot(ComboSerObj_top2K, group.by = "RNA_Leiden_res.1.2", label = T, label.size = 8) + NoLegend()
DimPlot(ComboSerObj_top2K, group.by = "RNA_Louv_refined_res.1.2", label = T, label.size = 8) + NoLegend()
DimPlot(ComboSerObj_top2K, group.by = "RNA_SLM_res.1.2", label = T, label.size = 8) + NoLegend()




# UMAP visualizations
p1 <- DimPlot(ComboSerObj, group.by = "Population", label = F, label.size = 6) + NoLegend() + ggtitle("175 Supervized Gene Space")
p2 <- DimPlot(ComboSerObj, group.by = "Population", reduction = "umap_sc", label = F, label.size = 6) + NoLegend() + ggtitle("2k Unsupervized Scanpy")
p3 <- DimPlot(ComboSerObj_top2K, group.by = "Population", reduction = "umap_2k", label = F, label.size = 6) + ggtitle("2k Unsupervized Seurat")

p1 | p2 | p3


DimPlot(ComboSerObj, group.by = "Pop2", label = F, label.size = 8, cols = col_vector)  + ggtitle("FACS sorted")



# Visualize the Scanpy clustering results
DimPlot(ComboSerObj, group.by = "leiden", reduction = "umap_sc", label = TRUE) +
  ggtitle("Scanpy Leiden Clustering")

DimPlot(ComboSerObj, group.by = "Population", label = T, label.size = 8) + NoLegend() |
DimPlot(ComboSerObj, group.by = "Population", reduction = "umap_sc", label = T, label.size = 8) + NoLegend()

DimPlot(ComboSerObj, group.by = "Pop2", label = T, label.size = 8) + NoLegend() + facet_wrap(~Pop2)  + ggtitle("FACS sorted")
DimPlot(ComboSerObj, group.by = "Population", label = F, label.size = 6) + facet_wrap(~Population) + NoLegend()+ ggtitle("FACS sorted")
