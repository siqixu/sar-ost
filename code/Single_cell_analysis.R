library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)  


single_cell_data <- readRDS('/singlecell/SKM_human_pp/SKM_human_pp.rds')
single_cell_data <- subset(single_cell_data, batch == 'cells')

# --------------------QC -----------------------------------------

single_cell_data[["percent.mt"]] <- PercentageFeatureSet(single_cell_data, pattern = "^MT-")
single_cell_data[["percent.rb"]] <- PercentageFeatureSet(single_cell_data, pattern = "^RP[SL]")

HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(single_cell_data))
single_cell_data[["percent.HB"]] <- PercentageFeatureSet(single_cell_data, features = HB.genes)


qc_table <- single_cell_data@meta.data[, c('percent.mt', "percent_mito", 'percent.rb',
                               "percent_ribo", "percent.HB")]

theme.set2 <- theme(axis.title.x = element_blank())
plot.features <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent.HB",
                   "percent_ribo")
group <- "SampleID"

plots_before <- list()
for (i in seq_along(plot.features)) {
  plots_before[[i]] <- VlnPlot(single_cell_data, group.by = group, pt.size = 0,
                               features = plot.features[i]) +
                       theme.set2 + NoLegend()
}
violin_before <- wrap_plots(plots = plots_before, nrow = 2)


minGene <- 300
maxGene <- 10000
minUMI  <- 600
pctMT   <- 10    
pctHB   <- 1     

single_cell_data <- subset(single_cell_data,
               subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
               nCount_RNA > minUMI & percent.mt < pctMT & percent.HB < pctHB)

plots_after <- list()
for (i in seq_along(plot.features)) {
  plots_after[[i]] <- VlnPlot(single_cell_data, group.by = group, pt.size = 0,
                              features = plot.features[i]) +
                      theme.set2 + NoLegend()
}
violin_after <- wrap_plots(plots = plots_after, nrow = 2)
print(violin_after)


single_cell_data <- subset(single_cell_data, batch == 'cells')

# ------------------ PCA  --------------------------
single_cell_data <- NormalizeData(single_cell_data)           
single_cell_data <- FindVariableFeatures(single_cell_data)   
single_cell_data <- ScaleData(single_cell_data)               
single_cell_data <- RunPCA(single_cell_data)                  

DimPlot(single_cell_data, reduction = 'pca', group.by = 'SampleID')
DimPlot(single_cell_data, reduction = 'pca', group.by = 'batch')
DimPlot(single_cell_data, reduction = 'pca', group.by = 'Age_bin')


ElbowPlot(single_cell_data, ndims = 50)
pct  <- single_cell_data[["pca"]]@stdev / sum(single_cell_data[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
pcs <- 1:41  

# ----------------- clustering -----------------------------------------
single_cell_data <- RunHarmony(single_cell_data,
                   group.by.vars = "SampleID",
                   assay.use = "RNA",
                   max.iter.harmony = 20)
 
resolutions <- seq(0.1, 2, by = 0.1)
single_cell_data <- FindNeighbors(single_cell_data, dims = pcs)

for (res in resolutions) {
  single_cell_data <- FindClusters(single_cell_data, resolution = res, verbose = FALSE)
}
 
single_cell_data <- FindNeighbors(single_cell_data, reduction = "harmony", dims = pcs) %>%
        FindClusters(resolution = 1)
single_cell_data <- RunUMAP(single_cell_data, reduction = "harmony", dims = pcs) %>%
        RunTSNE(dims = pcs, reduction = "harmony", check_duplicates = FALSE)

# ---------------------- annotation ---------------------------
FeaturePlot(single_cell_data, features = c('PAX7','DLK1','MYF5'), ncol = 2)   # Muscle stem cell
FeaturePlot(single_cell_data, features = c('MYH7','TNNT1','TNNI1'), ncol = 2)  # Type I muscle fibers
FeaturePlot(single_cell_data, features = c('TNNT3','MYH1','MYH2','TNNI2'), ncol = 3) # Type II muscle fibers
FeaturePlot(single_cell_data, features = c('CKM','MB','TTN','ACTA1'), ncol = 2)      # Skeletal muscles

FeaturePlot(single_cell_data, features = c('PTPRC'), ncol = 2)               # Immune cell
FeaturePlot(single_cell_data, features = c('CD3D','IL7R'), ncol = 2)         # T cell
FeaturePlot(single_cell_data, features = c('NKG7','PRF1'), ncol = 2)         # NK cell
FeaturePlot(single_cell_data, features = c('IGHM','CD79A','MS4A1'), ncol = 2) # B cell

FeaturePlot(single_cell_data, features = c('FCGR3B','CSF3R','SORL1'), ncol = 2)   # Neutrophils
FeaturePlot(single_cell_data, features = c('TPSB2','MS4A2'), ncol = 2)            # mast cell
FeaturePlot(single_cell_data, features = c('S100A8','S100A12','CD163'), ncol = 2) # monocytes
FeaturePlot(single_cell_data, features = c('CD163','C1QA','CD14'), ncol = 2)      # macrophages

FeaturePlot(single_cell_data, features = c('PDGFRA','DCN'), ncol = 2) # Stromal cell
FeaturePlot(single_cell_data, features = c('TAGLN','MYL9','NOTCH3','ACTA2'), ncol = 2) # Mural cell
FeaturePlot(single_cell_data, features = c('EMCN','CDH5','VWF'), ncol = 2) # Endothelial cell
FeaturePlot(single_cell_data, features = c('CD14','CD68','CD163','S100A9','S100A8'), ncol = 3) # Myeloid cell

# others
FeaturePlot(single_cell_data, features = c('PLP1','SFRP5','CDH19'), ncol = 2)     
FeaturePlot(single_cell_data, features = c('HPR'), ncol = 1)                      
 
DimPlot(single_cell_data, reduction = "umap", label = TRUE)
single_cell_data$celltype.2 <- recode(single_cell_data$seurat_clusters,
  "0"  = "MuSc",           "1"  = "Stromal cell",
  "2"  = "Type II",        "3"  = "Stromal cell",
  "4"  = "Type I",         "5"  = "Mural cell",
  "6"  = "MuSc",           "7"  = "Endothelial cell",
  "8"  = "T cell",         "9"  = "Stromal cell",
  "10" = "NK cell",        "11" = "Myeloid cell",
  "12" = "Mural cell",     "13" = "Endothelial cell",
  "14" = "Myeloid cell",   "15" = "T cell",
  "16" = "Mural cell",     "17" = "B cell",
  "18" = "Stromal cell",   "19" = "Endothelial cell",
  "20" = "Myeloid cell",   "21" = "Myeloid cell",
  "22" = "Endothelial cell","23" = "Other",
  "24" = "Type II",        "26" = "Mast cell",
  "25" = "Mural cell",     "27" = "Stromal cell",
  "28" = "B cell",         "29" = "Other",
  "30" = "Myeloid cell",   "31" = "Myeloid cell",
  "32" = "Stromal cell"
)

DimPlot(single_cell_data, reduction = "umap", label = TRUE, group.by = "celltype.2")
saveRDS(single_cell_data, '/single_cell_re/nature_aing_celltype.rds')


# --------------- Expression of mediators -----------------------------

single_cell_data <- SetIdent(single_cell_data, value = single_cell_data@meta.data$celltype.2)

protein <- fread('/single_cell_re/UWP_OP.csv')
gene_names <- protein$Protein

available_genes <- intersect(gene_names, rownames(single_cell_data))
gene_names <- available_genes
cells <- unique(single_cell_data@meta.data$celltype.2)

for (cell in cells) {
  results_df <- data.frame(
    gene_name = character(),
    target_cell = character(),
    expression_ratio = numeric(),
    target_mean = numeric(),
    other_mean = numeric(),
    diff_mean = numeric(),
    stringsAsFactors = FALSE
  )
  
  
  target_cell_names <- WhichCells(single_cell_data, ident = cell)
  other_cell_names <- WhichCells(single_cell_data, ident = setdiff(unique(single_cell_data@active.ident), cell))
  
  for (gene_name in gene_names) {
    
    target_expr <- FetchData(single_cell_data, vars = gene_name, cells = target_cell_names)
    expr_ratio <- sum(target_expr[, 1] > 0, na.rm = TRUE) / nrow(target_expr)
    target_mean <- mean(target_expr[, 1], na.rm = TRUE)
    
    other_expr <- FetchData(single_cell_data, vars = gene_name, cells = other_cell_names)
    other_mean <- mean(other_expr[, 1], na.rm = TRUE)
    
    results_df <- rbind(results_df, data.frame(
      gene_name = gene_name,
      target_cell = cell,
      expression_ratio = expr_ratio,
      target_mean = target_mean,
      other_mean = other_mean,
      stringsAsFactors = FALSE
    ))
  }
  results_df$diff_mean <- results_df$target_mean / results_df$other_mean
  
  write.csv(results_df, 
            paste0('/single_cell_re/', cell, '_exp_UWP.csv'),
            row.names = FALSE)
}
