
sf_tcr_ful       <- read.delim("data/scRNAseq+TCRseq/Sade_Feldman_tcrab_full.txt") %>% mutate(clonotype_name = paste(CDR3..AA....alpha.or.gamma.chain, CDR3..AA....beta.or.delta.chain, sep = "_"))

## Download the TPM-file from GEO (GSE120575) and wrangle it accodringly
sf               <- fread("data/scRNASeq/raw/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt", fill = T)
cellnames        <- sf[1,-1] %>% as.character()
genenames        <- sf[-c(1:2),1] %>% pull(V1)
sf               <- sf[-c(1:2),-1]
sf               <- sf %>% as.data.frame()

colnames(sf)     <- cellnames
rownames(sf)     <- genenames
sf               <- Matrix::Matrix(as.matrix(sf), sparse = TRUE)  

sf_meta           <- fread("data/sade-feldman/clusters.txt")[,1:2]
colnames(sf_meta) <- make.names(colnames(sf_meta))

sf_sample_meta    <-  fread("data/scRNAseq+TCRseq/sf_celldata.txt") %>% select(V1:V2, V5:V6)
colnames(sf_sample_meta) <- c("id", "Cell.Name", "Sample", "Status")
sf_response <- fread("data/scRNAseq+TCRseq/Sade_Feldman_response.txt")
sf_meta <- sf_meta %>% left_join(sf_sample_meta) %>% left_join(sf_response, by = c("Sample" = "sample_name"))
sf_meta$cluster <- sf_meta$Cluster.number %>% as.character() %>% as.numeric() %>% as.factor() %>% getSfCluster()

getSfCluster <- function(clusters){
  
  clusters <- plyr::revalue(clusters, replace = c(
    
    "1"  = "1 B-cells" ,
    "2"  = "2 Plasma cells" ,
    "3"  = "3 Monocytes/macrophages" ,
    "4"  = "4 Dendritic cells" ,
    "5"  = "5 Lymphocytes" ,
    "6"  = "6 Exhausted CD8+ T-cells" ,
    "7"  = "7 Regulatory T-cells" ,
    "8"  = "8 Cyotoxic" ,
    "9"  = "9 Exhausted/HS CD8+ T-cells" ,
    "10" = "10 Memory T-cells" ,
    "11" = "11 Lymphocytes exhausted/cycling" 
  ))
  
  return(clusters)
  
}
#### =============

# Create Seurat
sf_seurat = CreateSeuratObject(sf, meta.data = sf_meta)

## ===== Calculate cell percentages
cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                  "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ", 
                  "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA", 
                  "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B", 
                  "TUBB", "TYMS", "UBE2C")

sf_seurat[["percent.mt"]] <- PercentageFeatureSet(sf_seurat, pattern = "^MT-")
sf_seurat[["percent.rb"]] <- PercentageFeatureSet(sf_seurat, pattern = "^RP")
sf_seurat[["percent.cc"]] <- PercentageFeatureSet(sf_seurat, features = cycle.genes)


## No need for QC
# sf_seurat <- subset(sf_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20 & percent.rb < 50)
temp <- readRDS("results/scrnaseq/sf/sf_seurat.rds")
HVFInfo(object = temp)

sf_all_cells <- readRDS("results/scrnaseq/all_sf_cells.rds")
sf_all_cells <- UpdateSeuratObject(sf_all_cells)

sf_all_cells$Cell.Name <- colnames(sf_all_cells)
df <- sf_all_cells@meta.data %>% left_join(sf_meta)
sf_all_cells@meta.data <- df

getHVGscran <- function(seurat_object){
  
  ## @ params
  # object = seurat object
  
  # Using scran to pull the list of hvgs
  scran_object <- as.SingleCellExperiment(seurat_object)
  fit          <- scran::trendVar(scran_object, parametric = TRUE, use.spikes = FALSE)
  hvgs         <- scran::decomposeVar(scran_object, fit)
  return(hvgs)
  
}

## ==== Calculate HVGs with scran 
scran_hvg <- getHVGscran(sf_seurat)
scran_hvg <- scran_hvg %>% as.data.frame() %>% add_rownames(var = "gene")
scran_hvg <- scran_hvg[!scran_hvg$gene %in% clonality_genes, ]
scran_hvg <- scran_hvg[!scran_hvg$gene %in% unwanted_genes, ]

plotScranHVG(scran_hvg, top_n = 20)
ggsave(paste0(output_dir, "scran_hvg.pdf"), width = 8, height = 6)


scran_hvg_sigf <- scran_hvg %>% arrange(FDR) %>% filter(total > 0.5 & FDR < 0.05)

write.table(scran_hvg, paste0(output_dir, "scran_hvg.txt"), sep = "\t", quote = F, row.names = F)
write.table(scran_hvg_sigf, paste0(output_dir, "scran_hvg_sigf.txt"), sep = "\t", quote = F, row.names = F)




## ===== do not normalize, find HVGs and scale 
# sf_seurat <- NormalizeData(sf_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
sf_seurat <- FindVariableFeatures(sf_seurat)

sf_seurat <- FindVariableFeatures(sf_all_cells, selection.method = "vst", nfeatures = 1000, mean.cutoff = c(1, Inf))
hvg.info  <- HVFInfo(object = sf_seurat)  %>% add_rownames(var = "gene") %>% filter(mean > 1 & variance.standardized > 1)

## Remove clonality genes
hvg <- HVFInfo(object = sf_seurat)  %>% add_rownames(var = "gene") %>% filter(mean > 1 & variance.standardized > 1)
hvg <- hvg[!hvg$gene %in% getClonalityGenes(sf_seurat), ]
hvg <- hvg[!hvg$gene %in% getUnwantedGenes(sf_seurat), ]
VariableFeatures(sf_seurat) <- hvg$gene

## Scale
sf_seurat <- ScaleData(sf_seurat, features = rownames(sf_seurat))

## PCA
sf_seurat <- RunPCA(sf_seurat, features = hvg$gene)
ElbowPlot(sf_seurat)
nPcs = sum(sf_seurat[["pca"]]@stdev > 2.5)
nPcs

## UMAP
sf_seurat <- RunUMAP(sf_seurat, dims = 1:nPcs)

sf_seurat@reductions$umap@cell.embeddings %>% plot()
Idents(sf_seurat) <- sf_seurat$clu
sf_seurat <- sf_seurat %>% fixSeurat()
DimPlot(sf_seurat) + theme(legend.position = "none")
DimPlot(sf_seurat, reduction = "umap", group.by = "Cluster.number")

table(is.na(Idents(sf_seurat)))
Idents(sf_seurat) <- Idents(sf_seurat) %>% as.character()
sf_seurat$cluster <- sf_seurat$cluster %>% as.character()
sf_seurat$cluster[is.na(sf_seurat$cluster)] <- "None"
sf_seurat$cluster <- as.factor(sf_seurat$cluster)
Idents(sf_seurat) <- sf_seurat$cluster

DimPlot(sf_seurat)





## Melanoma targeting cells
sf_seurat$cluster <- sf_seurat$Cluster.number %>% as.character() %>% as.numeric() %>% as.factor() %>% getSfCluster()
Idents(sf_seurat) <- sf_seurat$cluster

DimPlot(sf_seurat, group.by = "cluster")


temp$target

sf_seurat$mart1_motif <- ifelse(sf_seurat$trb_cdr3aa %in% matches, T, F)
table(sf_seurat$mart1_motif, sf_seurat$target)

DimPlot(sf_seurat, reduction = "umap", label = T, split.by = "melanoma", group.by = "orig2", cols = getPalette(2))
DimPlot(sf_seurat, reduction = "umap", label = T, split.by = "mart1_motif", group.by = "orig2", cols = getPalette(2))
DimPlot(sf_seurat, reduction = "umap", label = T, split.by = "mart1_motif", group.by = "orig6", cols = getPalette(6))

DimPlot(sf_seurat, reduction = "umap", label = F, group.by = "orig2", cols = getPalette(2)) + theme_void() + add_guide
ggsave("results/manuscript/figure3/umap_sf.pdf", width = 6, height = 4)

DimPlot(sf_seurat, reduction = "umap", split.by = "orig2", label = F, group.by = "target", cols = getPalette(3)) + theme_void() + add_guide



sf_seurat@meta.data %>% group_by(hla_a2_pos, mart1_motif) %>% summarise(n = n())



# ==== Done
saveRDS(sf_seurat, file = "results/scrnaseq/sf/sf_seurat.rds")

sf_seurat_loaded <- readRDS("results/scrnaseq/sf/sf_seurat.rds")
sf_seurat <- readRDS("results/scrnaseq/sf/sf_seurat2.rds") #.rds")
DimPlot(sf_seurat, reduction = "umap") + theme(legend.position = "none")

sf_seurat = sf_seurat_loaded




