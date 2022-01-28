
## Download and process the Jerby-Arnon samples
ja_counts = fread("data/scRNASeq/jerbyarnon_counts.csv")
ja_counts = ja_counts %>% as.data.frame() %>% `row.names<-`(ja_counts$V1) %>% dplyr::select(-V1) %>% as.matrix() %>%
  Matrix::Matrix(sparse = TRUE)

ja_annot = fread("data/scRNASeq/jerbyarnon_cell.annotations.csv")
rownames(ja_annot) = ja_annot$cells

## CreateSeurat
ja_seurat = CreateSeuratObject(ja_counts, meta.data = ja_annot, min.cells = 3, min.features = 200)

## Calculate cell percentages
cycle.genes  <- c("ANLN", "ASPM","BIRC5","CCNA2","CCNB1","CCNB2","CCND1","CD63","CDC20","CDCA8","CDKN3","CENPE","CENPF",
                  "CEP55","CKAP2L","DLGAP5","FOXM1","GTSE1","H2AFZ","HIST1H1B", "HIST1H1C", "HIST1H1D", "HIST1H1E", "HIST1H2AJ",
                  "HIST1H4C", "HJURP", "HMGB1", "HMGB2", "HMMR", "KIF11", "KIF14", "KIF15", "KIF2C", "LMNA",
                  "MCM3", "MKI67", "NCAPG", "NUSAP1", "PCNA", "PLK1", "PRC1", "RRM2", "SMC4", "STMN1", "TK1", "TOP2A", "TPX2", "TUBA1B",
                  "TUBB", "TYMS", "UBE2C")

ja_seurat[["percent.mt"]] <- PercentageFeatureSet(ja_seurat, pattern = "^MT-")
ja_seurat[["percent.rb"]] <- PercentageFeatureSet(ja_seurat, pattern = "^RP")
ja_seurat[["percent.cc"]] <- PercentageFeatureSet(ja_seurat, features = cycle.genes)

## Normalize, find HVGs and scale
ja_seurat <- NormalizeData(ja_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
ja_seurat <- FindVariableFeatures(ja_seurat, selection.method = "vst", nfeatures = 1000, mean.cutoff = c(1, Inf))

## Scale
ja_seurat <- ScaleData(ja_seurat, features = rownames(ja_seurat))

## PCA
ja_seurat <- RunPCA(ja_seurat, features = VariableFeatures(ja_seurat))
ja_seurat <- RunUMAP(ja_seurat, dims = 1:nPcs)

Idents(ja_seurat) = ja_seurat@meta.data$cell.types
