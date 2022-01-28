
sf_tcr_ful       <- read.delim("data/scRNAseq+TCRseq/Sade_Feldman_tcrab_full.txt") %>% mutate(clonotype_name = paste(CDR3..AA....alpha.or.gamma.chain, CDR3..AA....beta.or.delta.chain, sep = "_"))

## Download the TPM-file from GEO (GSE120575) and wrangle it accodringly
sf               <- fread("data/scRNASeq/raw/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt", fill = T)
cellnames        <- sf[1,-1] %>% as.character()
genenames        <- sf[-c(1:2),1] %>% pull(V1)
sf               <- sf[-c(1:2),-1]
sf               <- sf %>% as.data.frame()

sf               <- sf[,which(cellnames %in% sf_tcr_ful$cell_name)]
colnames(sf)     <- cellnames[which(cellnames %in% sf_tcr_ful$cell_name)]
rownames(sf)     <- genenames
sf               <- Matrix::Matrix(as.matrix(sf), sparse = TRUE)  


## Download the TCRab data and predictions
sf_tcr_tot       <- read.delim("data/scRNAseq+TCRseq/Sade_Feldman_tcrab_clean.txt") 
sf_tcr_pred1     <- getPreds_singlecell(tcrgp_filename = "results/tcrgp/raw/predictions_ver1/scRNAseq+TCRseq/Sade_Feldman_tcrab_clean.txt", thresholds = thresholds_ver1, fdr = 0.05)
sf_tcr_pred2     <- getPreds_singlecell("results/tcrgp/raw/predictions_ver2/scRNAseq+TCRseq/Sade_Feldman_tcrab_clean.txt", thresholds = thresholds_ver2, fdr = 0.05)
sf_tcr_pred      <- cbind(sf_tcr_pred1, dplyr::select(sf_tcr_pred2, ELAGIGILTV_cdr3ab:AMFWSVPTV_cdr3b))

sf_tcr_pred_summ <- summarise_predictions(sf_tcr_tot, sf_tcr_pred) %>% mutate(species = getSpecies(pred_epitope)) %>% mutate(clonotype_name = paste(tra_cdr3aa, trb_cdr3aa, sep = "_"))
sf_tcr_pred_summ <- sf_tcr_pred_summ %>% dplyr::select(trb_cdr3aa,tra_cdr3aa, NLVPMVATV_cdr3ab:clonotype_name)

## Add the cellname information, that is only in the TCRab-file.
sf_tcrgp <- merge(sf_tcr_ful, sf_tcr_pred_summ, all.x = T, by = "clonotype_name")
sf_tcrgp <- sf_tcrgp[!duplicated(sf_tcrgp), ]
sf_tcrgp <- sf_tcrgp %>% mutate(cell_name = as.character(cell_name))

## Add the sure hits from VDJdb
sf_vdjdb  = read.delim("results/vdjdb/sf_tcrb.sf_vdjdb.annot.txt", stringsAsFactors = F) %>% preprocessMultiVDJdb %>% dplyr::select(cdr3, species, antigen.epitope, antigen.gene, antigen.species, vdjdb.score)
colnames(sf_vdjdb)[-1] = paste0("vdjdb.", colnames(sf_vdjdb)[-1])
sf_tcrgp = merge(sf_tcrgp, sf_vdjdb, by.x = "CDR3..AA....beta.or.delta.chain", by.y = "cdr3", all.x = T)


## Select only cells with TCRab
sf_tcrgp_filt = sf_tcrgp %>% filter(cell_name %in% colnames(sf))
sf_tcrgp_filt = sf_tcrgp_filt[match(colnames(sf), sf_tcrgp_filt$cell_name), ]
rownames(sf_tcrgp_filt) = sf_tcrgp_filt$cell_name


#### =============

# Create Seurat
sf_seurat = CreateSeuratObject(sf, meta.data = sf_tcrgp_filt, min.cells = 3)

## ======  Add the proper meta
sf_meta = sf_seurat@meta.data

## Get anti-melanoma predicted cells
viral_species_tcgrp <- c("CMV", "EBV", "INFa")
anti_viral    = do.call("c", lapply(viral_species_tcgrp, function(x) grep(x, sf_meta$species)))
anti_melanoma = do.call("c", lapply(melanoma_epitopes, function(x) grep(x, sf_meta$pred_epitope)))

sf_meta$melanoma_target <- ifelse(sf_meta$pred_epitope %in% c("melana_cdr3b", "ELAGIGILTV_cdr3b_comb", "mart1_cdr3b"), "mart1", "non_melanoma")
sf_meta$melanoma_target <- ifelse(sf_meta$pred_epitope %in% "meloe1_cdr3b", "meloe1", sf_meta$melanoma_target)
sf_meta$melanoma_target <- ifelse(sf_meta$pred_epitope %in% "AMFWSVPTV_cdr3b", "AMFWSVPTV", sf_meta$melanoma_target)
sf_meta$melanoma_target <- ifelse(sf_meta$pred_epitope %in% "FLYNLLTRV_cdr3b", "FLYNLLTRV", sf_meta$melanoma_target)

sf_meta$target = "None"
sf_meta$target[anti_viral] = "anti-viral"
sf_meta$target[anti_melanoma] = "anti-melanoma"

## Add meta for only the HLA-A2 pos cells
hla_a2 = c(grep("hla_a_02_01_01_01", sf_meta$HLA.A...allele.1), grep("hla_a_02_01_01_01", sf_meta$HLA.A...allele.2))
sf_meta$hla_a2_pos = "No"
sf_meta$hla_a2_pos[hla_a2] = "Yes"

## Change the names for clusters
sf_meta$CD8.cluster...2.clusters[sf_meta$CD8.cluster...2.clusters == 1] = "Naive"
sf_meta$CD8.cluster...2.clusters[sf_meta$CD8.cluster...2.clusters == 2] = "Exhausted"

sf_meta$CD8.cluster...6.clusters[sf_meta$CD8.cluster...6.clusters == 1] = "Exhaustion/cycling"
sf_meta$CD8.cluster...6.clusters[sf_meta$CD8.cluster...6.clusters == 2] = "Exhaustion/HSP"
sf_meta$CD8.cluster...6.clusters[sf_meta$CD8.cluster...6.clusters == 3] = "Exhaustion"
sf_meta$CD8.cluster...6.clusters[sf_meta$CD8.cluster...6.clusters == 4] = "Memory/effector1"
sf_meta$CD8.cluster...6.clusters[sf_meta$CD8.cluster...6.clusters == 5] = "EMRA"
sf_meta$CD8.cluster...6.clusters[sf_meta$CD8.cluster...6.clusters == 6] = "Memory/effector2"

colnames(sf_meta)[23] = "orig6"
colnames(sf_meta)[24] = "orig2"

## Add time point
sf_meta$timepoint <- ifelse(grepl("Pre", sf_meta$sample_name), "Pre", "Post")

## Add response
sf_response <- fread("data/scRNAseq+TCRseq/Sade_Feldman_response.txt")
sf_meta <- sf_meta %>% left_join(sf_response)

sf_meta$overall[sf_meta$sample_name == "Post_P4_B"] <- "R"
sf_meta$overall[sf_meta$sample_name == "Pre_P8_B"]  <- "R"
sf_meta$overall[sf_meta$sample_name == "Pre_P33_B"] <- "R"
sf_seurat@meta.data = sf_meta


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


## ===== Normalize, find HVGs and scale 
# sf_seurat <- NormalizeData(sf_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
sf_seurat <- FindVariableFeatures(sf_seurat, selection.method = "vst", nfeatures = 1000, mean.cutoff = c(1, Inf))
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

DimPlot(sf_seurat, reduction = "pca", group.by = "orig2" )
ggsave("results/scrnaseq/sf/pca/orig_2_clusters.pdf", width = 6, height = 4) 

DimPlot(sf_seurat, reduction = "pca", group.by = "orig6" )
ggsave("results/scrnaseq/sf/pca/orig_6_clusters.pdf", width = 6, height = 4)

ElbowPlot(sf_seurat)
nPcs = sum(sf_seurat[["pca"]]@stdev > 2)
nPcs

## UMAP
sf_seurat <- RunUMAP(sf_seurat, dims = 1:nPcs)

table(sf_seurat$orig6)

sf_seurat@reductions$umap@cell.embeddings %>% as.data.frame() %>% bind_cols(sf_seurat@meta.data) %>% 
  ggplot(aes(UMAP_1, UMAP_2, color = orig2)) + geom_point()

sf_seurat@reductions$umap@cell.embeddings %>% as.data.frame() %>% bind_cols(sf_seurat@meta.data) %>% 
  ggplot(aes(UMAP_1, UMAP_2, color = orig6)) + geom_point()

sf_seurat <- sf_seurat %>% getClu
DimPlot(sf_seurat, reduction = "umap", group.by = "orig2")
gsave("results/scrnaseq/sf/umap/orig_2_clusters.pdf", width = 6, height = 4) 

DimPlot(sf_seurat, reduction = "umap", group.by = "orig6" )
ggsave("results/scrnaseq/sf/umap/orig_6_clusters.pdf", width = 6, height = 4) 

print(sf_seurat[["pca"]], dims = 1, nfeatures = 5)
p <- FeaturePlot(sf_seurat, reduction = "pca", features = c("NKG7", "PRF1", "CD38", "HAVCR2", "PDCD1", "IL7R", "CXCR4", "CCR7", "LMNA", "S1PR1"), ncol = 5, cols = c("lightgrey", "salmon"))
ggsave(plot = p, "results/scrnaseq/sf/pca/feature_top_pca.pdf", width = 22, height = 8) 



## Melanoma targeting cells
Idents(sf_seurat) <- sf_seurat$orig6
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
saveRDS(sf_seurat, file = "results/scrnaseq/sf/sf_seurat2.rds")

sf_seurat_loaded <- readRDS("results/scrnaseq/sf/sf_seurat.rds")
sf_seurat <- readRDS("results/scrnaseq/sf/sf_seurat2.rds")

DimPlot(sf_seurat)
sf_seurat = sf_seurat_loaded





## 

sf_seurat@reductions$umap@cell.embeddings %>% as.data.frame() %>% bind_cols(sf_seurat@meta.data) %>% 
  ggplot(aes(UMAP_1, UMAP_2, color = orig2)) + geom_point(size = 0.01) + theme_void() + theme(legend.position = "none") + scale_color_manual(values = c("darkgreen", "darkred"))
ggsave("results/manuscript/figure5/umap_sf.png", width = 3.5, height = 3)


sf_seurat@reductions$umap@cell.embeddings %>% as.data.frame() %>% bind_cols(sf_seurat@meta.data) %>% 
  filter(target != "None") %>% 
  ggplot(aes(UMAP_1, UMAP_2, color = orig2)) + geom_point(size = 0.01) + theme_void() + theme(legend.position = "none") + scale_color_manual(values = c("darkgreen", "darkred")) + facet_grid(timepoint~target+overall)


# FeaturePlot(sf_seurat, features = c("PDCD1", "LAG3", "TCF7", "NKG7", "GZMB", "MKI67"), ncol = 3)
FeaturePlot(sf_seurat, features = c("NKG7", "PRF1", "TCF7", "LAG3", "PDCD1", "MKI67"), ncol = 3, order = T, cols = c("gray90", "steelblue4"), combine = TRUE) & NoAxes()
ggsave("results/manuscript/figure5/feature_umap_sf.png", width = 10, height = 6)


## de genes
sf_seurat$barcode <- colnames(sf_seurat)
cells.to.keep    <- sf_seurat@meta.data %>% filter(target == "anti-melanoma" & overall == "R") %>% pull(barcode)
cells.to.keep1   <- sf_seurat@meta.data %>% filter(target == "anti-melanoma" & overall == "R" & timepoint == "Pre") %>% pull(barcode)
cells.to.keep2   <- sf_seurat@meta.data %>% filter(target == "anti-melanoma" & overall == "R" & timepoint == "Post") %>% pull(barcode)

sf_maa_r_seurat  <- subset(sf_seurat, cells = cells.to.keep)
sf_meta <- sf_seurat@meta.data %>% filter(target == "anti-melanoma" & overall == "R")

rownames(sf_meta) <- make.unique(sf_meta$barcode)
sf_maa_r_seurat@meta.data <- sf_meta

Idents(sf_maa_r_seurat) <- sf_maa_r_seurat$timepoint
sf_maa_prepost_deg <- FindMarkers(sf_maa_r_seurat, ident.1 = "Pre", ident.2 = "Post", test.use = "t", logfc.threshold = 0.01)
sf_maa_prepost_deg$p.adj <- p.adjust(sf_maa_prepost_deg$p_val, method = "BH")
sf_maa_prepost_deg$gene  <- rownames(sf_maa_prepost_deg)

clonality_genes <- getClonalityGenes(sf_maa_prepost_deg)
unwanted_genes <- getUnwantedGenes(sf_maa_prepost_deg)

ggplot(sf_maa_prepost_deg, aes(avg_logFC, -log2(p.adj), color = avg_logFC)) + geom_point(alpha=0.5) + 
  ggrepel::geom_text_repel(data = sf_maa_prepost_deg %>% filter(p.adj < 0.05 & avg_logFC < 0 & !gene %in% c(clonality_genes, unwanted_genes)), aes(avg_logFC, -log2(p.adj), label = gene), fontface = "italic", size = 3.5) + 
  ggrepel::geom_text_repel(data = sf_maa_prepost_deg %>% filter(p.adj < 0.05 & avg_logFC > 0 & !gene %in% c(clonality_genes, unwanted_genes)), aes(avg_logFC, -log2(p.adj), label = gene), fontface = "italic", size = 3.5) + 
  
  geom_hline(yintercept = -log2(0.05), linetype = "dotted") +
  theme_classic(base_size = 18) + 
  # xlim(values = c(-1,1)) +
  
  scale_colour_gradient2(trans = 'reverse') + theme(legend.position = "none") + labs(x = "average logFC", y = "-log2(Padj)")
ggsave("results/scrnaseq/sf/volcano_genes_maa_prepost.png", width = 6, height = 5)
fwrite(sf_maa_prepost_deg, "results/scrnaseq/sf/volcano_genes_maa_prepost.txt", sep = "\t", quote = F, row.names = F)


table(sf_maa_r_seurat$orig2, sf_maa_r_seurat$timepoint)

sf_maa_r_seurat <- subset(sf_seurat, target %in% "anti-melanoma")
sf_seurat$target
table(sf_maa_r_seurat$timepoint)


DimPlot(sf_seurat)


sf_maa_r <- FindMarkers(sf_seurat, cells.1 = cells.to.keep1, cells.2 = cells.to.keep2, test.use = "t")

Seurat::FindMarkers(object = sf_seurat, cells.1 = cells.to.keep1, cells.2. = cells.to.keep2)



sf_seurat@meta.data %>% group_by(target == "anti-melanoma", orig2) %>% summarise(n = n())
sf_seurat@meta.data %>% group_by(target == "anti-viral", orig2) %>% summarise(n = n())

matrix(c(180,104,187,151), ncol = 2) %>% fisher.test()

a <- matrix(c(180,104,1552,1960), ncol = 2) %>% fisher.test() %>% broom::tidy() %>% mutate(cluster = "Exhausted", target = "anti-maa")
b <- matrix(c(187,151,1953,1505), ncol = 2) %>% fisher.test() %>% broom::tidy() %>% mutate(cluster = "Exhausted", target = "anti-viral")

rbind(a,b) %>% 
  ggplot(aes(cluster,-log10(p.value), fill = target)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(values = c("darkred", "darkblue")) + coord_flip() + 
  theme_classic(base_size = 15) + theme(legend.position = "none") + labs(x = "") + geom_hline(yintercept = -log10(0.05), linetype = "dotted")
ggsave("results/scrnaseq/sf/bar_fisher_exhausted.pdf", width = 4, height = 2)

sf_seurat@meta.data %>% group_by(maa = target == "anti-melanoma", orig2) %>% summarise(n = n()) %>% mutate(prop = n / sum(n)) %>% 
  ggplot(aes(maa,prop,fill=orig2)) + geom_bar(stat = "identity")
