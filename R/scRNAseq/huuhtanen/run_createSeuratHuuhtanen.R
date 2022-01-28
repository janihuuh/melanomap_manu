
## Run create seurat files

## huuhtanen
huuhtanen_dirs   <- list.dirs("data/scRNAseq/huuhtanen/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
folders        <- huuhtanen_dirs
scrnaseq_files <- lapply(folders, function(x){message(x); Read10X(data.dir = x) %>% CreateSeuratObject(project = extractSeuratName(x), min.cells = 3, min.features = 200)})
huuhtanen_seurat  <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = extractSeuratName(folders))

## Basic QC
huuhtanen_seurat  <- PercentageFeatureSet(huuhtanen_seurat, pattern = "^MT-", col.name = "percent.mt")
huuhtanen_seurat  <- PercentageFeatureSet(huuhtanen_seurat, pattern = "^RP", col.name = "percent.ribo")
huuhtanen_seurat  <- PercentageFeatureSet(huuhtanen_seurat, features = cycle.genes, col.name = "percent.cycle")
huuhtanen_seurat@meta.data$barcode   <- colnames(huuhtanen_seurat)
huuhtanen_seurat  <- huuhtanen_seurat %>% getQC()

## Put TCR
tot_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/huuhtanen_barcode.txt")
huuhtanen_seurat <- mergeTCRtoSeurat(seurat_object = huuhtanen_seurat, tcr_df = tot_barcode)

## Get latents, use them for UMAPs, clustering
huuhtanen_seurat %>% getScviInput()
latents <- fread("results/scvi/latents_huuhtanen.txt")

## Now run python/run_scvi_example.py to obtain the latent dimensions
huuhtanen_seurat <- huuhtanen_seurat %>% putLatentsSeurat()
huuhtanen_seurat <- huuhtanen_seurat %>% getLatentUMAP() %>% getLatentClustering()

## Scale data
clonality_genes <- getClonalityGenes(huuhtanen_seurat)
unwanted_genes  <- getUnwantedGenes(huuhtanen_seurat)
huuhtanen_seurat <- huuhtanen_seurat %>% preprocessSeurat(cells.to.use = colnames(huuhtanen_seurat))

Idents(huuhtanen_seurat) <- huuhtanen_seurat$RNA_snn_res.0.2
huuhtanen_seurat$cluster <- Idents(huuhtanen_seurat)

## Get singleR
huuhtanen_seurat <- huuhtanen_seurat %>% getSingler()

## Get DE-genes
huuhtanen_cluster_markers <- FindAllMarkers(huuhtanen_seurat, test.use = "t", max.cells.per.ident = 3e3) %>% filter(p_val_adj < 0.05 & ave_logFC > 0)
huuhtanen_deg_markers     <- lapply(unique(huuhtanen_seurat$cluster), getDEGbyClusterLGL)
