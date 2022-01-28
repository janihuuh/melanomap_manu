
## Run create seurat files

## Durante
durante_dirs   <- list.dirs("data/scRNAseq/durante/", recursive = T) %>% grep(pattern = "filtered_feature_bc_matrix", value = T)
folders        <- durante_dirs
scrnaseq_files <- lapply(folders, function(x){message(x); Read10X(data.dir = x) %>% CreateSeuratObject(project = extractSeuratName(x), min.cells = 3, min.features = 200)})
durante_seurat  <- merge(scrnaseq_files[[1]], scrnaseq_files[-1], add.cell.ids = extractSeuratName(folders))

## Basic QC
durante_seurat  <- PercentageFeatureSet(durante_seurat, pattern = "^MT-", col.name = "percent.mt")
durante_seurat  <- PercentageFeatureSet(durante_seurat, pattern = "^RP", col.name = "percent.ribo")
durante_seurat  <- PercentageFeatureSet(durante_seurat, features = cycle.genes, col.name = "percent.cycle")
durante_seurat@meta.data$barcode   <- colnames(durante_seurat)
durante_seurat  <- durante_seurat %>% getQC()

## Put TCR
tot_barcode <- fread("data/scRNAseq+TCRseq/preprocessed/durante_barcode.txt")
durante_seurat <- mergeTCRtoSeurat(seurat_object = durante_seurat, tcr_df = tot_barcode)

## Get latents, use them for UMAPs, clustering
durante_seurat %>% getScviInput()
latents <- fread("results/scvi/latents_durante.txt")

## Now run python/run_scvi_example.py to obtain the latent dimensions
durante_seurat <- durante_seurat %>% putLatentsSeurat()
durante_seurat <- durante_seurat %>% getLatentUMAP() %>% getLatentClustering()

## Scale data
clonality_genes <- getClonalityGenes(durante_seurat)
unwanted_genes  <- getUnwantedGenes(durante_seurat)
durante_seurat <- durante_seurat %>% preprocessSeurat(cells.to.use = colnames(durante_seurat))

Idents(durante_seurat) <- durante_seurat$RNA_snn_res.0.2
durante_seurat$cluster <- Idents(durante_seurat)

## Get singleR
durante_seurat <- durante_seurat %>% getSingler()

## Get DE-genes
durante_cluster_markers <- FindAllMarkers(durante_seurat, test.use = "t", max.cells.per.ident = 3e3) %>% filter(p_val_adj < 0.05 & ave_logFC > 0)
durante_deg_markers     <- lapply(unique(durante_seurat$cluster), getDEGbyClusterLGL)
