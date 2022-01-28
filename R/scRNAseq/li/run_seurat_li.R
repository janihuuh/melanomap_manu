
li_counts <- readMM("data/scRNASeq/Li2018_umi_counts/Li2018_matrix.mtx")
li_genes  <- fread("data/scRNASeq/Li2018_umi_counts/Li2018_genes.tsv", header = F)
li_cells  <- fread("data/scRNASeq/Li2018_umi_counts/Li2018_cells.tsv", header = F)
li_meta   <- fread("data/scRNASeq/Li2018_umi_counts/Li2018_metadata.tsv")
rownames(li_meta) <- li_cells$V1

rownames(li_counts) <- li_genes$V1
colnames(li_counts) <- li_cells$V1
li_seurat <- CreateSeuratObject(counts = li_counts, meta.data = li_meta, project = "li")
li_seurat$barcode <- colnames(li_seurat)

## Make QC, seurat object
li_seurat <- li_seurat %>% preprocessSeurat(cells.to.use = colnames(li_seurat))
li_seurat <- li_seurat %>% getClustering()
li_seurat %>% plotClustering

Idents(li_seurat) <- as.factor(as.numeric(as.character(li_seurat$RNA_snn_res.0.3))) %>% getClusterPhenotypesLiTotal
li_seurat$cluster <- Idents(li_seurat)

## Add tcrgp predictions
li_tcrgp <- fread("results/scrnaseq/li/li_tcrgp.txt")
df <- li_tcrgp %>% dplyr::select(Well_ID, pred_epitope, target)
df <- li_seurat@meta.data %>% left_join(df, by = c("barcode" = "Well_ID"))
rownames(df) <- df$barcode
li_seurat@meta.data <- df
