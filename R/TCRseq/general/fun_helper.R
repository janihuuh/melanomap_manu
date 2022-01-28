
getUnwantedGenes <- function(object){
  
  unwanted_variation <- c(grep("^LINC", rownames(object), value = T), grep("^AC", rownames(object), value = T),
                          grep("^AL", rownames(object), value = T),
                          grep("^MT-", rownames(object), value = T), grep("^RP", rownames(object), value = T))
  
}

getClonalityGenes <- function(object){
  
  clonality_genes <- c(grep("^TRAV", rownames(object), value = T), grep("^TRBV", rownames(object), value = T),
                       grep("^TRGV", rownames(object), value = T), grep("^TRDV", rownames(object), value = T),
                       grep("^IGLV", rownames(object), value = T), grep("^IGLC", rownames(object), value = T),
                       grep("^IGLL", rownames(object), value = T), grep("^IGKV", rownames(object), value = T),
                       grep("^IGHV", rownames(object), value = T), grep("^IGKC", rownames(object), value = T),
                       grep("^IGH", rownames(object), value = T),  grep("^IGK", rownames(object), value = T))
  
}

preprocessSeurat <- function(orig_object){
  
  ## Subset object
  cells.to.use <- colnames(orig_object)
  object <- subset(orig_object, cells = cells.to.use)
  clonality_genes <- getClonalityGenes(object)
  unwanted_genes <- getUnwantedGenes(object)
  
  # orig_object@meta.data$barcode
  temp_meta <- orig_object@meta.data[as.character(orig_object@meta.data$barcode) %in% cells.to.use, ]
  temp_meta <- temp_meta[match(colnames(object), temp_meta$barcode), ]
  temp_meta$barcode == colnames(object)
  object@meta.data <- temp_meta
  
  ## Normalize and find HVGs
  object  <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object  <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000, clip.max = 10)
  
  ## Remove clonality genes
  hvg     <- VariableFeatures(object)
  too_hvg <- HVFInfo(object = object) %>% add_rownames(var = "gene") %>% filter(variance.standardized > 10) %>% pull("gene") %>% as.character()
  hvg     <- hvg[!hvg %in% too_hvg]
  hvg     <- hvg[!hvg %in% clonality_genes]
  hvg     <- hvg[!hvg %in% unwanted_genes]
  
  VariableFeatures(object) <- hvg
  # plotHVG(object, 30) #+ ylim(values = c(0,10))
  
  ## Scale data
  object <- ScaleData(object, features = hvg)
  
  ## PCA data
  object <- RunPCA(object, features = hvg, npcs = 50)
  nPCs   <- sum(object[["pca"]]@stdev > 2)
  print(paste("nPCs:", nPCs))
  
  ## RunUMAP does not work
  object <- RunUMAP(object, dims = 1:nPCs, learning.rate = 1)
  
  # Meanwhile try something hacky-ish
  # umap_df <- object[["pca"]]@cell.embeddings[,1:nPCs] %>% umapr::umap() %>% select(UMAP1:UMAP2)
  # umap_df <- CreateDimReducObject(key = "umap", embeddings = as.matrix(x = umap_df))
  # object[["umap"]] <- umap_df
  
  return(object)
  
}

extractClusterNumber <- function(strs){
  
  p <- NULL
  i <- 1
  for(str1 in strs){
    p[[i]] <- strsplit(str1, "[ ]")[[1]][1]
    i <- i + 1
  }
  
  return(p)
  
}


getSingler <- function(seurat_object, cluster = NULL, method = NULL, sample = NULL){
  
  hpca.se   <- SingleR::HumanPrimaryCellAtlasData()
  blueprint <- SingleR::BlueprintEncodeData()
  ## @ params
  ## cluster = possible cluster vec, if not provided, tries to find in meta.data$cluster
  ## method = if "cluster", then performs preds based on clusters, not cells
  ## sample = to subsample or not
  
  if(!is.null(sample)){
    
    set.seed(123)
    seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[sample(1:ncol(seurat_object), sample)])
    
  }
  
  sce       <- as.SingleCellExperiment(seurat_object)
  
  ## Predictions
  if(is.null(method)){
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine)
    
    if(is.null(sample)){
      seurat_object$singler_hpca_pred      <- pred.hca$first.labels
      seurat_object$singler_blueprint_pred <- pred.blu$first.labels
      return(seurat_object)
    }
    
    else{
      df <- data.frame(barcode = rownames(pred.hca), cluster = seurat_object$cluster, singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
      return(df)
    }
    
  }
  
  
  if(method == "cluster"){
    if(is.null(cluster)){
      cluster=seurat_object$cluster
    }
    pred.hca <- SingleR::SingleR(test = sce, ref = hpca.se, assay.type.test = 1,   labels = hpca.se$label.fine, method = "cluster", clusters = cluster)
    pred.blu <- SingleR::SingleR(test = sce, ref = blueprint, assay.type.test = 1, labels = blueprint$label.fine, method = "cluster", clusters = cluster)
    df <- data.frame(cluster = rownames(pred.hca), singler_hpca_pred = pred.hca$labels, singler_blueprint_pred = pred.blu$labels)
    return(df)
  }
}


getClonalityGenes <- function(object){
  
  clonality_genes <- c(grep("^TRAV", rownames(object), value = T), grep("^TRBV", rownames(object), value = T),
                       grep("^TRGV", rownames(object), value = T), grep("^TRDV", rownames(object), value = T),
                       
                       grep("^TRAJ", rownames(object), value = T), grep("^TRBJ", rownames(object), value = T),
                       grep("^TRGJ", rownames(object), value = T), grep("^TRDJ", rownames(object), value = T),
                       
                       grep("^IGLV", rownames(object), value = T), grep("^IGLC", rownames(object), value = T),
                       grep("^IGLL", rownames(object), value = T), grep("^IGKV", rownames(object), value = T),
                       grep("^IGHV", rownames(object), value = T), grep("^IGKC", rownames(object), value = T),
                       grep("^IGH", rownames(object), value = T),  grep("^IGK", rownames(object), value = T))
  
}



getUnwantedGenes <- function(object){
  
  unwanted_variation <- c(grep("^LINC", rownames(object), value = T), grep("^AC", rownames(object), value = T),
                          grep("^AL", rownames(object), value = T),
                          grep("^MT-", rownames(object), value = T), grep("^RP", rownames(object), value = T))
  
}




fixSeurat <- function(seurat_object){
  
  ## Fix meta data if it brokes
  
  meta.data           <- seurat_object@meta.data
  count.data          <- seurat_object@assays$RNA@counts
  scale.data          <- seurat_object@assays$RNA@scale.data
  # hvg                 <- VariableFeatures(seurat_object)
  
  # pca_dimred          <- seurat_object[["pca"]]
  # umap_dimred         <- seurat_object[["umap"]]
  latent_dimred       <- seurat_object[["latent"]]
  latent_umap_dimred  <- seurat_object[["latent_umap"]]
  
  rownames(meta.data) <- meta.data$barcode
  
  old_idents <- Idents(seurat_object)
  new_seurat <- CreateSeuratObject(counts = count.data)
  
  new_seurat@meta.data             <- meta.data
  new_seurat@assays$RNA@counts     <- count.data
  new_seurat@assays$RNA@scale.data <- scale.data
  # VariableFeatures(seurat_object)  <- hvg
  
  # new_seurat[["pca"]]              <- pca_dimred
  # new_seurat[["umap"]]             <- umap_dimred
  new_seurat[["latent"]]           <- latent_dimred
  new_seurat[["latent_umap"]]      <- latent_umap_dimred
  Idents(new_seurat) <- old_idents
  return(new_seurat)
  
}

getLatentUMAP <- function(seurat_object){
  
  umap_df           <- seurat_object[["latent"]]@cell.embeddings %>% uwot::umap()
  colnames(umap_df) <- c("latent_umap1", "latent_umap2")
  rownames(umap_df) <- colnames(seurat_object)
  umap_df           <- CreateDimReducObject(key = "latent_umap", embeddings = as.matrix(x = umap_df))
  seurat_object[['latent_umap']] <- umap_df
  
  return(seurat_object)
  
}


## Helper function to break names into clinical data
breakName <- function(name_temp){

  ## Breaks the filename into relevant information

  # "2121_PR_PB_LAG3_R_0m_MNC.txt"

  ## Emerson samples
  # if(substr(name_temp, 1, 2) %in% c("HI", "Ke")){
  #   
  #   for(i in 1:length(name_temp)){
  # 
  #     name         <- name_temp
  #     response[i]  <- "CT"
  #     type[i]      <- "PB"
  #     regimen[i]   <- "CTRL"
  #     overall[i]   <- "C"
  #     timepoint[i] <- "0m"
  #     celltype[i]  <- "MNC"
  #     io_stat[i]   <- "C"
  #     
  #   }
  #   
  #   return(data.frame(name, response, type, regimen, overall, timepoint, celltype, io_stat, stringsAsFactors = F))
  #   
  # }
  
  
  name      <- substr(name_temp, 1, 4)
  response  <- substr(name_temp, 6, 7)
  type      <- substr(name_temp, 9, 10)
  regimen   <- substr(name_temp, 12, 15)
  overall   <- substr(name_temp, 17, 17)
  timepoint <- substr(name_temp, 19, 20)
  celltype  <- substr(name_temp, 22, 24)
  io_stat   <- substr(name_temp, 26, 27)


# 
  
  for(i in 1:length(type)){
    if(type[i] == "PB"){type[i] = "Blood"}
    if(type[i] == "TX"){type[i] = "Tumor"}
    if(substr(name[i], 1, 2) %in% c("HI", "Ke")){type[i] = "Blood"}
  }

  for(i in 1:length(regimen)){
    if(regimen[i] == "CTLA"){regimen[i] = "antiCTLA4"}
    if(regimen[i] == "iPD1"){regimen[i] = "antiPD1"}
    if(regimen[i] == "NiIp"){regimen[i] = "antiPD1+antiCTLA4"}
    if(regimen[i] == "IpNi"){regimen[i] = "antiCTLA4+antiPD1"}
    if(regimen[i] == "LAG3"){regimen[i] = "antiPD1+antiLAG3"}
    if(regimen[i] == "NANA"){regimen[i] = "TreatmentNA"}
    if(substr(name[i], 1, 2) %in% c("HI", "Ke")){regimen[i] = "Ctrl"}
    
  }

  for(i in 1:length(overall)){
    if(overall[i] == "U"){overall[i] = "NA"}
    if(substr(name[i], 1, 2) %in% c("HI", "Ke")){overall[i] = "Ctrl"}
    
  }

  for(i in 1:length(io_stat)){
    if(io_stat[i] == "pr"){io_stat[i] = "Prior.IO"}
    if(io_stat[i] == "nv"){io_stat[i] = "IO.naive"}
    if(substr(name[i], 1, 2) %in% c("HI", "Ke")){io_stat[i] = "Ctrl"}

  }

  return(data.frame(name, response, type, regimen, overall, timepoint, celltype, io_stat, stringsAsFactors = F))

}


UpdateFilename <- function(filename){

  ## If filename doesn't have Prior IO status

  ## All robert files all io naive and start with "GA"
  if(seqinr::s2c(filename)[1:2] == "GA"){return(filename)}

  filename <- gsub(".txt", "", filename)

  if(nchar(filename) != 24){
    # print("File already contains prior info status")
    return(filename)
  }

  ## Add prior IO-status if not found
  new_name_df <- merge(breakName(filename), meta_clinical, by.x = "name", by.y = "Study ID")

  if(new_name_df$`Previous IO` == "IO naive"){io_stat = "nv"}
  if(new_name_df$`Previous IO` == "Prior IO"){io_stat = "pr"}

  new_name <- paste0(filename, "_", io_stat)
  return(new_name)

}


extractName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub("\\_.*", "", str1)
}

extractFileName = function(str1){
  # strsplit(str1, "[_]")[[1]][1]
  sub(".*\\/", "", str1)
}





updateRiazNames <- function(filename){

  ## Update the prior IO status of the Riaz patients
  filenames_df <- read.delim("data/unselected_TCRseq/riaz_conversion.txt") %>%

    mutate(oldfiles = extractFileName(oldfiles), newfiles = extractFileName(newfiles)) %>%
    mutate(oldfiles = gsub(".txt", "", oldfiles), newfiles = gsub(".txt", "", newfiles))

  if(!filenames_df$oldfiles %in% filename){return(filename)}

  newname <- filenames_df[filenames_df$oldfiles == filename, 2] %>% as.character()
  return(newname)

}
