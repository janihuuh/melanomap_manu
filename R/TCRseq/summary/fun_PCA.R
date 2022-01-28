

doPCA <- function(df){
  
  df <- filterSummary(df)
  
  num.columns <- lapply(df, is.numeric) %>% unlist
  num.columns <- num.columns[num.columns == T] %>% names
  
  summary_clinical <- df %>% dplyr::select(-num.columns)
  summary_numeric  <- df %>% dplyr::select(num.columns) %>% select_if(~ !any(is.na(.))) 
  
  ## Do not take reads into account
  if("reads" %in% colnames(summary_numeric)){summary_numeric <- summary_numeric %>% dplyr::select(-reads)}

  summary_pca      <- prcomp(summary_numeric, scale. = T)
  summary_pca_df   <- cbind(summary_pca$x[,1:10], summary_clinical, summary_numeric)  
  
}


filterSummary <- function(df){
  
  ## Remove NA or "" columns or constant, change char to factor
  summary <- df %>% select_if(~ !all(is.na(.))) # %>% select_if(~ !all(. == "")) #%>% select_if(~ !any(is.na(.)))
  summary <- summary[, sapply(summary, function(col) length(unique(col))) > 1]
  
  summary <- summary %>% mutate_if(is.character, as.factor)
  
  ## Remove non-useful variables
  not_prediction     <- c("filename", "Best.response", "response", "regimen") 
  not_prediction     <- not_prediction[not_prediction %in% colnames(summary)]
  
  vdjdb.vars.full    <- colnames(summary) %>% grep(pattern = "vdjdb", value = T) 
  vdjdb.variables    <- sapply(strsplit(vdjdb.vars.full, "\\."), "[", 3)
  vdjdb.vars.to.rm   <- vdjdb.vars.full[!vdjdb.variables %in% vdjdb_antigens]
  
  tcrgp.vars.full    <- colnames(summary) %>% grep(pattern = "tcrgp", value = T) 
  tcrgp.variables    <- gsub("_cdr3b", "", sapply(strsplit(tcrgp.vars.full, "\\."), "[", 3))
  tcrgp.vars.to.rm   <- tcrgp.vars.full[!tcrgp.variables %in% c(vdjdb_antigens, "no_target", "multi_melanoma", extractName(melanoma_epitopes), "ELAGIGILTV_comb")]
  
  ## Remove also "amount" and "std" variables
  std.vars           <- grep("std", colnames(summary), value = T)
  amount.vars        <- grep("amount", colnames(summary), value = T)
  
  summary_classification  <- summary %>% dplyr::select(-not_prediction, -vdjdb.vars.to.rm, -amount.vars, -std.vars, -tcrgp.vars.to.rm)
  
  print(paste("Removed:", colnames(df)[!colnames(df) %in% colnames(summary_classification)]))
  return(summary_classification)
  
}

plotPCrotation <- function(pca_object, n = 10){
  
  rotation_df <- pca_object$rotation %>% as.data.frame() %>% add_rownames(var = "feature")
  
  rotation_df[,1:7] %>% melt(id = "feature") %>% group_by(variable) %>% top_n(n, wt = abs(value)) %>% 
    ggplot(aes(feature, value, fill = ifelse(value > 0, "yes", "no")), alpha = 0.8) + geom_bar(stat = "identity") + facet_wrap(~variable, scales = "free") + coord_flip() +
    scale_fill_manual(values = c("dodgerblue", "salmon")) + labs(x = "", y = "importance") + theme(legend.position = "none")
  
}
