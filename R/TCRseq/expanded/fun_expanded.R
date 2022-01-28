

da_analysis_aa <- function(row, folder, name = T){
  
  ## Takes in two files and counts two-sided Fisher's test for each of the clonotype.
  ## Clonotype is defined here only based on its AA-sequence, hence for each clonotype we can get many different V, D and J-genes
  ## We take only the biggest V, D and J in the output and provide the frequency for it 
  
  ## Params
  ## @input: paths of sample and the folder
  ## @output: df with Fisher p-value
  
  require(pbapply)
  require(data.table)
  require(dplyr)
  
  sample1 = row[1]
  sample2 = row[2]
  
  ## Gather sample name
  if(name){
    sample1_name <- paste(substr(sample1, 1, nchar(sample1) - 4), sep = "")
    sample2_name <- paste(substr(sample2, 1, nchar(sample2) - 4), sep = "")
  }
  
  else{
    sample1_name <- sample1
    sample2_name <- sample2
  }

  
  ## Read the sample
  sample1_tot <- fread(paste0(folder, sample1)) %>% as.data.frame() %>% dplyr::select(count, freq, cdr3aa, v, d, j)
  sample2_tot <- fread(paste0(folder, sample2)) %>% as.data.frame() %>% dplyr::select(count, freq, cdr3aa, v, d, j)
  
  ## Aggregate different clonotypes by aa 
  sample1 <- aggregate(count ~ cdr3aa, sample1_tot, FUN = sum)
  sample2 <- aggregate(count ~ cdr3aa, sample2_tot, FUN = sum)
  
  ## Consider only the cdr3aa that are found in both files
  total_df <- merge(sample1, sample2, by = "cdr3aa") %>% 
    mutate(Sample1_count =  count.x, Sample2_count =  count.y) %>%
    mutate(Sample1_neg = sum(sample1_tot$count) - Sample1_count, Sample2_neg = sum(sample2_tot$count) - Sample2_count) 
  
  fisher_df <- total_df %>% dplyr::select(Sample1_count, Sample1_neg, Sample2_count, Sample2_neg)
  
  # Set the Fisher-function, two-sided without CI:s 
  row_fisher <- function(row) {
    f <- fisher.test(matrix(row, nrow = 2), alternative = "two.sided", conf.level = 0.95)
    return(c(row, p_val = f$p.value))
  }
  
  
  # Run the code to count Fisher p-values
  df <- apply(fisher_df, 1, row_fisher) %>% t %>% as.data.frame() 
  
  df <- df %>%     
    mutate(cdr3aa        = total_df$cdr3aa,
           sigf          = ifelse(p_val < 0.05, "Sigf", "Not_sigf"),
           BH.pval       = p.adjust(p_val, method = "BH"),
           log2_FC_count = log2(df[,"Sample2_count"] / df[,"Sample1_count"]),
           count_means   = rowMeans(df[,c("Sample1_count", "Sample2_count")])) %>%
    mutate(BH.sigf = ifelse(BH.pval < 0.05, "Sigf", "Not_sigf"),
           direction = ifelse(log2_FC_count > 0, "Up", "Down"))
  
  
  
  
  
  
  ## Get VDJ info of the clones
  ## Recall: Clonotype is defined here only based on its AA-sequence, hence for each clonotype we can get many different V, D and J-genes
  get_vdj_da <- function(cdr3aa){
    
    cdr3aa_index1 <- which(as.character(sample1$cdr3aa) %in% cdr3aa)
    cdr3aa_index2 <- which(as.character(sample2$cdr3aa) %in% cdr3aa)
    
    conv_temp     <- rbind(sample1_tot[cdr3aa_index1, c("freq", "cdr3aa", "v", "d", "j")], 
                           sample2_tot[cdr3aa_index2, c("freq", "cdr3aa", "v", "d", "j")]) %>% arrange(desc(freq))
    
    v <- names(sort(table(conv_temp$v), decreasing = T)[1])
    d <- names(sort(table(conv_temp$d), decreasing = T)[1])
    j <- names(sort(table(conv_temp$j), decreasing = T)[1])
    
    conv_temp_2 <- conv_temp %>% dplyr::select(-"freq") 
    conv_temp_2 <- conv_temp_2[!duplicated(conv_temp_2), ]
    
    ## Count freqnecies
    v_freq = 1 / length(unique(conv_temp_2$v))
    d_freq = 1 / length(unique(conv_temp_2$d))
    j_freq = 1 / length(unique(conv_temp_2$j))
    
    
    ## Combine for consensus
    tot <- cbind(cdr3aa, v, d, j, v_freq, d_freq, j_freq) %>% as.data.frame()
    return(tot)
    
  }
  
  
  ## Get VDJ-info
  df_tot_vdj  <- lapply(as.character(df$cdr3aa), get_vdj_da)
  df_tot_vdj  <- do.call(rbind, df_tot_vdj)
  
  df_tot_vdj  <- cbind(df, select(df_tot_vdj, -cdr3aa)) %>% mutate(sample1_name = sample1_name, sample2_name = sample2_name) %>% 
    select(cdr3aa, p_val, sigf, BH.pval, BH.sigf, log2_FC_count, direction, count_means, Sample1_count, Sample1_neg, Sample2_count, Sample2_neg, v, d, j, v_freq, d_freq, j_freq, sample1_name, sample2_name)
  
  return(df_tot_vdj)
  
}


get_sample_matrix <- function(vdj_files, dataset){
  
  # Helper function to get the paired samples for downstream analyses
  
  ## Params
  # in:  filenames and dataset to know how many time points are there 
  # out: df containing the sample pairing
  
  # vdj_files = list.files("data/unselected_TCRseq/Yusko/MNC/", pattern = "MNC")
  # vdj_files = list.files("data/unselected_TCRseq/Yusko/CD8/", pattern = "CD8")
  
  sample_df   <- do.call(rbind, lapply(vdj_files, breakName))
  file_matrix <- table(sample_df$name, sample_df$timepoint)

  if(dataset %in% c("Tumeh", "Robert", 'Riaz')){
    
    for(i in 1:length(vdj_files)){
      
      for(j in 1:nrow(file_matrix)){
        
        if(substr(vdj_files[i], 1,4) == rownames(file_matrix)[j]){
          
          if(substr(vdj_files[i], 19, nchar(vdj_files[i]) - 11) == "0m"){
            file_matrix[j,1] <- vdj_files[i]
          }
          
          if(substr(vdj_files[i], 19, nchar(vdj_files[i]) - 11) == "3m"){
            file_matrix[j,2] <- vdj_files[i]
          }

          
        }}}
    
    file_matrix[file_matrix == 0] <- ""
    return(file_matrix)
    
  }

    
  if(dataset == "Helsinki"){
    
    for(i in 1:length(vdj_files)){
      
      for(j in 1:nrow(file_matrix)){
        
        if(substr(vdj_files[i], 1,4) == rownames(file_matrix)[j]){
  
          if(substr(vdj_files[i], 19, nchar(vdj_files[i]) - 11) == "0m"){
            file_matrix[j,1] <- vdj_files[i]
          }
          
          if(substr(vdj_files[i], 19, nchar(vdj_files[i]) - 11) == "1m"){
            file_matrix[j,2] <- vdj_files[i]
          }
          
          if(substr(vdj_files[i], 19, nchar(vdj_files[i]) - 11) == "3m"){
            file_matrix[j,3] <- vdj_files[i]
          }
          
       
          
        }}}
    
    file_matrix[file_matrix == 0] <- ""
    return(file_matrix)
    }
  
  if(dataset == "Yusko"){
    
    for(i in 1:length(vdj_files)){
      
      for(j in 1:nrow(file_matrix)){
        
        if(substr(vdj_files[i], 1,4) == rownames(file_matrix)[j]){
          
          if(substr(vdj_files[i], 19, nchar(vdj_files[i]) - 11) == "0m"){
            file_matrix[j,1] <- vdj_files[i]
          }
          
          if(substr(vdj_files[i], 19, nchar(vdj_files[i]) - 11) == "3m"){
            file_matrix[j,2] <- vdj_files[i]
          }
          
          if(substr(vdj_files[i], 19, nchar(vdj_files[i]) - 11) == "6m"){
            file_matrix[j,3] <- vdj_files[i]
          }
          
          
        }}}
    
    file_matrix[file_matrix == 0] <- ""
    return(file_matrix)
  }
  
  
}


# sample_matrix = helsinki_mat
sample_matrix_to_pairs <- function(sample_matrix){
  
  # Get the samples that are already TCR-sequeunced in the Lag3-project as pairs to do pairwise analyses
  # in: sample_matrix
  # out: df containing the samples
  
  # sample_matrix = yusko_cd8_mat
  timepoints = ncol(sample_matrix)
  pair_list <- list()
  
  
  ## If there's only 2 time points it's fairly straightforward
  if(timepoints == 2){
    sample_matrix <- sample_matrix[!sample_matrix[,1] == "", ]
    sample_matrix <- sample_matrix[!sample_matrix[,2] == "", ]
    return(sample_matrix)
    
  }

  
    
  ## If time points are more than 2
  
  for(i in 1:nrow(sample_matrix)){
    
    row <- data.frame(lapply(sample_matrix[i,], as.character), stringsAsFactors=FALSE)
    
    pair1 <- c(row[,1], row[,2])
    pair2 <- c(row[,1], row[,3])
    pair3 <- c(row[,2], row[,3])
    
    if(timepoints == 3){
      pair_df <- rbind(pair1, pair2, pair3)
      pair_df <- pair_df[complete.cases(pair_df), ]
    }
    
    if(timepoints == 4){
      
      pair4 <- c(row[,1], row[,4])
      pair5 <- c(row[,2], row[,4])
      pair5 <- c(row[,3], row[,4])
      
      pair_df <- rbind(pair1, pair2, pair3, pair4, pair5, pair5)
      pair_df <- pair_df[complete.cases(pair_df), ]
    }
    
    pair_list[[i]] <- pair_df
    
  }
  
  pair_df <- do.call(rbind, pair_list)
  colnames(pair_df) <- c("Sample1", "Sample2")
  rownames(pair_df) <- 1:nrow(pair_df)
  
  ## Remove empty pairing
  pair_df <- pair_df[(rowSums(pair_df == "") ==  0), ]
  pair_df <- pair_df[!pair_df[,1] == "1", ]
  pair_df <- pair_df[!pair_df[,2] == "1", ]
  
  return(pair_df)
  
}





