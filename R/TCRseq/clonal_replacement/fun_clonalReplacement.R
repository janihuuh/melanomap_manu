

## Calculate own tumour turnover
getStayAtHomeClones = function(row, folder){
  
  ## Outputs a summarised df
  file1 = row[1]
  file2 = row[2]
  
  df1 = paste0(folder, file1) %>% fread %>% mutate(id = paste0(cdr3aa, v, j))
  df2 = paste0(folder, file2) %>% fread %>% mutate(id = paste0(cdr3aa, v, j))
  
  ## Clonotypes found in df2 and df2
  df1 <- df1[df1$cdr3aa %in% df2$cdr3aa, ] %>% dplyr::select(count,freq,cdr3aa) %>% dplyr::rename(count.pre  = count, freq.pre  = freq)
  df2 <- df2[df2$cdr3aa %in% df1$cdr3aa, ] %>% dplyr::select(count,freq,cdr3aa) %>% dplyr::rename(count.post = count, freq.post = freq)
   
  df1 <- df1 %>% group_by(cdr3aa) %>% summarise(count.pre  = sum(count.pre),  freq.pre  = sum(freq.pre))
  df2 <- df2 %>% group_by(cdr3aa) %>% summarise(count.post = sum(count.post), freq.post = sum(freq.post))
  
  df1 %>% left_join(df2) %>% mutate(name = extractName(file1)) %>% return()

  }

getClonalReplacementClones = function(row, folder){
  
  ## Outputs a summarised df
  file1 = row[1]
  file2 = row[2]
  
  df1 = paste0(folder, file1) %>% fread %>% mutate(id = paste0(cdr3aa, v, j))
  df2 = paste0(folder, file2) %>% fread %>% mutate(id = paste0(cdr3aa, v, j))
  
  ## Clonotypes found in df2, not found in df1
  df2 <- df2[!df2$cdr3aa %in% df1$cdr3aa, ] %>% mutate(name = extractName(file1))
  return(df2)
  
}



## Calculate own tumour turnover
calculateOverlap = function(row, folder){

  ## Outputs a summarised df

  file1 = row[1]
  file2 = row[2]

  df1 = paste0(folder, file1) %>% fread %>% mutate(id = paste0(cdr3aa, v, j))
  df2 = paste0(folder, file2) %>% fread %>% mutate(id = paste0(cdr3aa, v, j))

  ## Because MNC files don't always contain v or j or the nt are weird, let's use aa's instead
  df1 = aggregate(freq ~ cdr3aa, df1, sum)
  df2 = aggregate(freq ~ cdr3aa, df2, sum)

  clonotypes1 = nrow(df1)
  clonotypes2 = nrow(df2)

  overlap = merge(df1, df2, by = "cdr3aa")

  a = data.frame(baseline            = sum(overlap$freq.x),
                 followup            = sum(overlap$freq.y),
                 baseline_clonotypes = clonotypes1,
                 followup_clonotypes = clonotypes2,
                 filename            = gsub('.txt', '', file1))
  return(a)

}


## Calculate own tumour turnover
calculateOverlapIndividual = function(row, folder){

  file1 = row[1]
  file2 = row[2]

  df1 = paste0(folder, file1) %>% fread %>% mutate(id = paste0(cdr3nt, v, j))
  df2 = paste0(folder, file2) %>% fread %>% mutate(id = paste0(cdr3nt, v, j))

  ## Because MNC files don't always contain v or j or the nt are weird, let's use aa's instead
  df1 = aggregate(freq ~ cdr3aa, df1, sum)
  df2 = aggregate(freq ~ cdr3aa, df2, sum)

  clonotypes1 = nrow(df1)
  clonotypes2 = nrow(df2)

  overlap = merge(df1, df2, by = "cdr3aa")
  return(overlap)

}







plotOverlapBox<- function(df){

  df %>%

    ggplot(aes(overall, 1 - followup, color = overall, fill = overall)) + geom_boxplot(outlier.shape = NA) +

    stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
    geom_jitter(size = 0.5, color = "black") +
    ggsignif::geom_signif(comparisons = list(c("N", "R")), color = "black") +
    labs(x = "", y = "new cells after therapy") +
    response_col +
    response_fill +
    theme(legend.position = "none")

}


plotExpandedBox<- function(df){

  df %>%

    ggplot(aes(overall, followup, color = overall, fill = overall)) + geom_boxplot(outlier.shape = NA) +

    stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
    geom_jitter(size = 0.5, color = "black") +
    ggsignif::geom_signif(comparisons = list(c("N", "R")), color = "black") +
    labs(x = "", y = "expanded cells after therapy") +
    response_col +
    response_fill +
    theme(legend.position = "none")

}
