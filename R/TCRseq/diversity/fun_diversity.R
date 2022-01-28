
## === Baseline

fancy_boxplot              <- function(df, variable){

  df %>%
    ggplot(aes_string('overall', variable, color = 'overall', fill = 'overall', label = 'name')) +

    geom_boxplot(outlier.shape = NA) +
    stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
    geom_jitter(size = 0.5, color = 'black') +
    response_fill +
    response_col +
    theme(legend.position = 'none') +
    ggsignif::geom_signif(comparisons = list(c('N', 'R')), test = 'wilcox.test', color = 'black')
}

wilcox.test.list           <- function(df){

  wilcox.test_fun = function(df, variable){

    if(variable == 'overall'){return()}

    df = df %>% dplyr::select(variable, overall)
    colnames(df) = c('variable', 'overall')

    mean_r = df %>% filter(overall == 'R') %>% dplyr::select(variable) %>% colMeans(na.rm = T)
    mean_n = df %>% filter(overall == 'N') %>% dplyr::select(variable) %>% colMeans(na.rm = T)
    a = wilcox.test(variable ~ overall, df, na.action = "na.omit") %>% broom::tidy() %>% mutate(name = variable, mean_r, mean_n, up = ifelse(mean_r > mean_n, 'R', 'N'))
    return(a)

  }

  wilcox_df = lapply(colnames(df), wilcox.test_fun, df = df)
  wilcox_df = do.call(rbind, wilcox_df)
  return(wilcox_df)

}

corr.test.list             <- function(df){

  if(variable != 'overall'){

    df   = df %>% dplyr::select(-overall) %>% as.matrix() %>% rcorr
    df   = melt(df$P) %>% filter(Var1 != Var2) %>% filter(row_number() %% 2 != 0) %>%
      mutate(p.adj = p.adjust(value, method = 'BH')) %>% filter(value < 0.05) %>% arrange(value)

    # df_r = try(subset(cd8_infusion_diversity, overall == 'R') %>% dplyr::select(-overall) %>% as.matrix() %>% rcorr, silent = T)
    #
    # if(class(df_r) != "try-error"){
    #   df_r = melt(df_r$P) %>% filter(Var1 != Var2) %>% filter(row_number() %% 2 != 0) %>%
    #     mutate(p.adj = p.adjust(value, method = 'BH')) %>% filter(value < 0.05) %>% arrange(value)
    # }
    # else(df_r = NA)
    #
    #
    # df_N = try(subset(cd8_infusion_diversity, overall == 'N') %>% dplyr::select(-overall) %>% as.matrix() %>% rcorr, silent = T)
    #
    # if(class(df_N) != "try-error"){
    #   df_N = melt(df_N$P)%>% filter(Var1 != Var2) %>% filter(row_number() %% 2 != 0) %>%
    #     mutate(p.adj = p.adjust(value, method = 'BH')) %>% filter(p.adj < 0.05) %>% arrange(value)
    # }
    # else(df_N = NA)

    return(df)

  }


}

preprocess_for_wilcox.test <- function(df, clin_df){

  ## Remove unnecessary variables
  std.vars   = grep("std", colnames(df), value = T)
  reads.vars = grep("extrapolate_reads", colnames(df), value = T)

  df = df %>%
    mutate(name = extractName(sample_id)) %>%
    left_join(clin_df, by = 'name')

  if('overall.x' %in% colnames(df)){df <- df %>% dplyr::rename(overall = overall.x)}

  df = df %>%

    select(-std.vars, -reads.vars) %>%

    # select_if(~ !any(is.na(.))) %>%
    select_if(~ !all(. == 0)) %>%
    select_if(~ !any(is.character(.))) %>%
    select_if(~ !any(is.factor(.))) %>%
    mutate(overall = df$overall)

  return(df)

}

polishDiversityResults     <- function(df){

  df %>%
    mutate(p.adj = p.adjust(p.value, method = 'BH')) %>% arrange(p.value) %>%
      select(p.value, p.adj, name, mean_r, mean_n, up) %>%
      mutate(p.value = round(p.value, 3), p.adj = round(p.adj, 3),   mean_r = round(mean_r, 3),
             mean_n = round(mean_n, 3))


}



## === Fold-change

preprocessFoldchange       <- function(df1, df2, clin_df){

  std.vars   = grep("std", colnames(df1), value = T)

  df1_melt = df1 %>% dplyr::select(reads:name) %>% dplyr::select(-std.vars) %>% melt(id = "name") %>% mutate(merge_id = paste(name, variable, sep = "_"))
  df2_melt = df2 %>% dplyr::select(reads:name) %>% dplyr::select(-std.vars) %>% melt(id = "name") %>% mutate(merge_id = paste(name, variable, sep = "_"))

  df_melt = merge(df1_melt, select(df2_melt, -name, -variable), by = "merge_id") %>%
    mutate(foldchange = log2(value.y / value.x)) %>%
    left_join(clin_df, by = "name")

  return(df_melt)

}
