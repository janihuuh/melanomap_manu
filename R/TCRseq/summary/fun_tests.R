
getTestableColumns <- function(df){

  ## Remove non-numeric and constant columns

  nums  <- unlist(lapply(df, is.numeric))
  df    <- df[,nums]
  const <- names(df[, sapply(df, function(v) var(v, na.rm=TRUE)!=0)])
  df    <- df[,const]
  return(colnames(df))

}



getTestableFactors<- function(df){

  ## Get test.factors

  nums  <- unlist(lapply(df, is.numeric))
  df    <- df[,!nums]
  df    <- sapply(df, function(x) as.factor(as.character(x))) %>% as.data.frame()

  two.levels <- lapply(df, function(x) levels(x) %>% length) %>% unlist == 2
  two.levels <- two.levels[two.levels == T] %>% names

  return(two.levels)

}



dfToTestableFactors<- function(df){

  ## Remove non-numeric and constant columns

  nums  <- unlist(lapply(df, is.numeric))
  df[,!nums]    <- sapply(df[,!nums], function(x) as.factor(as.character(x))) %>% as.data.frame()

  return(df)

}



wilcox.test.list           <- function(df){

  ## Counts wilcoxon-test for all variables in the given df

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

t.test.list <- function(df, test.colum, test.factor){

  if(test.colum == "extrapolate_reads"){return()}

  df_temp <- data.frame(df[,test.colum], df[,test.factor])
  df_temp <- df_temp[!is.na(df_temp[,1]), ]
  df_temp[,2] <- droplevels(df_temp[,2])

  if(length(levels(df_temp[,2])) == 2){

    t.test(df[,test.colum] ~ df[,test.factor], var.equal = TRUE) %>% broom::tidy() %>% mutate(test = test.colum)

  }
}

wilcox.test.list <- function(df, test.colum, test.factor){

  if(test.colum == "extrapolate_reads"){return()}

  df_temp <- data.frame(df[,test.colum], df[,test.factor])
  df_temp <- df_temp[!is.na(df_temp[,1]), ]
  df_temp[,2] <- droplevels(df_temp[,2])

  if(length(levels(df_temp[,2])) == 2){

    wilcox.test(df[,test.colum] ~ df[,test.factor], var.equal = TRUE) %>% broom::tidy() %>% mutate(test = test.colum)

  }
}




testDF <- function(df){

  test.columns <- getTestableColumns(df)
  test.factors <- getTestableFactors(df)
  df           <- df %>% dfToTestableFactors

  test.list <- NULL
  i <- 1

  for(test.factor in test.factors){

    print(paste("testing:", test.factor))
    test.list[[i]] <- lapply(test.columns, t.test.list, df = df, test.factor = test.factor) %>% do.call(what = rbind) %>% mutate(test.factor = test.factor)
    i <- i + 1

  }

  test.df <- do.call(rbind, test.list)
  return(test.df)

}


modifyTestResults <- function(res_df){

  res_df %>% mutate(p.adj = p.adjust(p.value, method = "BH")) %>% arrange(p.value)

}
