

## Load thresholds used for TCRhp predictions
thresholds_ver1 <- read.delim("results/tcrgp/threshold_table.txt", stringsAsFactors = F)
thresholds_ver2 <- read.delim("results/tcrgp/threshold_table.csv", stringsAsFactors = F, sep = ',')
thresholds_ver2 <- thresholds_ver2[,-6]

colnames(thresholds_ver2) = colnames(thresholds_ver1)
thresholds = rbind(thresholds_ver1, thresholds_ver2)

## Generate a function that understands the output
getPreds <- function(tcrgp_filename, thresholds = thresholds, fdr = 0.05){


  ## Process only the TCRGP-file
  require(dplyr)

  ## Select the fdr for the look-up table, thresholds
  if(fdr == 0){fdr = 1}
  if(fdr == 0.05){fdr = 2}
  if(fdr == 0.1){fdr = 3}
  if(fdr == 0.2){fdr = 4}

  ## Read in the prediction data
  tcrgp_df <- tcrgp_filename %>% data.table::fread() %>%
    mutate(cdr3aa = gsub("-", "", cdr3aa)) %>%
    mutate(clonotype_name = cdr3aa)

  tcrgp_df[tcrgp_df == "NaN"] <- NA

  ## Select only the epitopes found in this TCRGP prediction
  thresholds_temp   <- thresholds %>% dplyr::slice(which(thresholds$model %in% colnames(tcrgp_df)))
  thresholds_to_use <- thresholds_temp %>% dplyr::select(fdr) %>% t %>% as.vector()

  ## In the TCRGP-tcrgp_df, get rid of meta data
  tcrgp_df_meta <- as.data.frame(tcrgp_df)[ , which(!colnames(tcrgp_df) %in% thresholds_temp$model)]

  ## Select only the columns with the model and reorder accodingly
  tcrgp_df_temp <- as.data.frame(tcrgp_df)[ , which(colnames(tcrgp_df) %in% thresholds_temp$model)]
#  tcrgp_df_temp <- tcrgp_df_temp[ ,match(colnames(tcrgp_df_temp), thresholds_temp$model)]
  tcrgp_df_temp <- tcrgp_df_temp[ ,match(thresholds_temp$model, colnames(tcrgp_df_temp))]

  if(colnames(tcrgp_df_temp) == thresholds_temp$model){

    ## Main part: filter the cdr3s with low predictions
    antigen_names <- colnames(tcrgp_df_temp)
    tcrgp_df_temp <- apply(tcrgp_df_temp, 1, function(x){ x[x < thresholds_to_use] <- 0; return(data.frame(matrix(x,nrow=1)))}) %>% rbindlist()
    colnames(tcrgp_df_temp) <- antigen_names
    tcrgp_df_temp <- cbind(tcrgp_df_meta, tcrgp_df_temp)

    return(tcrgp_df_temp)

  }

  else{message("No matching models")}

}

## Summarise predictions
summarise_predictions  <- function(orig_file, pred_file){

  pred_file = pred_file[pred_file$cdr3aa %in% orig_file$cdr3aa, ]
  pred_file = pred_file[match(orig_file$cdr3aa, pred_file$cdr3aa), ]

  ## Check that the files match (at least on number)
  if(nrow(orig_file) != nrow(pred_file)){
    message("Error: files do not match")
    return(NULL)
  }


  ## Get the models from the threshold look-up table
  pred_file_meta <- pred_file %>% dplyr::select(which(!colnames(pred_file) %in% thresholds$model))
  pred_file      <- pred_file %>% dplyr::select(which(colnames(pred_file) %in% thresholds$model))
  pred_file[is.na(pred_file)] <- 0

  ## To find predictions, first define which cells deviate from 0
  getPredictions <- function(row){

    ## If no prediction is found:
    if(sum(row) == 0){
      pred_epitope    <- ""
      n_pred_epitopes <- 0
      return(data.frame(pred_epitope = pred_epitope, n_pred_epitopes = n_pred_epitopes))
    }

    antigen_names     = names(row)
    specific_antigens = which(row > 0)
    pred_epitope      <- paste(antigen_names[specific_antigens], collapse = ", ")
    n_pred_epitopes   <- length(specific_antigens)
    return(data.frame(pred_epitope = pred_epitope, n_pred_epitopes = n_pred_epitopes))

  }


  summ_file      <- do.call(rbind, pbapply::pbapply(pred_file, 1, FUN = getPredictions))

  ## Summarise
  summ_file$target <- ifelse(summ_file$pred_epitope %in% c(melanoma_epitopes, "multi_melanoma"), "anti-melanoma", "anti-viral")
  summ_file$target[summ_file$n_pred_epitopes > 1]  <- "Multi"
  summ_file$target[summ_file$n_pred_epitopes == 0] <- "None"

  df <- cbind(orig_file, pred_file, summ_file)
  return(df)

}

summarise_predictions2 <- function(pred_file){

  ## Only pred file

  ## Get the models from the threshold look-up table
  pred_file_meta <- pred_file %>% dplyr::select(which(!colnames(pred_file) %in% thresholds$model))
  pred_file      <- pred_file %>% dplyr::select(which(colnames(pred_file) %in% thresholds$model))
  pred_file[is.na(pred_file)] <- 0

  ## To find predictions, first define which cells deviate from 0
  getPredictions <- function(row){

    ## If no prediction is found:
    if(sum(row) == 0){
      pred_epitope    <- ""
      n_pred_epitopes <- 0
      return(data.frame(pred_epitope = pred_epitope, n_pred_epitopes = n_pred_epitopes))
    }

    antigen_names     = names(row)
    specific_antigens = which(row > 0)
    pred_epitope      <- paste(antigen_names[specific_antigens], collapse = ", ")
    n_pred_epitopes   <- length(specific_antigens)
    return(data.frame(pred_epitope = pred_epitope, n_pred_epitopes = n_pred_epitopes))

  }

  summ_file      <- do.call(rbind, pbapply::pbapply(pred_file, 1, FUN = getPredictions))
  # summ_file      <- do.call(rbind, apply(pred_file, 1, FUN = getPredictions))


  ## Summarise
  summ_file$target <- ifelse(summ_file$pred_epitope %in% melanoma_epitopes, "anti-melanoma", "anti-viral")
  summ_file$target[summ_file$n_pred_epitopes > 1]  <- "Multi"
  summ_file$target[summ_file$n_pred_epitopes == 0] <- "None"

  df <- cbind(pred_file_meta, pred_file, summ_file)
  return(df)

}



##########


# row = emerson_df[1,] %>% as.character() %>% as.vector() %>% as.character()
# 
# orig_file <- emerson_df[1,1] %>% as.character()
# pred_file <- emerson_df[1,2] %>% as.character()
# tcrgp_filename = pred_file
# df2 <- fread(pred_file)


##########


## Summarise predictions, special to TCRb
summariseTCRb <- function(row){

  ## Summarise the TCRb predictions from TCRGP
  # @ input: row of a df where first column is the original filename and second is the TCRgp result
  # @ output: summarised form of the predictions

  row <- c( "data/unselected_TCRseq/ACT//TIL/autologous/166_FFPEQ1.tsv", "results/tcrgp/raw/summary/ACT//TIL/autologous/166_FFPEQ1.txt")
  
  orig_file = row[1] %>% as.character()
  pred_file = row[2] %>% as.character()

  message(orig_file)

  df <- fread(orig_file)

  ## Get predictions
  pred     <- getPreds(tcrgp_filename = pred_file, thresholds = thresholds, fdr = 0.05)
  tcrgp_df <- summarise_predictions(orig_file = df, pred_file = pred) %>% mutate(filename =  extractFileName(orig_file))

  tcrgp_df <- cbind(breakName(tcrgp_df$filename), tcrgp_df)
  return(tcrgp_df)

}

## TCRGP produces files that are cut in multiple files. Get them in one file.
concatenatePredFiles <- function(folder){

  filenames = list.files(folder)

  for(filename in filenames){

    if(filename %in% list.files(folder)){

      # filename <- substr(filename, 39, nchar(filename))
      filename <- gsub( ".txt.*$", "", as.character(filename))

      if(nchar(filename) > 24){

        filename_short <- substr(filename, 1, 24)
        file_list <- grep(substr(filename, 1, 24), filenames, value = T)

        df <- NULL; i <- 1
        for(file in file_list){
          df[[i]] <- data.table::fread(paste0(folder, file))
          system(paste("rm", paste0(folder, file)))
          i <- i + 1
        }

        df <- do.call(rbind, df)
        write.table(df, paste0(folder, filename_short, ".txt"), sep = ",", row.names = F, quote = F)

        print(paste("Wrote", filename_short))

      }
    }
  }
}

## Translate epitopes to species
getSpecies <- function(models){

  # models = as.character(thresholds$model)
  epitopes <- gsub( "_.*$", "", as.character(models))

  species = plyr::revalue(as.factor(epitopes), c("ATDALMTGY"            = "HCV",
                                                 "CINGVCWTV"            = "HCV",
                                                 "KLVALGINAV"           = "HCV",

                                                 "EIYKRWII"             = "HIV1",
                                                 "FLKEKGGL"             = "HIV1",
                                                 "FRDYVDRFYKTLRAEQASQE" = "HIV1",
                                                 "KAFSPEVIPMF"          = "HIV1",
                                                 "GPGHKARVL"            = "HIV1",
                                                 "KRWIILGLNK"           = "HIV1",
                                                 "TPQDLNTML"            = "HIV1",

                                                 "GLCTLVAML"            = "EBV",
                                                 "RAKFKQLL"             = "EBV",
                                                 "YVLDHLIVV"            = "EBV",

                                                 "GILGFVFTL"            = "INFa",
                                                 "PKYVKQNTLKLAT"        = "INFa",

                                                 "GTSGSPIINR"           = "DENV",
                                                 "GTSGSPIVNR"           = "DENV",

                                                 "LLWNGPMAV"            = "YEF",
                                                 "melana"               = "MELANOMA",
                                                 "meloe1"               = "MELANOMA",

                                                 "NLVPMVATV"            = "CMV",
                                                 "IPSINVHHY"            = "CMV",
                                                 "TPRVTGGGAM"           = "CMV",

                                                 "RPRGEVRFL"            = "HSV2"))

  species = as.character(species)
  return(species)

}

## Correct multi-category in the TCRGP-predictions to melanoma, if some of the targets are melanoma
## (which is the most relevant in this cohort)
multiToMelanoma <- function(df){

  multi_df       <- df %>% filter(target == "Multi")
  single_df      <- df %>% filter(target != "Multi")

  melanoma_multi <- do.call("c", lapply(melanoma_epitopes, function(x) grep(x, multi_df$pred_epitope)))

  multi_df[melanoma_multi, "target"]       <- "anti-melanoma"
  multi_df[melanoma_multi, "species"]      <- "MELANOMA"
  multi_df[melanoma_multi, "pred_epitope"] <- "multi_melanoma"

  df <- rbind(single_df, multi_df)
  return(df)

}


## Merge ver1 and ver2 TCRGP result files together
summariseTCRGP = function(row){

  ## Merge ver1 and ver2 TCRGP result files together
  df1      = row[1] %>% as.character() %>% fread
  df2      = row[2] %>% as.character() %>% fread
  filename = row[3] %>% as.character()


  print(paste0('############################################################'))

  print(paste0('#####################', filename, '#####################'))

  print(paste0('############################################################'))


  ## Merge res1 and res2 files together
  df1 = df1 %>% mutate(clonotype = paste(cdr3aa, v, d, j, sep = '_')) %>% dplyr::select(clonotype, melana_cdr3b:TPQDLNTML_cdr3b, filename)
  df2 = df2 %>% mutate(clonotype = paste(cdr3aa, v, d, j, sep = '_')) %>% dplyr::select(-filename)

  df1 = df1[df1$clonotype %in% df2$clonotype, ]
  df1 = df1[match(df2$clonotype, df1$clonotype), ]

  df = cbind(df1, df2)
  df = df[,-47]
  df = df %>% dplyr::select(name:TPQDLNTML_cdr3b, ELAGIGILTV_cdr3b_comb:AMFWSVPTV_cdr3b, pred_epitope:clonotype)

  colnames(df)
  breakName(df$filename) %>% head

  ## Produce new summarisations; takes a long time so divide the cdr3s into batches of 50e3
  n       = 50e3
  nr      = nrow(df)
  df_list = split(df, rep(1:ceiling(nr/n), each = n, length.out = nr))

  pred_file_list = lapply(df_list, summarise_predictions2)
  pred_file      = do.call(rbind, pred_file_list)

  # pred_file = summarise_predictions2(pred_file = df)
  write.table(pred_file, paste0('results/tcrgp/summary/meta/', filename), sep = '\t', quote = F, row.names = F)

}





## Correct multi-category in the TCRGP-predictions to melanoma, if some of the targets are melanoma
## (which is the most relevant in this cohort)

CorrectMelanomaSpecies <- function(df){

  df$species = ifelse(df$pred_epitope %in% melanoma_epitopes, 'MELANOMA', df$species)
  return(df)

}



filterTCRGP <- function(df){
  
  
  # models = as.character(thresholds$model)
  epitopes <- gsub("_.*$", "", colnames(df))
  
  species = plyr::revalue(as.factor(epitopes), c(
                                                 
                                                 "GLCTLVAML"            = "EBV",
                                                 "RAKFKQLL"             = "EBV",
                                                 "YVLDHLIVV"            = "EBV",
                                                 
                                                 "GILGFVFTL"            = "INFa",
                                                 "PKYVKQNTLKLAT"        = "INFa",
                                                 
                                                 
                                                 "melana"               = "MELANOMA",
                                                 "meloe1"               = "MELANOMA",
                                                 "mar"                  = "MELANOMA",
                                                 "meloe1"               = "MELANOMA",
                                                 
                                                 
                                                 "NLVPMVATV"            = "CMV",
                                                 "IPSINVHHY"            = "CMV",
                                                 "TPRVTGGGAM"           = "CMV",
                                                 
                                                 "RPRGEVRFL"            = "HSV2"))
  
  species = as.character(species)
  return(species)
  
  
  ## Get rid of data what we don't know how to handle
  df <- df[df$antigen.species %in% c("CMV", "EBV", "HomoSapiens", "HSV-1", "HSV-2", "InfluenzaA"), ]
  return(df)
  
}

multiToMelanoma <- function(df){

  multi_df       <- df %>% filter(target == "Multi")
  single_df      <- df %>% filter(target != "Multi")

  melanoma_multi <- do.call("c", lapply(melanoma_epitopes, function(x) grep(x, multi_df$pred_epitope)))

  multi_df[melanoma_multi, "target"]       <- "anti-melanoma"
  multi_df[melanoma_multi, "species"]      <- "MELANOMA"
  multi_df[melanoma_multi, "pred_epitope"] <- "multi_melanoma"

  df <- rbind(single_df, multi_df)
  return(df)

}

correctColnames <- function(df){

  dupColnames = colnames(df)[duplicated(colnames(df))]

  colnames(df) = colnames(df) %>% make.unique
  df = df %>% dplyr::select(-dupColnames)
  colnames(df) = gsub('.1', '', colnames(df))
  return(df)

}



mergeTCRGP <- function(row){

  ## Merge ver1 and ver2 TCRGP result files together



  df1      = row[1] %>% as.character() %>% fread
  df2      = row[2] %>% as.character() %>% fread
  filename = row[3] %>% as.character()


  print(paste0('############################################################'))

  print(paste0('#####################', filename, '#####################'))

  print(paste0('############################################################'))


  ## Merge res1 and res2 files together
  df1 <- df1 %>% mutate(clonotype = paste(filename, cdr3aa, v, d, j, sep = '_')) %>% dplyr::select(clonotype, melana_cdr3b:TPQDLNTML_cdr3b, filename)
  df2 <- df2 %>% mutate(clonotype = paste(filename, cdr3aa, v, d, j, sep = '_')) %>% dplyr::select(-filename)

  df1 <- df1[df1$clonotype %in% df2$clonotype, ]
  df1 <- df1[match(df2$clonotype, df1$clonotype), ]

  df = cbind(df1, df2)
  df = df[,-47]
  df = df %>% dplyr::select(name:TPQDLNTML_cdr3b, ELAGIGILTV_cdr3b_comb:AMFWSVPTV_cdr3b, pred_epitope:clonotype)

  patient_info <- colnames(breakName(df$filename))
  df <- df %>% dplyr::select(-patient_info)
  df <- df %>% bind_cols(breakName(df$filename))

  ## Produce new summarisations; takes a long time so divide the cdr3s into batches of 50e3
  n       = 50e3
  nr      = nrow(df)
  df_list = split(df, rep(1:ceiling(nr/n), each = n, length.out = nr))

  pred_file_list = lapply(df_list, summarise_predictions2)
  pred_file      = do.call(rbind, pred_file_list)

  aggregate(freq ~ filename, pred_file, sum)

  # pred_file = summarise_predictions2(pred_file = df)
  write.table(pred_file, paste0('results/tcrgp/summary/meta/', filename), sep = '\t', quote = F, row.names = F)

}














plot_line_individual    <- function(tot_df){

  ggplot() +
    geom_area(data = tot_df, aes(timepoint.x,freq,fill=species,group=species), alpha = 0.3) +
    labs(x = "", y = "fraction of cells") + facet_wrap(~species,ncol=4) + facets_nice

}

summarisePatientSpecies <- function(patient, clinical_file, prediction_file, outputdir){

  dir.create(outputdir, showWarnings = F)

  patient_pred    = prediction_file %>% filter(name == patient)
  clinical_data   = clinical_file %>% filter(name == patient)
  filename        = paste0(paste(as.character(clinical_data), collapse = "_"), ".pdf")
  print(filename)

  patient_summary <- aggregate(freq ~ species + name + timepoint,
                               filter(patient_pred, target != "Multi" & species %in% interesting_species_tcgrp),
                               sum, drop = F) %>%
    tidyr::replace_na(list(freq = 0)) %>%
    left_join(clinical_file, by = "name")

  if(length(unique(patient_summary$timepoint.x)) == 1){return()}

  plot_line_individual(patient_summary)
  ggsave(paste0(outputdir, filename), width = 8, height = 3)

}

calculateFoldchange     <- function(clinical_file, prediction_file){

  patient_summary_0m <- aggregate(freq ~ species + name,
                                  filter(prediction_file, target != "Multi" & species %in% interesting_species_tcgrp & timepoint == '0m'),
                                  sum, drop = F) %>%
    tidyr::replace_na(list(freq = 0)) %>%
    left_join(clinical_file, by = "name") %>% dplyr::select(-timepoint)

  patient_summary_3m <- aggregate(freq ~ species + name,
                                  filter(prediction_file, target != "Multi" & species %in% interesting_species_tcgrp & timepoint == '3m'),
                                  sum, drop = F) %>%
    tidyr::replace_na(list(freq = 0))

  ## Consider only the patients found in both of the time points
  cmn.names          = intersect(patient_summary_0m$name, patient_summary_3m$name)
  patient_summary_0m = patient_summary_0m %>% filter(name %in% cmn.names) %>% mutate(merge_id = paste(name, species, sep = "_")) # %>% dplyr::select(-timepoint)
  patient_summary_3m = patient_summary_3m %>% filter(name %in% cmn.names) %>% mutate(merge_id = paste(name, species, sep = "_")) # %>% dplyr::select(-timepoint)

  patient_summary    = merge(patient_summary_0m, patient_summary_3m, by = "merge_id") %>% mutate(foldchange = log2(freq.y / freq.x))

  return(patient_summary)

}

plotFoldchange          <- function(foldchange_file, plotSigf = T){

  ## Calculate different comparisons
  if(plotSigf){

    regimens    = combn(unique(foldchange_file$regimen), 2) %>% t %>% as.data.frame()
    comparisons = NULL
    for(i in 1:nrow(regimens)){
      comparisons[[i]] <- regimens[i,] %>% sapply(as.character) %>% as.vector()
    }

    ## Where the significance lines should be plotted
    y_positions = c(limit-0.2, limit-1.5,limit-3)

    if(calculate_y_pos){
      nComparisons  = length(comparisons)
      y_positions   = seq(1, limit, length.out = nComparisons)
    }

  }



  limit = foldchange_file$foldchange[is.finite(foldchange_file$foldchange)] %>% abs %>% max %>% round(0)

  foldchange_file %>%
    ggplot(aes(overall,foldchange, fill = overall)) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = 3, color = 'black') +
    geom_jitter(size = 0.5) +
    facet_wrap(~species.x, ncol = 4) + facets_nice +
    ggsignif::geom_signif(comparisons = list(c('N', 'R'))) +
    response_fill + facets_nice + theme(legend.position = 'none') +
    labs(x = "", y = "log2(fc)") +
    ylim(-limit, limit)

}




plotFoldchangeComparisons <- function(foldchange_file, x_axis = "regimen", fill = x_axis, calculate_y_pos = F, map_signif_level = F, plotSigf = T, sigfVariable = "regimen"){

  ## Set upper and lower limits so that the figures are easier to interpret
  limit = foldchange_file$foldchange[is.finite(foldchange_file$foldchange)] %>% abs %>% max %>% round(0)

  ## Calculate different comparisons
  if(plotSigf){

    regimens    = combn(unique(foldchange_file$regimen), 2) %>% t %>% as.data.frame()
    comparisons = NULL
    for(i in 1:nrow(regimens)){
      comparisons[[i]] <- regimens[i,] %>% sapply(as.character) %>% as.vector()
    }

    ## Where the significance lines should be plotted
    y_positions = c(limit-0.2, limit-1.5,limit-3)

    if(calculate_y_pos){
      nComparisons  = length(comparisons)
      y_positions   = seq(1, limit, length.out = nComparisons)
    }

  } else{comparisons = NULL
  y_positions = NULL}



  ## Actual plot
  foldchange_file %>%
    ggplot(aes_string(x_axis, 'foldchange', fill = fill)) +
    geom_hline(yintercept = 0, linetype = 3, color = 'black') +

    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.5) +
    ggsignif::geom_signif(comparisons = comparisons, na.rm = T, y_position = y_positions, map_signif_level = map_signif_level) +
    facets_nice + theme(legend.position = 'none') +
    labs(x = "", y = "log2(fc)") +
    scale_fill_manual(values = brewer.pal(9, 'Pastel2')) +
    ylim(-limit, limit) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

}


correctMelanomaTarget <- function(df){
  
  ## Summarise
  df$target <- ifelse(df$pred_epitope %in% c(melanoma_epitopes, "multi_melanoma"), "anti-melanoma", df$target)
  return(df)
}
