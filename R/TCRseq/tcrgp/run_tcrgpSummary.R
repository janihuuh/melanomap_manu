
## Durante ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

## Make one tcrgp-file
## Read in the original file (on which the predictions were made on)
orig_files <- list.files("results/tcrgp/raw/durante/", full.names = T) 
orig_df    <- fread("data/scRNAseq+TCRseq/durante_tcrb.txt")  %>% mutate(cdr3aa = cdr3)

## Summarise predictions
df             <- lapply(orig_files, FUN = function(x) x %>% fread()) %>% do.call(what = cbind)
colnames(df)   <- lapply(orig_files, FUN = function(x) x %>% extractFileName %>% gsub(pattern = "durante_", replacement = "") %>% gsub(pattern = ".csv", replacement = "")) %>% do.call(what = "c")

## Combine
df_tot <- cbind(orig_df, df) %>% mutate(cdr3aa = cdr3)
write.table(df_tot, paste0("results/tcrgp/summary/meta/durante_predictions_new.txt"), sep = "\t", quote = F, row.names = F)

## Get predictions and summarise
pred     <- getPreds(tcrgp_filename = "results/tcrgp/summary/meta/durante_predictions_new.txt", thresholds = thresholds, fdr = 0.05)
tcrgp_df <- summarise_predictions(orig_file = orig_df, pred_file = pred)
write.table(tcrgp_df, paste0("results/tcrgp/summary/meta/durante_predictions_all.txt"), sep = "\t", quote = F, row.names = F)



### ACT
list.files("data/unselected_TCRseq/ACT/", full.names = T, recursive = T)

## Create folder structure
folders_to_create <- list.dirs(tcrb_data_folder)
folder_to_create <- gsub(folders_to_create, pattern = "data/unselected_TCRseq/ACT/", replacement = "results/tcrgp/raw/summary/ACT/")

for(folder in folder_to_create){
  # message(folder)
  dir.create(folder, showWarnings = F)
}

folder="//Infusion/CD8_btlaneg_infusion"                    
name="237_CD8BTLAneg.tsv"       

folder="/Infusion/CD4_infusion/"
name="518_CD4pos.txt" 

folders <- list.dirs("data/unselected_TCRseq/ACT/", full.names = F, recursive = T)[-1]
folders <- grep("TIL", folders, invert = T, value = T)

for(folder in folders){
  
  message(folder)
  
  ## For one folder at a time
  tcrb_data_folder = "data/unselected_TCRseq/ACT/"
  tcrgp_folder = "results/tcrgp/raw/ACT/"
  
  ## First, pool samples to one file
  pred_files <- list.files(paste0(tcrgp_folder, folder), full.names = T)
  names      <- pred_files %>% extractFileName() %>% extractName() %>% unique()
  
  dir.create(paste0("results/tcrgp/raw/summary/ACT/", folder), showWarnings = F)
  
  for(name in names){
    
    name_index     <- pred_files %>% extractFileName() %>% extractName() == name
    files_to_merge <- pred_files[name_index]
    
    message(name)
    
    ## Read in the original file (on which the predictions were made on)
    orig_files <- list.files(paste0(tcrb_data_folder, folder), full.names = T)
    orig_names <- orig_files %>% extractFileName() %>% extractName() %>% gsub(pattern = ".txt", replacement = "")
    
    f <- orig_files[orig_names %in% name]
    
    if(file.exists(f) && !dir.exists(f)){
      
      orig_df     <- orig_files[orig_names %in% name] %>% fread()
      output_name <- orig_files[orig_names %in% name] %>% extractFileName()
      
      filename    <- gsub(output_name, pattern = ".tsv", replacement = "")
      filename    <- gsub(filename, pattern = ".txt", replacement = "")
      
      ## Summarise predictions
      df           <- lapply(files_to_merge, FUN = function(x) x %>% fread()) %>% do.call(what = cbind)
      colnames(df) <- lapply(files_to_merge, FUN = function(x) x %>% extractFileName %>% gsub(pattern = paste0(filename, "_"), replacement = "") %>% gsub(pattern = ".csv", replacement = "")) %>% do.call(what = "c")
      
      # x <- lapply(files_to_merge, FUN = function(x) x %>% extractFileName)  %>% gsub(pattern = ".csv", replacement = "") #%>% gsub(pattern = paste0(filename, "_"), replacement = "") %>% gsub(pattern = ".csv", replacement = ""))
      # x %>% gsub(pattern = paste0(filename, "_"), replacement = "")
      # grep(x, pattern = paste0(filename, "_"))
      
      
      ## Combine
      df <- cbind(orig_df, df)
      write.table(df, paste0("results/tcrgp/raw/summary/ACT/", folder, "/", output_name), sep = "\t", quote = F, row.names = F)
      
    }
  }
}


## TIL
folders <- list.dirs("data/unselected_TCRseq/ACT/", full.names = F, recursive = T)[-1]
folders <- grep("TIL", folders, invert = F, value = T)
# folders <- grep("Infusion", folders, invert = F, value = T)

for(folder in folders){
  
  # folder = "On_therapy/PBMC/wk0/"
  message(folder)
  
  ## For one folder at a time
  tcrb_data_folder = "data/unselected_TCRseq/ACT/"
  tcrgp_folder = "results/tcrgp/raw/ACT/"
  
  ## First, pool samples to one file
  pred_files <- list.files(paste0(tcrgp_folder, folder), full.names = T)
  names      <- list.files(paste0(tcrb_data_folder, folder), full.names = T) %>% extractFileName() %>% gsub(pattern = ".tsv", replacement = "") %>% gsub(pattern = ".txt", replacement = "") %>% unique()
  
  dir.create(paste0("results/tcrgp/raw/summary/ACT/", folder), showWarnings = F)
  
  for(name in names){
    
    name_index <- pred_files %>% extractFileName() %>% gsub(pattern = ".tsv", replacement = "") 
    name_index <- grep(name, name_index)
    files_to_merge <- pred_files[name_index]
    
    message(name)
    
    ## Read in the original file (on which the predictions were made on)
    orig_files <- list.files(paste0(tcrb_data_folder, folder), full.names = T)
    orig_names <- orig_files %>% extractFileName() %>% gsub(pattern = ".tsv", replacement = "") %>% gsub(pattern = ".txt", replacement = "")
    
    f <- orig_files[orig_names %in% name]
    
    if(file.exists(f) && !dir.exists(f)){
      
      orig_df    <- orig_files[orig_names %in% name] %>% fread()
      
      ## Summarise predictions
      df             <- lapply(files_to_merge, FUN = function(x) x %>% fread()) %>% do.call(what = cbind)
      colnames(df)   <- lapply(files_to_merge, FUN = function(x) x %>% extractFileName %>% gsub(pattern = paste0(name, "_"), replacement = "") %>% gsub(pattern = ".csv", replacement = "")) %>% do.call(what = "c")
      
      ## Combine
      df <- cbind(orig_df, df)
      write.table(df, paste0("results/tcrgp/raw/summary/ACT/", folder, "/", name, ".txt"), sep = "\t", quote = F, row.names = F)
      
    }
  }
}



## Combine the results
orig_act_files <- list.files("data/unselected_TCRseq/ACT/", full.names = T, recursive = T) %>% grep(pattern = "TIL", invert = T, value = T)
pred_act_files <- list.files("results/tcrgp/raw/summary/ACT/", full.names = T, recursive = T) %>% grep(pattern = "TIL", invert = T, value = T)

orig_names <- orig_act_files %>% extractFileName() %>% extractName() %>% gsub(pattern = ".txt", replacement = "")
pred_names <- pred_act_files %>% extractFileName() %>% extractName() %>% gsub(pattern = ".txt", replacement = "")
act_df     <- data.frame(orig_act_files, pred_act_files)

act_predictions <- pbapply::pbapply(act_df, 1, summariseTCRb) %>% rbindlist()
fwrite(act_predictions, "results/tcrgp/summary/meta/act_predictions.txt", sep = "\t", quote = F, row.names = F)



orig_act_til_files <- list.files("data/unselected_TCRseq/ACT/", full.names = T, recursive = T) %>% grep(pattern = "TIL", invert = F, value = T)
pred_act_til_files <- list.files("results/tcrgp/raw/summary/ACT/", full.names = T, recursive = T) %>% grep(pattern = "TIL", invert = F, value = T)
act_til_df <- data.frame(orig_act_til_files, pred_act_til_files)

act_til_predictions <- pbapply::pbapply(act_til_df, 1, summariseTCRb) %>% rbindlist()
fwrite(act_til_predictions, "results/tcrgp/summary/meta/act_til_predictions.txt", sep = "\t", quote = F, row.names = F)






### Emerson ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
list.files("results/tcrgp/raw/Emerson/Users/janihuuh/Dropbox/Emerson/data/vdjt/", full.names = T, pattern = "ELAGIG") %>% length()

## First, pool samples to one file
pred_emerson_files <- list.files("results/tcrgp/raw/Emerson/Users/janihuuh/Dropbox/Emerson/data/vdjt/", full.names = T)
names              <- pred_emerson_files %>% extractFileName() %>% extractName() %>% unique()
dir.create("results/tcrgp/raw/Emerson/summary/", showWarnings = F)

already_summarised_samples <- list.files("results/tcrgp/raw/Emerson/summary/") %>% extractName()
already_summarised_samples %>% length()

for(name in names){
  
  ## Don't recalculate samples
  if(!names %in% gsub(".txt", "", already_summarised_samples)){
  
    name_index     <- pred_emerson_files %>% extractFileName() %>% extractName() == name
    files_to_merge <- pred_emerson_files[name_index]
    
    ## If TCRGP hasn't finished
    if(files_to_merge %>% length() == 14){

      message(name)
      
      ## Read in the original file (on which the predictions were made on)
      orig_files <- list.files("/Users/janihuuh/Dropbox/Emerson/data/vdjt/", full.names = T)
      orig_names <- orig_files %>% extractFileName() %>% extractName() %>% gsub(pattern = ".txt", replacement = "")
      orig_df    <- orig_files[orig_names %in% name] %>% fread()
      
      ## Summarise predictions
      df             <- lapply(files_to_merge, FUN = function(x) x %>% fread()) %>% do.call(what = cbind)
      colnames(df)   <- lapply(files_to_merge, FUN = function(x) x %>% extractFileName %>% gsub(pattern = paste0(name, "_"), replacement = "") %>% gsub(pattern = ".csv", replacement = "")) %>% do.call(what = "c")
      
      ## Combine
      df <- cbind(orig_df, df)
      write.table(df, paste0("results/tcrgp/raw/Emerson/summary/", name, ".txt"), sep = "\t", quote = F, row.names = F)
      
    }
  }
}


## Combine the results
orig_emerson_files <- list.files("/Users/janihuuh/Dropbox/Emerson/data/vdjt/", full.names = T)
pred_emerson_files <- list.files("results/tcrgp/raw/Emerson/summary/", full.names = T)

orig_names <- orig_emerson_files %>% extractFileName() %>% extractName() %>% gsub(pattern = ".txt", replacement = "")
pred_names <- pred_emerson_files %>% extractFileName() %>% extractName() %>% gsub(pattern = ".txt", replacement = "")
emerson_df <- data.frame(orig_emerson_files[orig_names %in% pred_names], pred_emerson_files)

## First 50
emerson_predictions1 <- pbapply::pbapply(emerson_df[1:50,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions1, "results/tcrgp/summary/meta/emerson_predictions_1.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions2 <- pbapply::pbapply(emerson_df[51:100,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions2, "results/tcrgp/summary/meta/emerson_predictions_2.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions3 <- pbapply::pbapply(emerson_df[101:150,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions3, "results/tcrgp/summary/meta/emerson_predictions_3.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions4 <- pbapply::pbapply(emerson_df[151:200,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions4, "results/tcrgp/summary/meta/emerson_predictions_4.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions5 <- pbapply::pbapply(emerson_df[201:250,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions5, "results/tcrgp/summary/meta/emerson_predictions_5.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions6 <- pbapply::pbapply(emerson_df[251:300,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions6, "results/tcrgp/summary/meta/emerson_predictions_6.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions7 <- pbapply::pbapply(emerson_df[301:350,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions7, "results/tcrgp/summary/meta/emerson_predictions_7.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions8 <- pbapply::pbapply(emerson_df[351:400,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions8, "results/tcrgp/summary/meta/emerson_predictions_8.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions9 <- pbapply::pbapply(emerson_df[351:400,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions9, "results/tcrgp/summary/meta/emerson_predictions_9.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions10 <- pbapply::pbapply(emerson_df[451:500,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions10, "results/tcrgp/summary/meta/emerson_predictions_10.txt", sep = "\t", quote = F, row.names = F)



emerson_predictions11 <- pbapply::pbapply(emerson_df[501:551,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions11, "results/tcrgp/summary/meta/emerson_predictions_11.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions12 <- pbapply::pbapply(emerson_df[551:600,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions12, "results/tcrgp/summary/meta/emerson_predictions_12.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions13 <- pbapply::pbapply(emerson_df[601:651,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions13, "results/tcrgp/summary/meta/emerson_predictions_13.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions14 <- pbapply::pbapply(emerson_df[651:700,], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions14, "results/tcrgp/summary/meta/emerson_predictions_14.txt", sep = "\t", quote = F, row.names = F)

emerson_predictions15 <- pbapply::pbapply(emerson_df[701:nrow(emerson_df),], 1, summariseTCRb) %>% rbindlist()
fwrite(emerson_predictions15, "results/tcrgp/summary/meta/emerson_predictions_15.txt", sep = "\t", quote = F, row.names = F)



# === New

## First, pool samples to one file
pred_helsinki_files <- list.files("results/tcrgp/raw/helsinki/patient_wise/", full.names = T)
names               <- pred_helsinki_files %>% extractFileName() %>% substr(start = 1, stop = 27) %>% unique()
dir.create("results/tcrgp/raw/helsinki/summary/", showWarnings = F)

for(name in names){

  name_index     <- pred_helsinki_files %>% extractFileName() %>% substr(start = 1, stop = 27) == name
  files_to_merge <- pred_helsinki_files[name_index]
  
  ## If TCRGP hasn't finished
  if(files_to_merge %>% length() == 14){
    
    message(name)
    
    ## Read in the original file (on which the predictions were made on)
    orig_files <- list.files("data/unselected_TCRseq/Helsinki/", full.names = T)
    orig_names <- orig_files %>% extractFileName() %>% substr(start = 1, stop = 27) %>% gsub(pattern = ".txt", replacement = "")
    orig_df    <- orig_files[orig_names %in% name] %>% fread()
    
    ## Summarise predictions
    df             <- lapply(files_to_merge, FUN = function(x) x %>% fread()) %>% do.call(what = cbind)
    colnames(df)   <- lapply(files_to_merge, FUN = function(x) x %>% extractFileName %>% gsub(pattern = paste0(name, "_"), replacement = "") %>% gsub(pattern = ".csv", replacement = "")) %>% do.call(what = "c")
    
    ## Combine
    df <- cbind(orig_df, df)
    write.table(df, paste0("results/tcrgp/raw/helsinki/summary/", name, ".txt"), sep = "\t", quote = F, row.names = F)
    
  }
  
}


## Helsinki cohort
orig_helsinki_files <- list.files("data/unselected_TCRseq/Helsinki/", full.names = T)
pred_helsinki_files <- list.files("results/tcrgp/raw/helsinki/summary/", full.names = T)
helsinki_df         <- data.frame(orig_helsinki_files, pred_helsinki_files)

helsinki_predictions <- pbapply::pbapply(helsinki_df, 1, summariseTCRb) %>% rbindlist()
write.table(helsinki_predictions, "results/tcrgp/summary/helsinki_predictions_new.txt", sep = "\t", quote = F, row.names = F)





# === Ver1

## These results do not contain all available epitopes, which were produced in ver2
## Note that these results also contains non-functional cdr3aa, which were filtered in ver2

## Helsinki cohort
orig_helsinki_files <- list.files("data/unselected_TCRseq/Helsinki/", full.names = T)
pred_helsinki_files <- list.files("results/tcrgp/predictions_ver1/unselected_tcrb/Helsinki/", full.names = T)
helsinki_df         <- data.frame(orig_helsinki_files, pred_helsinki_files)

helsinki_predictions <- pbapply::pbapply(helsinki_df, 1, summariseTCRb)
helsinki_predictions <- do.call(rbind, helsinki_predictions)
write.table(helsinki_predictions, "results/tcrgp/summary/ver1/helsinki_predictions.txt", sep = "\t", quote = F, row.names = F)


## Tumeh cohort
orig_tumeh_files <- list.files("data/unselected_TCRseq/Tumeh", full.names = T)
pred_tumeh_files <- list.files("results/tcrgp/predictions_ver1/unselected_tcrb/Tumeh", full.names = T)
tumeh_df         <- data.frame(orig_tumeh_files, pred_tumeh_files) 

tumeh_predictions <- pbapply::pbapply(tumeh_df, 1, summariseTCRb)
tumeh_predictions <- do.call(rbind, tumeh_predictions)
write.table(tumeh_predictions, "results/tcrgp/summary/ver1/tumeh_predictions.txt", sep = "\t", quote = F, row.names = F)


## Yusko MNC cohort
orig_yusko_mnc_files <- list.files("data/unselected_TCRseq/Yusko/MNC/", full.names = T)
pred_yusko_mnc_files <- list.files("results/tcrgp/predictions_ver1/unselected_tcrb/yusko/MNC/", full.names = T)
yusko_mnc_df         <- data.frame(orig_yusko_mnc_files, pred_yusko_mnc_files)

yusko_mnc_predictions <- pbapply::pbapply(yusko_mnc_df, 1, summariseTCRb)
yusko_mnc_predictions <- do.call(rbind, yusko_mnc_predictions)
write.table(yusko_mnc_predictions, "results/tcrgp/summary/ver1/tyusko_mnc_predictions.txt", sep = "\t", quote = F, row.names = F)


## Yusko CD8 cohort
orig_yusko_cd8_files <- list.files("data/unselected_TCRseq/Yusko/CD8/", full.names = T)
pred_yusko_cd8_files <- list.files("results/tcrgp/predictions_ver1/unselected_tcrb/yusko/CD8/", full.names = T)
yusko_cd8_df         <- data.frame(orig_yusko_cd8_files, pred_yusko_cd8_files)

yusko_cd8_predictions <- pbapply::pbapply(yusko_cd8_df, 1, summariseTCRb)
yusko_cd8_predictions <- do.call(rbind, yusko_cd8_predictions)
write.table(yusko_cd8_predictions, "results/tcrgp/summary/ver1/yusko_cd8_predictions.txt", sep = "\t", quote = F, row.names = F)


## Riaz cohort
orig_Riaz_files <- list.files("data/unselected_TCRseq/Riaz/", full.names = T)
pred_Riaz_files <- list.files("results/tcrgp/predictions_ver1/unselected_tcrb/Riaz/", full.names = T)
Riaz_df         <- data.frame(orig_Riaz_files, pred_Riaz_files)

Riaz_predictions <- pbapply::pbapply(Riaz_df, 1, summariseTCRb)
Riaz_predictions <- do.call(rbind, Riaz_predictions)
write.table(Riaz_predictions, "results/tcrgp/summary/ver1/Riaz_predictions.txt", sep = "\t", quote = F, row.names = F)


## Robert cohort
orig_robert_files <- list.files("data/unselected_TCRseq/Robert/", full.names = T)
pred_robert_files <- list.files("results/tcrgp/predictions_ver1/unselected_tcrb/Robert/", full.names = T)
robert_df         <- data.frame(orig_robert_files, pred_robert_files)

robert_predictions <- pbapply::pbapply(robert_df, 1, summariseTCRb)
robert_predictions <- do.call(rbind, robert_predictions)
write.table(robert_predictions, "results/tcrgp/summary/ver1/robert_predictions.txt", sep = "\t", quote = F, row.names = F)















## ========= Ver2

## Adds new epitopes, e.g. Mart1 (AAGIGILTV)

## Helsinki cohort
orig_helsinki_files <- list.files("data/unselected_TCRseq/Helsinki/", full.names = T)
pred_helsinki_files <- list.files("results/tcrgp/predictions_ver2/unselected_tcrb/Helsinki/", full.names = T)
helsinki_df         <- data.frame(orig_helsinki_files, pred_helsinki_files)

helsinki_predictions <- pbapply::pbapply(helsinki_df, 1, summariseTCRb)
helsinki_predictions <- do.call(rbind, helsinki_predictions)
write.table(helsinki_predictions, "results/tcrgp/summary/ver2/helsinki_predictions.txt", sep = "\t", quote = F, row.names = F)


## Tumeh cohort
orig_tumeh_files <- list.files("data/unselected_TCRseq/Tumeh", full.names = T)
pred_tumeh_files <- list.files("results/tcrgp/predictions_ver2/unselected_tcrb/Tumeh", full.names = T)
tumeh_df         <- data.frame(orig_tumeh_files, pred_tumeh_files)

tumeh_predictions <- pbapply::pbapply(tumeh_df, 1, summariseTCRb)
tumeh_predictions <- do.call(rbind, tumeh_predictions)
write.table(tumeh_predictions, "results/tcrgp/summary/ver2/tumeh_predictions.txt", sep = "\t", quote = F, row.names = F)

## Yusko MNC cohort
orig_yusko_mnc_files <- list.files("data/unselected_TCRseq/Yusko/MNC/", full.names = T)
pred_yusko_mnc_files <- list.files("results/tcrgp/predictions_ver2/unselected_tcrb/yusko/MNC/", full.names = T)
yusko_mnc_df         <- data.frame(orig_yusko_mnc_files, pred_yusko_mnc_files)

yusko_mnc_predictions <- pbapply::pbapply(yusko_mnc_df, 1, summariseTCRb)
yusko_mnc_predictions <- do.call(rbind, yusko_mnc_predictions)
write.table(yusko_mnc_predictions, "results/tcrgp/summary/ver2/tyusko_mnc_predictions.txt", sep = "\t", quote = F, row.names = F)


## Yusko CD8 cohort
orig_yusko_cd8_files <- list.files("data/unselected_TCRseq/Yusko/CD8/", full.names = T)
pred_yusko_cd8_files <- list.files("results/tcrgp/predictions_ver2/unselected_tcrb/yusko/CD8/", full.names = T)
yusko_cd8_df         <- data.frame(orig_yusko_cd8_files, pred_yusko_cd8_files)

yusko_cd8_predictions <- pbapply::pbapply(yusko_cd8_df, 1, summariseTCRb)
yusko_cd8_predictions <- do.call(rbind, yusko_cd8_predictions)
write.table(yusko_cd8_predictions, "results/tcrgp/summary/ver2/yusko_cd8_predictions.txt", sep = "\t", quote = F, row.names = F)


## Riaz cohort
orig_Riaz_files <- list.files("data/unselected_TCRseq/Riaz/", full.names = T)
pred_Riaz_files <- list.files("results/tcrgp/predictions_ver2/unselected_tcrb/Riaz/", full.names = T)
Riaz_df         <- data.frame(orig_Riaz_files, pred_Riaz_files)

Riaz_predictions <- pbapply::pbapply(Riaz_df, 1, summariseTCRb)
Riaz_predictions <- do.call(rbind, Riaz_predictions)
write.table(Riaz_predictions, "results/tcrgp/summary/ver2/Riaz_predictions.txt", sep = "\t", quote = F, row.names = F)


## Robert cohort
orig_robert_files <- list.files("data/unselected_TCRseq/Robert/", full.names = T)
pred_robert_files <- list.files("results/tcrgp/raw/predictions_ver2/unselected_tcrb/Robert/", full.names = T)
robert_df         <- data.frame(orig_robert_files, pred_robert_files)

robert_predictions <- pbapply::pbapply(robert_df, 1, summariseTCRb) %>% do.call(what = rbind)
write.table(robert_predictions, "results/tcrgp/summary/ver2/robert_predictions.txt", sep = "\t", quote = F, row.names = F)









## ACT
input_folder   <- "data/unselected_TCRseq/ACT/"
pred_folder    <- "results/tcrgp/raw/ACT/"
results_folder <- "results/tcrgp/summary/ACT/"

folder = list.dirs(input_folder)[3]

for(folder in list.dirs(input_folder)){
  
  message(folder)
  folder_to_work <- gsub("data/unselected_TCRseq/ACT/", "", folder)
  
  orig_act_files <- list.files(paste0(input_folder, folder_to_work), full.names = T)
  pred_act_files <- list.files(paste0(pred_folder, folder_to_work), full.names = T)
  act_df         <- data.frame(orig_act_files, pred_act_files)
  
  act_predictions <- pbapply::pbapply(act_df, 1, summariseTCRb)
  act_predictions <- do.call(rbind, act_predictions)
  
  dir.create(paste0(results_folder, folder_to_work), showWarnings = F)
  write.table(act_predictions, paste0(results_folder, folder_to_work, "precidctions.txt"),  sep = "\t", quote = F, row.names = F)
  
}





## Pruessmann
pruessmann_all <- list.files("data/unselected_TCRseq/Pruessmann/", full.names = T) %>% lapply(FUN = function(x) fread(x) %>% mutate(sample_id = extractFileName(x))) %>% rbindlist
fwrite(pruessmann_all, "results/pooled/pruessmann_all.txt", sep = "\t", quote = F, row.names = F)

pred_df <- list.files("results/tcrgp/raw/pruessmann/", full.names = T) %>% lapply(FUN = function(x) fread(x)) %>% Laurae::cbindlist()
colnames(pred_df) <- list.files("results/tcrgp/raw/pruessmann/")
colnames(pred_df) <- gsub("\\.csv", "", colnames(pred_df))
colnames(pred_df) <- gsub("pruessmann_", "", colnames(pred_df))

pred_df <- cbind(pruessmann_all, pred_df)
fwrite(pred_df, "results/tcrgp/raw/pruessmann/summary.txt", sep = "\t", quote = F, row.names = F)
# pred_df <- fread("results/tcrgp/raw/pruessmann/summary.txt")

fdr=0.05
fdr="threshold1"
# fdr="threshold"

## Read in the prediction data
tcrgp_df <- pred_df %>%
  mutate(cdr3aa = gsub("-", "", cdr3aa)) %>%
  mutate(clonotype_name = cdr3aa)

tcrgp_df[tcrgp_df == "NaN"] <- NA

## Select only the epitopes found in this TCRGP prediction
thresholds_temp   <- thresholds %>% slice(which(thresholds$model %in% colnames(tcrgp_df)))
thresholds_to_use <- thresholds_temp %>% dplyr::select(fdr) %>% t %>% as.vector()

## In the TCRGP-tcrgp_df, get rid of meta data
tcrgp_df_meta <- as.data.frame(tcrgp_df)[ , which(!colnames(tcrgp_df) %in% thresholds_temp$model)]

## Select only the columns with the model and reorder accodingly
tcrgp_df_temp <- as.data.frame(tcrgp_df)[ , which(colnames(tcrgp_df) %in% thresholds_temp$model)]
tcrgp_df_temp <- tcrgp_df_temp[ ,match(colnames(tcrgp_df_temp), thresholds_temp$model)]

tcrgp_df_temp <- tcrgp_df_temp[ ,match(colnames(tcrgp_df_temp), thresholds_temp$model)] # %>% colnames()

  
## Main part: filter the cdr3s with low predictions
antigen_names <- colnames(tcrgp_df_temp)
tcrgp_df_temp <- pbapply::pbapply(tcrgp_df_temp, 1, function(x){ x[x < thresholds_to_use] <- 0; return(data.frame(matrix(x,nrow=1)))})
tcrgp_df_temp <- do.call(rbind, tcrgp_df_temp)
colnames(tcrgp_df_temp) <- antigen_names
# tcrgp_df_temp[tcrgp_df_temp < thresholds_to_use] <- 0
tcrgp_df_temp <- cbind(tcrgp_df_meta, tcrgp_df_temp)

pred = tcrgp_df_temp

tcrgp_df <- summarise_predictions(orig_file = pruessmann_all, pred_file = pred) 
fwrite(tcrgp_df, "results/tcrgp/summary/meta/pruessmann_predictions.txt", sep = "\t", quote = F, row.names = F)
# fwrite(tcrgp_df, "results/tcrgp/summary/meta/pruessmann_predictions_fpr000.txt", sep = "\t", quote = F, row.names = F)








## Helsinki

## Make one tcrgp-file
## Read in the original file (on which the predictions were made on)
orig_files <- list.files("results/tcrgp/raw/helsinki/", full.names = T)
orig_df    <- fread("results/pooled/helsinki.join.strict.table.txt")

## Summarise predictions
df             <- lapply(orig_files, FUN = function(x) x %>% fread()) %>% do.call(what = cbind)
colnames(df)   <- lapply(orig_files, FUN = function(x) x %>% extractFileName %>% gsub(pattern = paste0(name, "_"), replacement = "") %>% gsub(pattern = ".csv", replacement = "")) %>% do.call(what = "c")

## Combine
df_tot <- cbind(orig_df, df)
write.table(df_tot, paste0("results/tcrgp/summary/meta/helsinki_predictions_new.txt"), sep = "\t", quote = F, row.names = F)


pred_df = df_tot
fdr=0.05

## Read in the prediction data
tcrgp_df <- pred_df %>%
  mutate(cdr3aa = gsub("-", "", cdr3aa)) %>%
  mutate(clonotype_name = cdr3aa)

tcrgp_df[tcrgp_df == "NaN"] <- NA

## Select only the epitopes found in this TCRGP prediction
thresholds_temp   <- thresholds %>% slice(which(thresholds$model %in% colnames(tcrgp_df)))
thresholds_to_use <- thresholds_temp %>% dplyr::select(fdr) %>% t %>% as.vector()

## In the TCRGP-tcrgp_df, get rid of meta data
tcrgp_df_meta <- as.data.frame(tcrgp_df)[ , which(!colnames(tcrgp_df) %in% thresholds_temp$model)]

## Select only the columns with the model and reorder accodingly
tcrgp_df_temp <- as.data.frame(tcrgp_df)[ , which(colnames(tcrgp_df) %in% thresholds_temp$model)]
tcrgp_df_temp <- tcrgp_df_temp[ ,match(colnames(tcrgp_df_temp), thresholds_temp$model)]



## Main part: filter the cdr3s with low predictions
antigen_names <- colnames(tcrgp_df_temp)
tcrgp_df_temp <- pbapply::pbapply(tcrgp_df_temp, 1, function(x){ x[x < thresholds_to_use] <- 0; return(data.frame(matrix(x,nrow=1)))})
tcrgp_df_temp <- do.call(rbind, tcrgp_df_temp)
colnames(tcrgp_df_temp) <- antigen_names
# tcrgp_df_temp[tcrgp_df_temp < thresholds_to_use] <- 0
tcrgp_df_temp <- cbind(tcrgp_df_meta, tcrgp_df_temp)

pred = tcrgp_df_temp

tcrgp_df <- summarise_predictions(orig_file = pruessmann_all, pred_file = pred) 
fwrite(tcrgp_df, "results/tcrgp/summary/meta/pruessmann_predictions.txt", sep = "\t", quote = F, row.names = F)