

createPatientSummary <- function(name){

  ## Create summary files for each patient, where is stored
  
  # i)    patient and clinical variables
  # ii)   number of expanded clonotypes
  # iii)  abundance and amount of TCRGP-predicted a) anti-melanoma and b) anti-viral clonotypes; in baseline and follow-up
  # iv)   abundance and amount of extended VDJdb-matched a) anti-melanoma and b) anti-viral clonotypes; in baseline and follow-up
  # v)    abundance and amount of GLIPH clusters; in baseline and follow-up
  # vi)   number of public clonotypes and generation probabilities; in baseline and follow-up
  # vii)  diversity statistics; at baseline and follow-up
  # viii) tumor turnover when applicable
  
  Name <- paste0(toupper(substr(name,1,1)), substr(name,2,nchar(name)))
  
  ## Read in the files based on the name. The folder structures could change in the future
  cat("reading files...")
  clinical    <- fread("data/clinical/clinical_meta.csv", sep = ",") %>% dplyr::rename(name = "Study ID") #%>% filter(Cohort == Name)
  
  raw_data    <- fread(paste0("results/pgen/", name, ".tsv")) 
  raw_data    <- raw_data %>% bind_cols(breakName(raw_data$name)) %>% mutate(clonotype_id = paste0(name1, timepoint, cdr3aa))
  expanded    <- fread(paste0("results/expansion/expanded/", name, "_sigf_expanded.txt")) %>% dplyr::select(cdr3aa, p_val, BH.pval, log2_FC_count, direction, sample1_name, Sample1_count, Sample2_count)
  
  tcrgp_data  <- fread(paste0("results/tcrgp/summary/meta/", name, "_predictions.txt")) %>% correctColnames %>% mutate(species = getSpecies(pred_epitope)) %>% multiToMelanoma %>% mutate(clonotypename = paste(name, timepoint, cdr3aa, sep = "_")) %>% CorrectMelanomaSpecies() %>% 
    select(cdr3aa, pred_epitope, n_pred_epitopes, target, filename, freq) 
  
  vdjdb_data  <- fread(paste0("results/vdjdb/summary/", name, "_vdjdb.txt")) %>% preprocessMultiVDJdb %>% dplyr::select(cdr3aa, score, mhc.a, antigen.epitope:antigen.species, vdjdb.score, filename, target)  %>% mutate(antigen.gene = plyr::revalue(antigen.gene, c("MelanA" = "MLNA"))) %>%
    mutate(filename = substr(filename, 1, nchar(filename) - 6)) 
  vdjdb_data <- vdjdb_data %>% bind_cols(breakName(vdjdb_data$filename)) %>%  mutate(clonotype_id = paste0(name, timepoint, cdr3aa)) 
  
  public_data <- fread(paste0("results/public/", Name, "/", name, "_baseline.join.aa.table.txt")) #%>% dplyr::select(cdr3aa, occurences, filename) 
  pgen_data   <- read.delim(paste0("results/pgen/", name, "_pgens.tsv"), header = F, stringsAsFactors = F) %>% 
    dplyr::rename(cdr3aa = V1, p.gen = V2) %>%
    mutate(p.gen = as.numeric(p.gen)) %>% 
    mutate(p.adj.gen = p.adjust(p.gen, method = "BH")) 

  diversity   <- fread(paste0("results/diversity/calculated/", name,".diversity.aa.exact.txt")) %>% dplyr::rename(filename = "sample_id")


  ## Aggregate to sample level
  cat("aggregate files...")
  expanded_abundance_sample1 <- aggregate(Sample1_count ~ sample1_name, expanded, sum) %>% dplyr::rename(filename = sample1_name,  expanded_counts1 = Sample1_count)
  expanded_abundance_sample2 <- aggregate(Sample2_count ~ sample1_name, expanded, sum) %>% dplyr::rename(filename = sample1_name,  expanded_counts2 = Sample2_count)
  expanded_amount            <- expanded %>% group_by(sample1_name) %>% dplyr::summarise(n = n()) %>% dplyr::rename(filename = sample1_name, expanded_clonotypes = n)
  
  tcrgp_abundance            <- aggregate(freq ~ filename + pred_epitope, subset(tcrgp_data, target != "Multi"), sum)
  tcrgp_abundance            <- cast(tcrgp_abundance, filename ~ pred_epitope, value = "freq") %>% dplyr::rename(no_target = V1)
  colnames(tcrgp_abundance)[-1] <- paste0("tcrgp.adundance.", colnames(tcrgp_abundance)[-1])
  
  tcrgp_amount               <- subset(tcrgp_data, target != "Multi") %>% group_by(filename, pred_epitope) %>% dplyr::summarise(n = n())
  tcrgp_amount               <- cast(tcrgp_amount, filename ~ pred_epitope, value = "n") %>% dplyr::rename(no_target = V1)
  colnames(tcrgp_amount)[-1] <- paste0("tcrgp.amount.", colnames(tcrgp_amount)[-1])
  
  vdjdb_df                   <- merge(vdjdb_data, raw_data, by = "clonotype_id", all.x = T)
  vdjdb_abundance            <- aggregate(freq ~ filename + antigen.epitope, subset(vdjdb_df, target != "Multi"), sum)
  vdjdb_abundance            <- cast(vdjdb_abundance, filename ~ antigen.epitope, value = "freq")# %>% dplyr::rename(no_target = V1)
  colnames(vdjdb_abundance)[-1] <- paste0("vdjdb.adundance.", colnames(vdjdb_abundance)[-1])
  
  vdjdb_amount               <- subset(subset(vdjdb_data, target != "Multi"), target != "Multi") %>% group_by(filename, antigen.epitope) %>% dplyr::summarise(n = n())
  vdjdb_amount               <- cast(vdjdb_amount, filename ~ antigen.epitope, value = "n") #%>% dplyr::rename(no_target = V1)
  colnames(vdjdb_amount)[-1] <- paste0("vdjdb.amount.", colnames(vdjdb_amount)[-1])
  
  public_abundance           <- colSums(public_data[,-c(1:14)]) %>% data.frame() %>% tibble::rownames_to_column("filename") %>% dplyr::rename(baseline.public.abundance = ".")
  public_amount              <- colSums(public_data[,-c(1:14)] > 0) %>% data.frame() %>% tibble::rownames_to_column("filename") %>% dplyr::rename(baseline.public.amount = ".")
  
  rare_public                <- merge(public_data, pgen_data, by = "cdr3aa", all.x = T) 
  rare_public                <- rare_public[!duplicated(rare_public), ]
  rare_public                <- rare_public %>% filter(p.gen < 1e-7) %>% dplyr::select(-p.gen, -p.adj.gen)
  
  rare_public_abundance      <- colSums(rare_public[,-c(1:14)]) %>% data.frame() %>% tibble::rownames_to_column("filename") %>% dplyr::rename(rare.baseline.public.abundance = ".")
  rare_public_amount         <- colSums(rare_public[,-c(1:14)] > 0) %>% data.frame() %>% tibble::rownames_to_column("filename") %>% dplyr::rename(rare.baseline.public.amount = ".")
  
  ## Mege all above
  cat("merge files...")
  merged_list <- list(tcrgp_abundance, tcrgp_amount,
                      vdjdb_abundance, vdjdb_amount,
                      public_abundance, public_amount, rare_public_abundance, rare_public_amount)
  
  tcrgpd                  <- merge(tcrgp_abundance, tcrgp_amount, by = "filename")
  vdjdbd                  <- merge(vdjdb_abundance, vdjdb_amount, by = "filename")
             
  publicd                 <- merge(public_abundance, public_amount, by = "filename")
  rare_publicd            <- merge(rare_public_abundance, rare_public_amount, by = "filename")
             
  expanded                <- merge(expanded_abundance_sample1, expanded_abundance_sample2, by = "filename")
  expanded                <- merge(expanded, expanded_amount, by = "filename")
  
  tcrgpd$filename         <- sapply(tcrgpd$filename, UpdateFilename)
  vdjdbd$filename         <- sapply(vdjdbd$filename, UpdateFilename)
  publicd$filename        <- sapply(publicd$filename, UpdateFilename)
  rare_publicd$filename   <- sapply(rare_publicd$filename, UpdateFilename)
  expanded$filename       <- sapply(expanded$filename, UpdateFilename)

  
  ## Handle correctly NAs and replace with 0s where applicable
  predictions  <- merge(tcrgpd, vdjdbd, all.x = T)        %>% replace(is.na(.), 0)
  publics      <- merge(publicd, rare_publicd, all.x = T) %>% replace(is.na(.), 0)
  
  
  
  ## Merge all together
  total <- merge(predictions, publics,  by = "filename", all.x = T) 
  total$filename <- sapply(total$filename, UpdateFilename)
  
  total <- merge(total, expanded, by = "filename", all.x = T)
  total$filename <- sapply(total$filename, UpdateFilename)
  
  total <- merge(total, diversity, by = "filename")
  total$filename <- sapply(total$filename, UpdateFilename)
  
  if(name == "riaz"){total <- total %>% mutate(filename = updateRiazNames(filename))}
  
  total <- total %>% bind_cols(breakName(total$filename))
  total <- merge(total, clinical, by = "name")
  
  colnames(total) <- make.names(colnames(total), unique=TRUE)
  
  return(total)

}


total_tumeh     <- createPatientSummary(name = "tumeh")
total_riaz      <- createPatientSummary(name = "riaz")
total_yusko_mnc <- createPatientSummary(name = "yusko_mnc")

total_helsinki  <- createPatientSummary(name = "helsinki")
total_robert    <- createPatientSummary(name = "robert")
total_yusko_cd8 <- createPatientSummary(name = "yusko_cd8")


fwrite(total_tumeh,     "results/summary/tumeh.txt", sep = "\t", quote = F, row.names = F)
fwrite(total_riaz,      "results/summary/riaz.txt", sep = "\t", quote = F, row.names = F)
fwrite(total_yusko_mnc, "results/summary/yusko_mnc.txt", sep = "\t", quote = F, row.names = F)

fwrite(total_helsinki,  "results/summary/helsinki.txt", sep = "\t", quote = F, row.names = F)
fwrite(total_robert,    "results/summary/robert.txt", sep = "\t", quote = F, row.names = F)
fwrite(total_yusko_cd8, "results/summary/yusko_cd8.txt", sep = "\t", quote = F, row.names = F)

