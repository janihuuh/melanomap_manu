
## Clinical data. Clinical data doesn't change over time, so selecting rows with only one time point will do
helsinki_clin  <- breakName(list.files("data/unselected_TCRseq/Helsinki/")) %>% mutate(name = as.numeric(name)) %>% filter(timepoint == "0m" & type == "Blood")
tumeh_clin     <- breakName(list.files("data/unselected_TCRseq/Tumeh/")) %>% filter(timepoint == "0m")
yusko_mnc_clin <- breakName(list.files("data/unselected_TCRseq/Yusko/MNC/")) %>% mutate(name = as.numeric(name)) %>% filter(timepoint == "0m")
yusko_cd8_clin <- breakName(list.files("data/unselected_TCRseq/Yusko/CD8/")) %>% filter(!is.na(name)) %>% filter(timepoint == "0m") %>% filter(!is.na(overall))
robert_clin    <- breakName(list.files("data/unselected_TCRseq/Robert/")) %>% filter(timepoint == "0m")
riaz_clin      <- breakName(list.files("data/unselected_TCRseq/Riaz/")) %>% filter(timepoint == "0m")

act_clin       <- fread("data/clinical/ACT/act_meta.csv") %>% mutate(name = gsub("TIL", "", Name), overall = Response) %>% dplyr::select(name, overall)
meta_clinical  <- fread("data/clinical/clinical_meta.csv") %>% 
  dplyr::rename(response = `Best response`) %>% 
  mutate(overall = ifelse(response != "PD", "R", "N")) %>% 
  mutate(response  = factor(as.character(response), levels = c("CR", "PR", "SD", "PD"))) %>% 
  mutate(objective = ifelse(response %in% c("CR", "PR"), "R", "N")) %>% 
  mutate(riaz      = ifelse(response %in% c("CR", "PR"), "CR/PR", as.character(response))) %>% 
  mutate(riaz      = factor(as.character(riaz), levels = c("CR/PR", "SD", "PD")))

colnames(meta_clinical) <- colnames(meta_clinical) %>% make.names()



## Other clinical data
yusko_clin_full <- read.delim("data/clinical/yusko_clinical.txt", stringsAsFactors = F) %>% mutate(name = substr(Subject, 10,13)) %>%
  mutate(Neoatingen.Load = as.numeric(gsub("([0-9]+).*$", "\\1", Neoatingen.Load)),
         Mutation.Load   = as.numeric(gsub("([0-9]+).*$", "\\1", Mutation.Load)),
         OS              = as.numeric(gsub("([0-9]+).*$", "\\1", Overall.Survival..Months.)))



yusko_clin_full %>% distinct(name, .keep_all = T) %>% write.table("data/clinical/yusko_clinical_individual.csv", sep = ",", quote = F, row.names = F)

riaz_clin_full <- read.delim("data/clinical/Riaz_clinical.csv", stringsAsFactors = F, sep = ',')
riaz_clin_full <- merge(riaz_clin, riaz_clin_full, by.x = 'name', by.y = 'Patient') %>% distinct(name, .keep_all = T)

## To allow left_join
helsinki_clin  <- helsinki_clin  %>% mutate(name = as.character(name))
yusko_mnc_clin <- yusko_mnc_clin %>% mutate(name = as.character(name))
yusko_cd8_clin <- yusko_cd8_clin %>% mutate(name = as.character(name))

# helsinki_clin  <- helsinki_clin  %>% mutate(name = as.numeric(name))
# yusko_mnc_clin <- yusko_mnc_clin %>% mutate(name = as.numeric(name))
# yusko_cd8_clin <- yusko_cd8_clin %>% mutate(name = as.numeric(name))
