


## Read in vdjdb-info
vdjdb_new       <- read.delim("data/selected_TCRseq/vdjdb_new.txt", stringsAsFactors = F) %>% 
  filter(species == "HomoSapiens") %>% 
  mutate(antigen.species = ifelse(antigen.gene %in% melanoma_antigens_vdjdb, "Melanoma", antigen.species))

vdjdb_antigens <- vdjdb_new %>% filter(antigen.species %in% interesting_species) %>% group_by(antigen.epitope) %>% dplyr::summarise(n = n()) %>% pull(antigen.epitope)
vdjdb_antigens <- c(vdjdb_antigens, "Melanoma_multi")


## Load summaries
total_tumeh     <- fread("results/summary/tumeh.txt") %>% as.data.frame()
total_riaz      <- fread("results/summary/riaz.txt") %>% as.data.frame()
total_yusko_mnc <- fread("results/summary/total_yusko_mnc.txt") %>% as.data.frame()

total_helsinki  <- fread("results/summary/helsinki.txt") %>% as.data.frame()
total_robert    <- fread("results/summary/robert.txt")  %>% as.data.frame()
total_yusko_cd8 <- fread("results/summary/yusko_cd8.txt") %>% as.data.frame()
