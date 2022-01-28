
source("src/jani/R/tcrb/expanded/fun_expanded.R")

## Analyse the expanded clonotypes

# Get sample matrix and sample matrix into sampling pairs
helsinki_mat    <- get_sample_matrix(vdj_files = list.files("data/unselected_TCRseq/Helsinki/", pattern = "MNC"), dataset = "Helsinki")
helsinki_pairs  <- sample_matrix_to_pairs(helsinki_mat)
helsinki_exp    <- pbapply::pbapply(helsinki_pairs, 1, da_analysis_aa, folder = "data/unselected_TCRseq/Helsinki/")
helsinki_exp    <- do.call(rbind, helsinki_exp)

write.table(helsinki_exp, paste0(local_folder," results/expansion/helsinki_expansion.txt"), sep = "\t", quote = F, row.names = F)
write.table(helsinki_pairs, "results/expansion/helsinki_sample_pairs.txt", sep = "\t", quote = F, row.names = F)
write.table(helsinki_mat,   "results/expansion/helsinki_sample_matrix.txt", sep = "\t", quote = F, row.names = F)


tumeh_mat       <- get_sample_matrix(list.files("data/unselected_TCRseq/Tumeh/", pattern = "MNC"), dataset = "Tumeh")
tumeh_pairs     <- sample_matrix_to_pairs(tumeh_mat)
tumeh_exp       <- pbapply::pbapply(tumeh_pairs, 1, da_analysis_aa, folder = "data/unselected_TCRseq/Tumeh/")
tumeh_exp       <- do.call(rbind, tumeh_exp)

write.table(tumeh_exp, paste0("results/expansion/tumeh_expansion.txt"), sep = "\t", quote = F, row.names = F)
write.table(tumeh_pairs, "results/expansion/tumeh_sample_pairs.txt", sep = "\t", quote = F, row.names = F)

yusko_mnc_mat       <- get_sample_matrix(list.files("data/unselected_TCRseq/Yusko/MNC/", pattern = "MNC"), dataset = "Yusko")
yusko_mnc_pairs     <- sample_matrix_to_pairs(yusko_mnc_mat)
yusko_mnc_exp       <- pbapply::pbapply(yusko_mnc_pairs, 1, da_analysis_aa, folder = "data/unselected_TCRseq/Yusko/MNC//")
yusko_mnc_exp       <- do.call(rbind, yusko_mnc_exp)

write.table(yusko_mnc_exp, paste0("results/expansion/yusko_mnc_expansion.txt"), sep = "\t", quote = F, row.names = F)
write.table(yusko_mnc_pairs, "results/expansion/yusko_mnc_sample_pairs.txt", sep = "\t", quote = F, row.names = F)
write.table(yusko_mnc_mat,   "results/expansion/yusko_mnc_sample_matrix.txt", sep = "\t", quote = F, row.names = F)


yusko_cd8_mat       <- get_sample_matrix(list.files("data/unselected_TCRseq/Yusko/CD8/", pattern = "CD8"), dataset = "Yusko")
yusko_cd8_pairs     <- sample_matrix_to_pairs(yusko_cd8_mat)
yusko_cd8_exp       <- pbapply::pbapply(yusko_cd8_pairs, 1, da_analysis_aa, folder = "data/unselected_TCRseq/Yusko/CD8/")
yusko_cd8_exp       <- do.call(rbind, yusko_cd8_exp)


write.table(yusko_cd8_exp, paste0("results/expansion/yusko_cd8_expansion.txt"), sep = "\t", quote = F, row.names = F)
write.table(yusko_cd8_pairs, "results/expansion/yusko_cd8_sample_pairs.txt", sep = "\t", quote = F, row.names = F)
write.table(yusko_cd8_mat,   "results/expansion/yusko_cd8_sample_matrix.txt", sep = "\t", quote = F, row.names = F)


robert_mat       <- get_sample_matrix(list.files("data/unselected_TCRseq/Robert/"), dataset = "Robert")
robert_pairs     <- sample_matrix_to_pairs(robert_mat)
robert_exp       <- pbapply::pbapply(robert_pairs, 1, da_analysis_aa, folder = "data/unselected_TCRseq/robert/")
robert_exp       <- do.call(rbind, robert_exp)

write.table(robert_exp, paste0("results/expansion/robert_expansion.txt"), sep = "\t", quote = F, row.names = F)
write.table(robert_pairs, "results/expansion/robert_sample_pairs.txt", sep = "\t", quote = F, row.names = F)



riaz_mat       <- get_sample_matrix(list.files("data/unselected_TCRseq/Riaz/"), dataset = "Riaz")
riaz_pairs     <- sample_matrix_to_pairs(riaz_mat)
riaz_exp       <- pbapply::pbapply(riaz_pairs, 1, da_analysis_aa, folder = "data/unselected_TCRseq/Riaz/")
riaz_exp       <- do.call(rbind, riaz_exp)

write.table(riaz_exp, paste0("results/expansion/riaz_expansion.txt"), sep = "\t", quote = F, row.names = F)
write.table(riaz_pairs, "results/expansion/riaz_sample_pairs.txt", sep = "\t", quote = F, row.names = F)



helsinki_exp    <- helsinki_exp   %>% filter(direction == "Up" & BH.sigf == "Sigf" & log2_FC_count > 1)
tumeh_exp       <- tumeh_exp      %>% filter(direction == "Up" & BH.sigf == "Sigf" & log2_FC_count > 1)
yusko_mnc_exp   <- yusko_mnc_exp  %>% filter(direction == "Up" & BH.sigf == "Sigf" & log2_FC_count > 1)
yusko_cd8_exp   <- yusko_cd8_exp  %>% filter(direction == "Up" & BH.sigf == "Sigf" & log2_FC_count > 1)
robert_expanded <- robert_exp     %>% filter(direction == "Up" & BH.sigf == "Sigf" & log2_FC_count > 1)
riaz_expanded   <- riaz_exp       %>% filter(direction == "Up" & BH.sigf == "Sigf" & log2_FC_count > 1)


helsinki_expanded       <- helsinki_expanded  %>% bind_cols(breakName(helsinki_expanded$sample1_name))  %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))
tumeh_expanded          <- tumeh_expanded     %>% bind_cols(breakName(tumeh_expanded$sample1_name))     %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))
yusko_mnc_expanded      <- yusko_mnc_expanded %>% bind_cols(breakName(yusko_mnc_expanded$sample1_name)) %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))
yusko_cd8_expanded      <- yusko_cd8_expanded %>% bind_cols(breakName(yusko_cd8_expanded$sample1_name)) %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))
robert_expanded         <- robert_expanded    %>% bind_cols(breakName(robert_expanded$sample1_name))    %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))
riaz_expanded           <- riaz_expanded      %>% bind_cols(breakName(riaz_expanded$sample1_name))    %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))

write.table(helsinki_expanded,  "results/expansion/expanded/helsinki_sigf_expanded.txt", sep = "\t", quote = F, row.names = F)
write.table(tumeh_expanded,     "results/expansion/expanded/tumeh_sigf_expanded.txt", sep = "\t", quote = F, row.names = F)
write.table(yusko_mnc_expanded, "results/expansion/expanded/yusko_mnc_sigf_expanded.txt", sep = "\t", quote = F, row.names = F)
write.table(yusko_cd8_expanded, "results/expansion/expanded/yusko_cd8_sigf_expanded.txt", sep = "\t", quote = F, row.names = F)
write.table(robert_expanded,    "results/expansion/expanded/robert_sigf_expanded.txt", sep = "\t", quote = F, row.names = F)
write.table(riaz_expanded,      "results/expansion/expanded/riaz_sigf_expanded.txt", sep = "\t", quote = F, row.names = F)
