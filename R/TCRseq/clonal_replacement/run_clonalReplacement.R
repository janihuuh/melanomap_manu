
## Tumeh
tumeh_turnover <- fread("results/expansion/tumeh_sample_pairs.txt") %>% apply(1, FUN = calculateOverlap, folder = "data/unselected_TCRseq/Tumeh/") %>% rbindlist()
tumeh_turnover <- tumeh_turnover %>% bind_cols(breakName(tumeh_turnover$filename))

tumeh_turnover %>% plotOverlapBox 
ggsave("results/turnover/plots/box_tumeh.pdf", width = 6, height = 4)

## Yusko MNC
yusko_mnc_turnover <- fread("results/expansion/yusko_mnc_sample_pairs.txt") %>% apply(1, FUN = calculateOverlap, folder = "data/unselected_TCRseq/Yusko/MNC/") %>% rbindlist()
yusko_mnc_turnover <- yusko_mnc_turnover %>% bind_cols(breakName(yusko_mnc_turnover$filename))

yusko_mnc_turnover %>% filter(regimen == "antiCTLA4+antiPD1") %>% plotOverlapBox
ggsave("results/turnover/plots/box_yusko_mnc_ipi_nivo.pdf", width = 6, height = 4)

yusko_mnc_turnover %>% filter(regimen == "antiPD1+antiCTLA4") %>% plotOverlapBox
ggsave("results/turnover/plots/box_yusko_mnc_nivo_ipi.pdf", width = 6, height = 4)


## Riaz
riaz_turnover <- fread("results/expansion/riaz_sample_pairs.txt") %>% apply(1, FUN = calculateOverlap, folder = "data/unselected_TCRseq/Riaz/") %>% rbindlist()
riaz_turnover <- riaz_turnover %>% bind_cols(breakName(riaz_turnover$filename))

riaz_turnover %>% filter(overall != "NA") %>% plotOverlapBox
ggsave("results/turnover/plots/box_riaz.pdf", width = 6, height = 4)



## Plot meta 
cmn.colnanes = intersect(intersect(colnames(yusko_mnc_turnover), colnames(tumeh_turnover)), colnames(riaz_turnover))

meta_df <- rbind(yusko_mnc_turnover %>% dplyr::select(cmn.colnanes),
                 tumeh_turnover %>% dplyr::select(cmn.colnanes),
                 riaz_turnover %>% dplyr::select(cmn.colnanes)) %>% filter(overall != "NA")

meta_df$log2fc_clonotypes <- log2(meta_df$followup_clonotypes / meta_df$baseline_clonotypes)
meta_df$objective <- ifelse(meta_df$response %in% c("CR", "PR"), "R", "N")
meta_df$riaz      <- ifelse(meta_df$response %in% c("CR", "PR"), "CR/PR", meta_df$response)
meta_df$riaz <- factor(meta_df$riaz, levels = c("CR/PR", "SD", "PD"))
meta_df$clonal_replacement <- 1 - meta_df$followup
meta_df$clonal_expansion   <- log2(meta_df$followup / meta_df$baseline)



## Clonal replacement (i.e., how big of freq clonotypes are novel after therapy)
meta_df %>% ggplot(aes(riaz, clonal_replacement, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("CR/PR", "PD"))) + facet_wrap(io_stat~.,ncol=4)
# meta_df %>% ggplot(aes(riaz, clonal_replacement, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggpubr::stat_compare_means(ref = "PD", label = "p.format") + facet_wrap(io_stat~.,ncol=4)
meta_df %>% ggplot(aes(objective, clonal_replacement, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(io_stat~.,ncol=4)
meta_df %>% ggplot(aes(overall, clonal_replacement, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(io_stat~.,ncol=4)

meta_df %>% ggplot(aes(riaz, clonal_replacement, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("CR/PR", "PD"), c("CR/PR", "SD"), c("SD", "PD")), step_increase = 0.05) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette3(4)) +facet_wrap(~io_stat)

meta_df %>% ggplot(aes(io_stat, clonal_replacement, fill = io_stat)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("IO.naive", "Prior.IO")), step_increase = 0.05) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette3(4)) 


meta_df %>% ggplot(aes(objective, clonal_replacement, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)
meta_df %>% ggplot(aes(overall, clonal_replacement, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(io_stat~.,ncol=4)


meta_df %>% filter(baseline_clonotypes > 500 & followup_clonotypes > 500) %>% ggplot(aes(riaz, clonal_replacement, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("CR/PR", "PD"), c("CR/PR", "SD"), c("SD", "PD")), step_increase = 0.05) + facet_wrap(io_stat~.,ncol=4)
meta_df %>% filter(baseline_clonotypes > 500) %>% ggplot(aes(objective, clonal_replacement, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)
meta_df %>% filter(baseline_clonotypes > 500) %>% ggplot(aes(overall, clonal_replacement, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)

meta_df %>% filter(baseline_clonotypes > 1000 & followup_clonotypes > 1000) %>% ggplot(aes(riaz, clonal_replacement, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("CR/PR", "PD"), c("CR/PR", "SD"), c("SD", "PD")), step_increase = 0.05) + facet_wrap(io_stat~.,ncol=4)
meta_df %>% filter(baseline_clonotypes > 1000) %>% ggplot(aes(objective, clonal_replacement, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)
meta_df %>% filter(baseline_clonotypes > 1000) %>% ggplot(aes(overall, clonal_replacement, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)


## Clonal replacement (i.e., how much the clonotypes that were in the TME expanded following PD-1)
meta_df %>% ggplot(aes(riaz, clonal_expansion, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(io_stat~.,ncol=4)
meta_df %>% ggplot(aes(objective, clonal_expansion, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(io_stat~.,ncol=4)
meta_df %>% ggplot(aes(overall, clonal_expansion, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(io_stat~.,ncol=4)

meta_df %>% ggplot(aes(riaz, clonal_expansion, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("CR/PR", "PD"), c("CR/PR", "SD"), c("SD", "PD")), step_increase = 0.05)# + facet_wrap(io_stat~.,ncol=4)
meta_df %>% ggplot(aes(objective, clonal_expansion, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)
meta_df %>% ggplot(aes(overall, clonal_expansion, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)

meta_df %>% filter(baseline_clonotypes > 500) %>% ggplot(aes(riaz, clonal_expansion, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("CR/PR", "PD"), c("CR/PR", "SD"), c("SD", "PD")), step_increase = 0.05) + facet_wrap(io_stat~.,ncol=4)
meta_df %>% filter(baseline_clonotypes > 500) %>% ggplot(aes(objective, clonal_expansion, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)
meta_df %>% filter(baseline_clonotypes > 500) %>% ggplot(aes(overall, clonal_expansion, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)

meta_df %>% filter(baseline_clonotypes > 1000) %>% ggplot(aes(riaz, clonal_expansion, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("CR/PR", "PD"), c("CR/PR", "SD"), c("SD", "PD")), step_increase = 0.05)# + facet_wrap(io_stat~.,ncol=4)
meta_df %>% filter(baseline_clonotypes > 1000) %>% ggplot(aes(objective, clonal_expansion, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)
meta_df %>% filter(baseline_clonotypes > 1000) %>% ggplot(aes(overall, clonal_expansion, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R")))# + facet_wrap(io_stat~.,ncol=4)


rbind(yusko_mnc_turnover, tumeh_turnover, riaz_turnover[,cmn.colnanes]) %>% 
  filter(overall != "NA") %>% 
  plotOverlapBox 
ggsave("results/turnover/plots/total_mnc.pdf", width = 4, height = 4)


rbind(tumeh_turnover, riaz_turnover[,cmn.colnanes]) %>% 
  ggplot(aes(io_stat, 1 - followup, fill = io_stat)) + 
  
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  
  geom_jitter(size = 0.5, color = "black") + 
  ggsignif::geom_signif(comparisons = list(c("IO.naive", "Prior.IO")), color = "black") +
  labs(x = "", y = "new cells after therapy") +
  scale_color_manual(values = brewer.pal(2, 'Set1')) + 
  scale_fill_manual(values = brewer.pal(2, 'Set1')) + 
  theme(legend.position = "none")

ggsave("results/turnover/plots/io_stat.pdf", width = 4, height = 4)

  








## Plot meta from expansion point-of-view
cmn.colnanes = intersect(intersect(colnames(yusko_mnc_turnover), colnames(tumeh_turnover)), colnames(riaz_turnover))

rbind(yusko_mnc_turnover[,cmn.colnanes], tumeh_turnover[,cmn.colnanes], riaz_turnover[,cmn.colnanes]) %>% 
  filter(overall != "NA") %>% 
  plotExpandedBox() + facet_wrap(io_stat~regimen,ncol=4) + facets_nice
ggsave("results/expansion/plots/expansion_space/total_mnc_regimen.pdf", width = 14, height = 4)

rbind(yusko_mnc_turnover, tumeh_turnover, riaz_turnover[,cmn.colnanes]) %>% 
  filter(overall != "NA") %>% 
  plotExpandedBox() 
ggsave("results/expansion/plots/expansion_space/total_mnc.pdf", width = 4, height = 4)


rbind(tumeh_turnover, riaz_turnover[,cmn.colnanes]) %>% 
  ggplot(aes(io_stat, followup, fill = io_stat)) + 
  
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  
  geom_jitter(size = 0.5, color = "black") + 
  ggsignif::geom_signif(comparisons = list(c("IO.naive", "Prior.IO")), color = "black") +
  labs(x = "", y = "expanded cells after therapy") +
  scale_color_manual(values = brewer.pal(2, 'Set1')) + 
  scale_fill_manual(values = brewer.pal(2, 'Set1')) + 
  theme(legend.position = "none")

ggsave("results/expansion/plots/expansion_space/io_stat.pdf", width = 4, height = 4)




#### Look for the TCRGP predictions for the clonotypes that stay
getOverlap <- function(row, folder){
  
  ## Outputs a summarised df 
  
  message(row)
  
  file1 = row[1]
  file2 = row[2]
  
  df1 = paste0(folder, file1) %>% fread %>% mutate(id = paste0(cdr3nt, v, j))
  df2 = paste0(folder, file2) %>% fread %>% mutate(id = paste0(cdr3nt, v, j))
  
  ## Because MNC files don't always contain v or j or the nt are weird, let's use aa's instead
  df1 = aggregate(freq ~ cdr3aa, df1, sum)
  df2 = aggregate(freq ~ cdr3aa, df2, sum)
  
  clonotypes1 = nrow(df1)
  clonotypes2 = nrow(df2)
  
  overlap = merge(df1, df2, by = "cdr3aa") %>% mutate(id = extractFileName(file1))
  return(overlap)
  
}

tumeh_turnover_df  <- fread("results/expansion/tumeh_sample_pairs.txt") %>% apply(1, FUN = getOverlap, folder = "data/unselected_TCRseq/Tumeh/") %>% rbindlist()
tumeh_turnover_df2 <- tumeh_turnover_df %>% mutate(temp = paste0(cdr3aa, id))
tumeh_predictions  <- tumeh_predictions %>% mutate(temp = paste0(cdr3aa, filename))

tumeh_remainers_tcrgp <- merge(tumeh_turnover_df2, tumeh_predictions, by = "temp")
tumeh_remainers_tcrgp <- tumeh_remainers_tcrgp[!duplicated(tumeh_remainers_tcrgp$temp), ]

riaz_turnover_df  <- fread("results/expansion/riaz_sample_pairs.txt") %>% apply(1, FUN = getOverlap, folder = "data/unselected_TCRseq/riaz/") %>% rbindlist()
riaz_turnover_df2 <- riaz_turnover_df %>% mutate(temp = paste0(cdr3aa, id))
riaz_predictions  <- riaz_predictions %>% mutate(temp = paste0(cdr3aa, filename))

riaz_remainers_tcrgp <- merge(riaz_turnover_df2, riaz_predictions, by = "temp")
riaz_remainers_tcrgp <- riaz_remainers_tcrgp[!duplicated(riaz_remainers_tcrgp$temp), ]

yusko_mnc_turnover_df  <- fread("results/expansion/yusko_mnc_sample_pairs.txt") %>% apply(1, FUN = getOverlap, folder = "data/unselected_TCRseq/Yusko/MNC_corrected/") %>% rbindlist()
yusko_mnc_turnover_df2 <- yusko_mnc_turnover_df %>% mutate(temp = paste0(cdr3aa, id))
yusko_mnc_predictions  <- yusko_mnc_predictions %>% mutate(temp = paste0(cdr3aa, filename))

yusko_mnc_remainers_tcrgp <- merge(yusko_mnc_turnover_df2, yusko_mnc_predictions, by = "temp")
yusko_mnc_remainers_tcrgp <- yusko_mnc_remainers_tcrgp[!duplicated(yusko_mnc_remainers_tcrgp$temp), ]


tota_remainers_tcrgp <- rbind(tumeh_remainers_tcrgp, riaz_remainers_tcrgp, yusko_mnc_remainers_tcrgp)

tota_remainers_tcrgp %>% group_by(name, species) %>% summarise(clonal_replacement = sum(freq.x)) %>% 
  filter(species %in% c("INFa", interesting_species)) %>% ungroup() %>% 
  tidyr::complete(name, species, fill = list(clonal_replacement = 0)) %>% 
  left_join(meta_clinical, by = c("name" = "Study.ID")) %>% filter(!is.na(response)) %>% 
  mutate(riaz = ifelse(response %in% c("CR", "PR"), "CR/PR", as.character(response))) %>% mutate(riaz = factor(as.character(riaz), levels = c("CR/PR", "SD", "PD"))) %>% 
  mutate(objective = ifelse(response %in% c("CR", "PR"), "R", "NR")) %>%  
  
  ggplot(aes(overall, clonal_replacement)) + geom_boxplot() + facet_wrap(~species) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~Previous.IO)




## Normalize, to have frequencies of clonally replaced cells
tota_remainers_tcrgp <- tota_remainers_tcrgp %>% group_by(id) %>% mutate(freq.x = freq.x / sum(freq.x))
tota_remainers_tcrgp <- tota_remainers_tcrgp %>% group_by(id) %>% mutate(freq.y = freq.y / sum(freq.y)) %>% 
  mutate(riaz = ifelse(response %in% c("CR", "PR"), "CR/PR", response)) %>% mutate(riaz = factor(as.character(riaz), levels = c("CR/PR", "SD", "PD")))

## Analyze the samples; freq of clonal replacement clonotypes
df_species <- tota_remainers_tcrgp %>% group_by(name, species) %>% summarise(clonal_replacement = sum(freq.y)) %>% 
  filter(species %in% c("INFa", interesting_species)) %>% ungroup() %>% 
  tidyr::complete(name, species, fill = list(clonal_replacement = 0)) %>% 
  left_join(meta_clinical, by = c("name" = "Study.ID")) %>% filter(!is.na(response)) 

df_species %>% ggplot(aes(species, clonal_replacement)) + geom_boxplot() + ggpubr::stat_compare_means(ref = "MELANOMA", label = "p.format")

df_species %>% ggplot(aes(riaz, clonal_replacement)) + geom_boxplot() + facet_wrap(~species, scales = "free") + ggpubr::stat_compare_means(ref = "PD", label = "p.signif")
df_species %>% ggplot(aes(objective, clonal_replacement)) + geom_boxplot() + facet_wrap(~species, scales = "free") + ggpubr::stat_compare_means(label = "p.format")
df_species %>% ggplot(aes(overall, clonal_replacement)) + geom_boxplot() + facet_wrap(~species, scales = "free") + ggpubr::stat_compare_means(label = "p.format")

## Analyze the samples; freq of stay-at-home clonotypes
df_species <- tota_remainers_tcrgp %>% group_by(name, species) %>% summarise(clonal_replacement = sum(freq.x)) %>% 
  filter(species %in% c("INFa", interesting_species)) %>% ungroup() %>% 
  tidyr::complete(name, species, fill = list(clonal_replacement = 0)) %>% 
  left_join(meta_clinical, by = c("name" = "Study.ID")) %>% filter(!is.na(response))

df_species %>% ggplot(aes(species, clonal_replacement)) + geom_boxplot() + ggpubr::stat_compare_means(ref = "MELANOMA", label = "p.format")

df_species %>% ggplot(aes(riaz, clonal_replacement)) + geom_boxplot() + facet_wrap(~species, scales = "free") + ggpubr::stat_compare_means(ref = "PD", label = "p.signif")
df_species %>% ggplot(aes(objective, clonal_replacement)) + geom_boxplot() + facet_wrap(~species, scales = "free") + ggpubr::stat_compare_means(label = "p.format")
df_species %>% ggplot(aes(overall, clonal_replacement)) + geom_boxplot() + facet_wrap(~species, scales = "free") + ggpubr::stat_compare_means(label = "p.format")


## Analyze the samples; n of stay-at-home clonotypes 
df_species <- tota_remainers_tcrgp %>% group_by(name, species) %>% summarise(clonal_replacement = n()) %>% 
  filter(species %in% c("INFa", interesting_species)) %>% ungroup() %>% 
  tidyr::complete(name, species, fill = list(clonal_replacement = 0)) %>% 
  left_join(meta_clinical, by = c("name" = "Study.ID")) %>% filter(!is.na(response))

df_species %>% ggplot(aes(species, clonal_replacement)) + geom_boxplot() + ggpubr::stat_compare_means(ref = "MELANOMA", label = "p.format")

df_species %>% ggplot(aes(riaz, clonal_replacement)) + geom_boxplot() + facet_wrap(~species, scales = "free") + ggpubr::stat_compare_means(ref = "PD", label = "p.signif")
df_species %>% ggplot(aes(objective, clonal_replacement)) + geom_boxplot() + facet_wrap(~species, scales = "free") + ggpubr::stat_compare_means(label = "p.format")
df_species %>% ggplot(aes(overall, clonal_replacement)) + geom_boxplot() + facet_wrap(~species, scales = "free") + ggpubr::stat_compare_means(label = "p.format")



tota_remainers_tcrgp %>% group_by(name, overall) %>% summarise(n = n()) %>% ggplot(aes(overall, n)) + geom_boxplot() + ggsignif::geom_signif(comparisons = list(c("N", "R")))
tota_remainers_tcrgp %>% group_by(name, riaz) %>% summarise(n = n()) %>% ggplot(aes(riaz, n)) + geom_boxplot() + ggsignif::geom_signif(comparisons = list(c("CR/PR", "SD"), c("SD", "PD"), c("CR/PR", "PD")), step_increase = 0.05)

tota_remainers_tcrgp %>% group_by(species, name, overall) %>% summarise(n = n()) %>% ggplot(aes(overall, n)) + geom_boxplot() + facet_wrap(~species, scales = "free_y", ncol = 4) + ggsignif::geom_signif(comparisons = list(c("N", "R")))
tota_remainers_tcrgp %>% mutate(overall = ifelse(response == "SD", "N", overall)) %>% group_by(species, name, overall) %>% summarise(n = n()) %>% ggplot(aes(overall, n)) + geom_boxplot() + facet_wrap(~species, scales = "free_y", ncol = 4) + ggsignif::geom_signif(comparisons = list(c("N", "R")))
tota_remainers_tcrgp %>% mutate(overall = ifelse(response == "SD", "N", overall)) %>% group_by(overall, species, name) %>% summarise(n = n()) %>% mutate(prop = n/sum(n)) %>% ggplot(aes(overall, prop)) + geom_boxplot() + facet_wrap(~species, scales = "free_y", ncol = 4) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + geom_jitter()
tota_remainers_tcrgp %>% group_by(species, name, response) %>% summarise(n = n()) %>% ggplot(aes(response, n)) + geom_boxplot() + facet_wrap(~species, scales = "free_y", ncol = 4) + ggsignif::geom_signif(comparisons = list(c("N", "R")))

tota_remainers_tcrgp %>% group_by(name, overall, species) %>% summarise(n = n()) %>% mutate(prop = n/sum(n)) %>% 
  ggplot(aes(overall, prop)) + geom_boxplot() + facet_wrap(~species, scales = "free_y", ncol = 4) + ggsignif::geom_signif(comparisons = list(c("N", "R")))
