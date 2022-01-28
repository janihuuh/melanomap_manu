
## Biomarker of expansion, freq and n, absolute and fold-change, of anti-MAA/anti-viral clonotypes in the TIL, total, stay-at-home and clonal replacement

## Clonal replacement cells: clonotypes that are not seen before treatment (i.e. 1 - cmn.clonotypes)

# tumeh_clonal_replacement  <- fread("results/expansion/tumeh_sample_pairs.txt") %>% apply(1, FUN = getClonalReplacementClones, folder = "data/unselected_TCRseq/Tumeh/") %>% rbindlist() %>% mutate(temp = paste0(cdr3aa, name))
tumeh_clonal_replacement  <- fread("results/expansion/tumeh_sample_pairs.txt")[-11,] %>% apply(1, FUN = getClonalReplacementClones, folder = "data/unselected_tcrb/resampled/Tumeh/") %>% rbindlist() %>% mutate(temp = paste0(cdr3aa, name))
tumeh_predictions  <- tumeh_predictions %>% mutate(temp = paste0(cdr3aa, name))

tumeh_clonal_replacement_tcrgp <- merge(tumeh_clonal_replacement, tumeh_predictions, by = "temp")
tumeh_clonal_replacement_tcrgp <- tumeh_clonal_replacement_tcrgp[!duplicated(tumeh_clonal_replacement_tcrgp$temp), ]

# riaz_clonal_replacement  <- fread("results/expansion/riaz_sample_pairs.txt")[-c(2,5,8,9,11,14),] %>% apply(1, FUN = getClonalReplacementClones, folder = "data/unselected_tcrb/resampled/riaz/") %>% rbindlist() %>% mutate(temp = paste0(cdr3aa, name))
riaz_clonal_replacement  <- fread("results/expansion/riaz_sample_pairs.txt")[-c(2,5,8,9,11,14),] %>% apply(1, FUN = getClonalReplacementClones, folder = "data/unselected_tcrb/resampled/riaz/") %>% rbindlist() %>% mutate(temp = paste0(cdr3aa, name))
riaz_predictions  <- riaz_predictions %>% mutate(temp = paste0(cdr3aa, name))

riaz_clonal_replacement_tcrgp <- merge(riaz_clonal_replacement, riaz_predictions, by = "temp")
riaz_clonal_replacement_tcrgp <- riaz_clonal_replacement_tcrgp[!duplicated(riaz_clonal_replacement_tcrgp$temp), ]

# yusko_mnc_clonal_replacement  <- fread("results/expansion/yusko_mnc_sample_pairs.txt") %>% apply(1, FUN = getClonalReplacementClones, folder = "data/unselected_TCRseq/Yusko/MNC_corrected//") %>% rbindlist() %>% mutate(temp = paste0(cdr3aa, name))
yusko_mnc_clonal_replacement  <- fread("results/expansion/yusko_mnc_sample_pairs.txt") %>% apply(1, FUN = getClonalReplacementClones, folder = "data/unselected_tcrb/resampled/Yusko_mnc/") %>% rbindlist() %>% mutate(temp = paste0(cdr3aa, name))
yusko_mnc_predictions  <- yusko_mnc_predictions %>% mutate(temp = paste0(cdr3aa, name))

yusko_mnc_clonal_replacement_tcrgp <- merge(yusko_mnc_clonal_replacement, yusko_mnc_predictions, by = "temp")
yusko_mnc_clonal_replacement_tcrgp <- yusko_mnc_clonal_replacement_tcrgp[!duplicated(yusko_mnc_clonal_replacement_tcrgp$temp), ]

tot_clonalReplacement_tcrgp <- rbind(tumeh_clonal_replacement_tcrgp, riaz_clonal_replacement_tcrgp, yusko_mnc_clonal_replacement_tcrgp)



###### Freq of clonal replacement cells

dir.create("results/clonal_replacement/thorough/", showWarnings = F)
dir.create("results/clonal_replacement/thorough/exact/", showWarnings = F)
dir.create("results/clonal_replacement/thorough/resampled/", showWarnings = F)

## Not corrected for io stat 
tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(riaz, clonal_replacement_freq, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("SD", "PD"), c("CR/PR", "SD"), c("CR/PR", "PD")), step_increase = 0.05) + 
   scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_freq_riaz.png", width = 5, height = 4)

tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(objective, clonal_replacement_freq, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_freq_objective.png", width = 5, height = 4)

tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(overall, clonal_replacement_freq, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_freq_overall.png", width = 5, height = 4)


## Corrected for io stat 
tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(riaz, clonal_replacement_freq, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("SD", "PD"), c("CR/PR", "SD"), c("CR/PR", "PD")), step_increase = 0.05) + facet_wrap(~Previous.IO) + 
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_freq_riaz_io.png", width = 6, height = 4.5)

tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(objective, clonal_replacement_freq, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~Previous.IO) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_freq_objective_io.png", width = 6, height = 4.5)

tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(overall, clonal_replacement_freq, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~Previous.IO) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_freq_overall_io.png", width = 6, height = 4.5)


## Only IO-stat
tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>%  ggplot(aes(Previous.IO, clonal_replacement_freq, fill = Previous.IO)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("IO naive", "Prior IO")), step_increase = 0.05) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette3(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_freq_io.png", width = 5, height = 4)


## By species, not corrected for io stat
tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(Cohort, clonal_replacement_freq, fill = Cohort)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("Yusko", "Riaz"))) + facet_wrap(~species) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_cohort.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(riaz, clonal_replacement_freq, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("SD", "PD"), c("CR/PR", "SD"), c("CR/PR", "PD")), step_increase = 0.05) + facet_wrap(~species) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_riaz.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(objective, clonal_replacement_freq, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~species)+ geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_objective.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(overall, clonal_replacement_freq, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~species) +
  geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_overall.png", width = 6, height = 6)



## By species, normalized for clonally replaced cells, not for io stat
tot_clonalReplacement_tcrgp %>% 
  group_by(name.x) %>% mutate(freq.x = freq.x / sum(freq.x)) %>% ungroup() %>% 
  group_by(name.x, species) %>% 
  summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(Cohort, clonal_replacement_freq, fill = Cohort)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("Yusko", "Riaz"))) + facet_wrap(~species)+ geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_cohort_norm.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% 
  group_by(name.x) %>% mutate(freq.x = freq.x / sum(freq.x)) %>% ungroup() %>% 
  group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(riaz, clonal_replacement_freq, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("SD", "PD"), c("CR/PR", "SD"), c("CR/PR", "PD")), step_increase = 0.05) + facet_wrap(~species) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_riaz_norm.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% 
  group_by(name.x) %>% mutate(freq.x = freq.x / sum(freq.x)) %>% ungroup() %>% 
  group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(objective, clonal_replacement_freq, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~species) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_objective_norm.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% 
  group_by(name.x) %>% mutate(freq.x = freq.x / sum(freq.x)) %>% ungroup() %>% 
  group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(overall, clonal_replacement_freq, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~species) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_overall_norm.png", width = 6, height = 6)





## By species, corrected for io stat
tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(Cohort, clonal_replacement_freq, fill = Cohort)) + geom_boxplot(outlier.shape = NA) + facet_wrap(Previous.IO~species, ncol = 4) + ggsignif::geom_signif(comparisons = list(c("Yusko", "Riaz"))) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_cohort_io.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(riaz, clonal_replacement_freq, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("SD", "PD"), c("CR/PR", "SD"), c("CR/PR", "PD")), step_increase = 0.05) + facet_wrap(Previous.IO~species, ncol = 4) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_riaz_io.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(objective, clonal_replacement_freq, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(Previous.IO~species, ncol = 4) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_objective_io.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(overall, clonal_replacement_freq, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(Previous.IO~species, ncol = 4) +
  geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_overall_io.png", width = 6, height = 6)



## By species, normalized for clonally replaced cells, not for io stat
tot_clonalReplacement_tcrgp %>% 
  group_by(name.x) %>% mutate(freq.x = freq.x / sum(freq.x)) %>% ungroup() %>% 
  group_by(name.x, species) %>% 
  summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(Cohort, clonal_replacement_freq, fill = Cohort)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("Yusko", "Riaz"))) + facet_wrap(Previous.IO~species, ncol = 4) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_cohort_norm_io.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% 
  group_by(name.x) %>% mutate(freq.x = freq.x / sum(freq.x)) %>% ungroup() %>% 
  group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(riaz, clonal_replacement_freq, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("SD", "PD"), c("CR/PR", "SD"), c("CR/PR", "PD")), step_increase = 0.05) + facet_wrap(Previous.IO~species, ncol = 4) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_riaz_norm_io.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% 
  group_by(name.x) %>% mutate(freq.x = freq.x / sum(freq.x)) %>% ungroup() %>% 
  group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(objective, clonal_replacement_freq, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(Previous.IO~species, ncol = 4) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_objective_norm_io.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% 
  group_by(name.x) %>% mutate(freq.x = freq.x / sum(freq.x)) %>% ungroup() %>% 
  group_by(name.x, species) %>% summarise(clonal_replacement_freq = sum(freq.x)) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(overall, clonal_replacement_freq, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(Previous.IO~species, ncol = 4) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement freq") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_freq_overall_norm_io.png", width = 6, height = 6)







###### N of clonal replacement cells



## Not corrected for io stat 
tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(riaz, clonal_replacement_n, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("SD", "PD"), c("CR/PR", "SD"), c("CR/PR", "PD")), step_increase = 0.05) + 
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_n_riaz.png", width = 5, height = 4)

tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(objective, clonal_replacement_n, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_n_objective.png", width = 5, height = 4)

tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(overall, clonal_replacement_n, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_n_overall.png", width = 5, height = 4)


## Corrected for io stat 
tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(riaz, clonal_replacement_n, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("SD", "PD"), c("CR/PR", "SD"), c("CR/PR", "PD")), step_increase = 0.05) + facet_wrap(~Previous.IO) + 
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_n_riaz_io.png", width = 6, height = 4.5)

tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(objective, clonal_replacement_n, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~Previous.IO) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_n_objective_io.png", width = 6, height = 4.5)

tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  ggplot(aes(overall, clonal_replacement_n, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~Previous.IO) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_n_overall_io.png", width = 6, height = 4.5)


## Only IO-stat
tot_clonalReplacement_tcrgp %>% group_by(name.x) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>%  ggplot(aes(Previous.IO, clonal_replacement_n, fill = Previous.IO)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("IO naive", "Prior IO")), step_increase = 0.05) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette3(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_n_io.png", width = 5, height = 4)


## By species, not corrected for io stat
tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(Cohort, clonal_replacement_n, fill = Cohort)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("Yusko", "Riaz"))) + facet_wrap(~species) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_species_n_cohort.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(riaz, clonal_replacement_n, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("SD", "PD"), c("CR/PR", "SD"), c("CR/PR", "PD")), step_increase = 0.05) + facet_wrap(~species) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_n_riaz.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(objective, clonal_replacement_n, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~species)+ geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_n_objective.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(overall, clonal_replacement_n, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(~species) +
  geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_n_overall.png", width = 6, height = 6)




## By species, corrected for io stat
tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(Cohort, clonal_replacement_n, fill = Cohort)) + geom_boxplot(outlier.shape = NA) + facet_wrap(Previous.IO~species, ncol = 4) + ggsignif::geom_signif(comparisons = list(c("Yusko", "Riaz"))) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5)
ggsave("results/clonal_replacement/thorough/resampled/box_species_n_cohort_io.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(riaz, clonal_replacement_n, fill = riaz)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("SD", "PD"), c("CR/PR", "SD"), c("CR/PR", "PD")), step_increase = 0.05) + facet_wrap(Previous.IO~species, ncol = 4) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_n_riaz_io.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(objective, clonal_replacement_n, fill = objective)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(Previous.IO~species, ncol = 4) + geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_n_objective_io.png", width = 6, height = 6)

tot_clonalReplacement_tcrgp %>% group_by(name.x, species) %>% summarise(clonal_replacement_n = n()) %>% left_join(meta_clinical, by = c("name.x" = "Study.ID")) %>% 
  filter(!is.na(response)) %>% 
  filter(species %in% c(interesting_species_tcgrp, "INFa")) %>% 
  ggplot(aes(overall, clonal_replacement_n, fill = overall)) + geom_boxplot(outlier.shape = NA) + ggsignif::geom_signif(comparisons = list(c("N", "R"))) + facet_wrap(Previous.IO~species, ncol = 4) +
  geom_jitter(size = 0.3) + scale_fill_manual(values = getPalette(4)) +
  scale_fill_manual(values = getPalette3(4)) + labs(x = "", y = "clonal replacement n") + theme_classic(base_size = 17) + theme(legend.position = "none") + geom_jitter(size = 0.5) + facets_nice
ggsave("results/clonal_replacement/thorough/resampled/box_species_n_overall_io.png", width = 6, height = 6)


## The expansion, freq and n, absolute and fold-change, of anti-MAA/anti-viral clonotypes in the TIL, total, stay-at-home and clonal replacement
## Following the simultaneous expansion in blood and recruitment to tumor site in Yusko cohort
