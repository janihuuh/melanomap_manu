
source("src/jani/R/tcrb/expanded/load_expanded.R")
source("src/jani/R/tcrb/expanded/fun_plotExpanded.R")

## Basic dotplot
helsinki_exp  %>% plot_dot
ggsave(paste0("results/expansion/plots/helsinki_exp.pdf"), width = 6, height = 4)

tumeh_exp     %>% plot_dot
ggsave(paste0("results/expansion/plots/tumeh_exp.pdf"), width = 6, height = 4)

yusko_mnc_exp %>% plot_dot
ggsave(paste0("results/expansion/plots/yusko_mnc_exp.pdf"), width = 6, height = 4)

yusko_cd8_exp %>% plot_dot
ggsave(paste0("results/expansion/plots/yusko_cd8_exp.pdf"), width = 6, height = 4)

robert_exp    %>% plot_dot
ggsave(paste0("results/expansion/plots/robert_exp.pdf"), width = 6, height = 4)

riaz_exp    %>% plot_dot
ggsave(paste0("results/expansion/plots/riaz_exp.pdf"), width = 6, height = 4)


## Plot the amount of expanded clonotypes per patient as a rank plot
hel_exp_rank <- helsinki_expanded %>% filter(type == "Blood") %>% distinct(clonotypename, .keep_all = TRUE) %>% 
  group_by(name, .drop = FALSE) %>% dplyr::summarise(n = n()) %>% left_join(helsinki_clin, by = "name") 
hel_exp_rank %>% plot_rank + facet_wrap(~io_stat)
ggsave(paste0("results/expansion/plots/n_expanded_helsinki.pdf"), width = 12, height = 4)

tum_exp_rank <- tumeh_expanded %>% 
  distinct(clonotypename, .keep_all = TRUE) %>% 
  group_by(name, .drop = FALSE) %>% dplyr::summarise(n = n()) %>% left_join(tumeh_clin, by = "name") 
tum_exp_rank %>% plot_rank + facet_wrap(~io_stat)
ggsave(paste0("results/expansion/plots/n_expanded_tumeh.pdf"), width = 12, height = 4)

yus_mnc_exp_rank <- yusko_mnc_expanded %>% filter(!is.na(overall) & overall != "NA") %>% distinct(clonotypename, .keep_all = TRUE) %>% 
  group_by(name, .drop = FALSE) %>% dplyr::summarise(n = n()) %>% left_join(yusko_mnc_clin, by = "name") 
yus_mnc_exp_rank %>% plot_rank + facet_wrap(~regimen, scales = 'free_x') + facets_nice
ggsave(paste0("results/expansion/plots/n_expanded_yusko_mnc.pdf"), width = 12, height = 4)

yus_cd8_exp_rank <- yusko_cd8_expanded %>% filter(!is.na(overall) & overall != "NA") %>% distinct(clonotypename, .keep_all = TRUE) %>% 
  group_by(name, .drop = FALSE) %>% dplyr::summarise(n = n()) %>% left_join(yusko_cd8_clin, by = "name")
yus_cd8_exp_rank %>% plot_rank + facet_wrap(~regimen, scales = 'free_x') + facets_nice
ggsave(paste0("results/expansion/plots/n_expanded_yusko_cd8.pdf"), width = 12, height = 4)

rob_exp_rank <- robert_expanded %>% distinct(clonotypename, .keep_all = TRUE) %>% 
  group_by(name, .drop = FALSE) %>% dplyr::summarise(n = n()) %>% left_join(robert_clin, by = "name") 
rob_exp_rank %>% plot_rank
ggsave(paste0("results/expansion/plots/n_expanded_robert.pdf"), width = 6, height = 4)

ria_exp_rank <- riaz_expanded %>% filter(overall != 'NA') %>% 
  distinct(clonotypename, .keep_all = TRUE) %>% 
  group_by(name, .drop = FALSE) %>% dplyr::summarise(n = n()) %>% left_join(riaz_clin, by = "name") 
ria_exp_rank %>% plot_rank + facet_wrap(~io_stat)
ggsave(paste0("results/expansion/plots/n_expanded_riaz.pdf"), width = 12, height = 4)


## Meta
rbind(hel_exp_rank, tum_exp_rank, yus_mnc_exp_rank) %>% plot_rank
ggsave(paste0("results/expansion/plots/n_expanded_meta_hel_tum_yus_mnc.pdf"), width = 6, height = 4)

rbind(tum_exp_rank, ria_exp_rank) %>% plot_rank + facet_wrap(regimen~io_stat, scales = 'free_x') + facets_nice
ggsave(paste0("results/expansion/plots/n_expanded_meta_tumor_mnc.pdf"), width = 8, height = 4)





## Basic boxplot
hel_exp_rank %>%  plot_box + facet_wrap(~io_stat) + facets_nice
ggsave(paste0("results/expansion/plots/box_expanded_helsinki.pdf"), width = 6, height = 4)

tum_exp_rank  %>% plot_box + facet_wrap(~io_stat) + facets_nice
ggsave(paste0("results/expansion/plots/box_expanded_tumeh.pdf"),  width = 6, height = 4)

yusko_mnc_expanded %>% group_by(name1,overall,regimen) %>% dplyr::summarise(n = n()) %>% 
  plot_box() + facet_wrap(~regimen, scales = "free") + facets_nice
ggsave(paste0("results/expansion/plots/box_expanded_yusko_mnc.pdf"),  width = 6, height = 4)

yusko_cd8_expanded %>% filter(!is.na(overall)) %>% group_by(name1,overall,regimen) %>% dplyr::summarise(n = n()) %>% 
  plot_box() + facet_wrap(~regimen, scales = "free") + facets_nice
ggsave(paste0("results/expansion/plots/box_expanded_yusko_cd8.pdf"),  width = 6, height = 4)

rob_exp_rank    %>% plot_box
ggsave(paste0("results/expansion/plots/box_expanded_robert.pdf"),  width = 3, height = 4)

ria_exp_rank    %>% plot_box + facet_wrap(~io_stat)  + facets_nice
ggsave(paste0("results/expansion/plots/box_expanded_riaz.pdf"),  width = 6, height = 4)


## Pool samples together
rbind(yus_mnc_exp_rank, tum_exp_rank, ria_exp_rank) %>% 
  plot_box() + facet_wrap(io_stat~regimen, ncol = 4) + facets_nice
ggsave(paste0("results/expansion/plots/box_expanded_tumor.pdf"),  width = 12, height = 4)

rbind(tum_exp_rank, ria_exp_rank) %>% filter(regimen == 'antiPD1') %>% 
  ggplot(aes(io_stat, n, fill = io_stat)) + geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_jitter(size = 0.3) + 
  ggsignif::geom_signif(comparisons = list(c('Prior.IO','IO.naive'))) +
  scale_fill_manual(values = brewer.pal(2, 'Set1')) + labs(fill = 'IO status', x = '', y = 'nExpanded clonotypes')
ggsave(paste0("results/expansion/plots/box_expanded_aPD1_iostat.pdf"),  width = 6, height = 4)

hel_exp_rank %>% 
  ggplot(aes(io_stat, n, fill = io_stat)) + geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_jitter(size = 0.3) + 
  ggsignif::geom_signif(comparisons = list(c('Prior.IO','IO.naive'))) +
  scale_fill_manual(values = brewer.pal(2, 'Set1')) + labs(fill = 'IO status', x = '', y = 'nExpanded clonotypes')
ggsave(paste0("results/expansion/plots/box_expanded_alag3_iostat.pdf"),  width = 6, height = 4)

rbind(yus_mnc_exp_rank, tum_exp_rank, ria_exp_rank) %>% 
  filter(io_stat == 'IO.naive') %>% 
  ggplot(aes(regimen, n, fill = regimen)) + geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_jitter(size = 0.3) + 
  ggsignif::geom_signif(comparisons = list(c('antiCTLA4+antiPD1','antiPD1'), 
                                           c('antiPD1+antiCTLA4','antiPD1'),
                                           c('antiCTLA4+antiPD1','antiPD1+antiCTLA4'))) +
  scale_fill_manual(values = brewer.pal(3, 'Set1')) + facet_wrap(~io_stat) + labs(x = '', y = 'nExpanded clonotypes') + theme(legend.position = 'none') + facets_nice
ggsave("results/expansion/plots/box_expanded_tumor_regimen.pdf",  width = 6, height = 4)



  
  
## Correlate to the amount neoantigens and B2M loss in the Yusko-datasets
yusko_clin_full_individual <- yusko_clin_full %>% distinct(name, .keep_all = TRUE)
yusko_clin_full_individual = merge(yusko_mnc_clin, yusko_clin_full_individual)

## Yusko CD8
yusko_cd8_neo_exp <- yusko_cd8_expanded %>% distinct(clonotypename, .keep_all = TRUE) %>% 
  group_by(name, .drop = FALSE) %>% dplyr::summarise(n = n()) %>% 
  left_join(yusko_clin_full_individual, by = "name") 
  
yusko_cd8_neo_exp %>% 
  ggplot(aes(Neoatingen.Load,n,, size = Mutation.Load, color = B2M.Loss)) + geom_point() + facet_wrap(~Treatment.Group, scales = "free") + scale_color_manual(values = brewer.pal(3,"Set1")) +
  geom_smooth(method = "lm", fill = NA)
ggsave(paste0("results/expansion/plots/neoantigen_yusko_cd8.pdf"), width = 12, height = 4)


## Yusko MNC
yusko_mnc_neo_exp <- yusko_mnc_expanded %>% distinct(clonotypename, .keep_all = TRUE) %>% 
  group_by(name, .drop = FALSE) %>% dplyr::summarise(n = n()) %>% 
  left_join(yusko_clin_full_individual, by = "name")
  
yusko_mnc_neo_exp %>% 
  ggplot(aes(Neoatingen.Load,n, size = Mutation.Load, color = B2M.Loss)) + geom_point() + facet_wrap(~Treatment.Group, scales = "free") + scale_color_manual(values = brewer.pal(3,"Set1")) +
  geom_smooth(method = "lm", fill = NA)
ggsave(paste0("results/expansion/plots/neoantigen_yusko_mnc.pdf"), width = 8, height = 4)



## Riaz MNC
riaz_mnc_neo_exp <- riaz_expanded %>% distinct(clonotypename, .keep_all = TRUE) %>% 
  group_by(name, .drop = FALSE) %>% dplyr::summarise(n = n()) %>% 
  right_join(riaz_clin_full, by = 'name') %>% dplyr::rename(Neoatingen.Load = Neo.antigen.Load)

riaz_mnc_neo_exp %>%
  ggplot(aes(Neoatingen.Load, n, size = Mutation.Load)) + geom_point() + scale_color_manual(values = brewer.pal(3,"Set1")) +
  geom_smooth(method = "lm", fill = NA) + labs(y = 'nExpanded clonotypes') + facet_wrap(~io_stat)
ggsave(paste0("results/expansion/plots/neoantigen_riaz.pdf"), width = 8, height = 4)



## Combine
cmn.colnames = intersect(colnames(yusko_mnc_neo_exp), colnames(riaz_mnc_neo_exp))
total_neo = rbind(riaz_mnc_neo_exp[,cmn.colnames], yusko_mnc_neo_exp[,cmn.colnames]) 
  
total_neo %>% 
  ggplot(aes(Neoatingen.Load, n, size = Mutation.Load)) + geom_point() + scale_color_manual(values = brewer.pal(3,"Set1")) +
  geom_smooth(method = "lm", fill = NA) + labs(y = 'nExpanded clonotypes')
ggsave('results/expansion/plots/neoantigen_total.pdf', width = 6, height = 4)


total_neo %>% 
  ggplot(aes(Neoatingen.Load, n)) + geom_point(pch = 21, fill = 'lightgrey', size = 2) + scale_color_manual(values = brewer.pal(3,"Set1")) +
  geom_smooth(method = "lm", fill = NA, color = 'darkred') + labs(y = 'nExpanded clonotypes') + facet_wrap(regimen~io_stat, ncol = 4) + facets_nice
ggsave('results/expansion/plots/neoantigen_total_regimen.pdf', width = 12, height = 4)




## Fitting linear models, excluding the B2M-mutations
a <- lm(Neoatingen.Load + Mutation.Load ~ n, data = subset(yusko_cd8_neo_exp, Treatment.Group == " Ipilimumab-Nivolumab" & B2M.Loss == " No B2M Loss")) %>% 
  broom::tidy() %>% mutate(name = "Yusko MNC, Ipi-Nivo")
b <- lm(Neoatingen.Load + Mutation.Load  ~ n, data = subset(yusko_cd8_neo_exp, Treatment.Group == " Nivolumab-Ipilimumab" & B2M.Loss == " No B2M Loss")) %>% 
  broom::tidy() %>% mutate(name = "Yusko MNC, Nivo-Ipi")

c <- lm(Neoatingen.Load + Mutation.Load ~ n, data = subset(yusko_mnc_neo_exp, Treatment.Group == " Ipilimumab-Nivolumab" & B2M.Loss == " No B2M Loss")) %>% 
  broom::tidy() %>% mutate(name = "Yusko CD8+, Ipi-Nivo")
d <- lm(Neoatingen.Load + Mutation.Load  ~ n, data = subset(yusko_mnc_neo_exp, Treatment.Group == " Nivolumab-Ipilimumab" & B2M.Loss == " No B2M Loss")) %>% 
  broom::tidy() %>% mutate(name = "Yusko CD8+, Nivo-Ipi")


lm(Neoatingen.Load ~ n +  regimen + io_stat, data = total_neo) %>%  broom::tidy() 




tot <- rbind(a,b,c,d) 
tot[,2:5] <- apply(tot[,2:5], 2, function(x)round(x,2))

tot %>% write.table(paste0("results/expansion/plots/neoantigen_yusko_lm.txt"), sep = "\t", quote = F, row.names = F)



