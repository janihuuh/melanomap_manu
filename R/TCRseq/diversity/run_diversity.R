
## Analyze the diversity in ACT samples

## CD8
cd8_infusion_diversity = read.delim('results/diversity/calculated/cd8_infusion.diversity.aa.exact.txt') %>% 
  mutate(name = extractName(sample_id)) %>% 
  left_join(act_clin_df) %>% 
  left_join(til_to_cd8_expanded_group) %>% 
  preprocess_for_wilcox.test

cd8_infusion_corr = corr.test.list(cd8_infusion_diversity)
write.table(cd8_infusion_diversity, 'results/diversity/correlation/cd8_infusion_diversity.txt', sep = '\t', quote = F, row.names = F)

cd8_infusion_diversity = lapply(colnames(cd8_infusion_diversity), wilcox.test.list, df = cd8_infusion_diversity)
cd8_infusion_diversity = do.call(rbind, cd8_infusion_diversity) %>% 
  mutate(p.adj = p.adjust(p.value, method = 'BH')) %>% arrange(p.value) %>% 
  select(p.value, p.adj, name, mean_r, mean_n, up) %>% 
  mutate(p.value = round(p.value, 3), p.adj = round(p.adj, 3),   mean_r = round(mean_r, 3),
         mean_n = round(mean_n, 3))
cd8_infusion_diversity
write.table(cd8_infusion_diversity, 'results/diversity/wilcox/cd8_infusion_diversity.txt', sep = '\t', quote = F, row.names = F)


## CD8, BTLA+
cd8_btlapos_infusion_diversity = read.delim('results/diversity/calculated/CD8_btlapos_infusion.diversity.aa.exact.txt') %>% 
  mutate(name = extractName(sample_id)) %>% 
  left_join(act_clin_df) %>% 
  left_join(til_to_cd8_expanded_group) %>% 
  preprocess_for_wilcox.test

cd8_btlapos_infusion_corr = corr.test.list(cd8_btlapos_infusion_diversity)
write.table(cd8_btlapos_infusion_corr, 'results/diversity/correlation/cd8_btlapos_infusion_corr.txt', sep = '\t', quote = F, row.names = F)

cd8_btlapos_infusion_diversity = lapply(colnames(cd8_btlapos_infusion_diversity), wilcox.test.list, df = cd8_btlapos_infusion_diversity)
cd8_btlapos_infusion_diversity = do.call(rbind, cd8_btlapos_infusion_diversity) %>% 
  mutate(p.adj = p.adjust(p.value, method = 'BH')) %>% arrange(p.value) %>% 
  select(p.value, p.adj, name, mean_r, mean_n, up) %>% 
  mutate(p.value = round(p.value, 3), p.adj = round(p.adj, 3),  mean_r = round(mean_r, 3),
         mean_n = round(mean_n, 3))
cd8_btlapos_infusion_diversity
write.table(cd8_btlapos_infusion_diversity, 'results/diversity/wilcox/cd8_btlapos_infusion_diversity.txt', sep = '\t', quote = F, row.names = F)


## CD8, BTLA-
cd8_btlaneg_infusion_diversity = read.delim('results/diversity/calculated/CD8_btlaneg_infusion.diversity.aa.exact.txt') %>% 
  mutate(name = extractName(sample_id)) %>% 
  left_join(act_clin_df) %>% 
  preprocess_for_wilcox.test

cd8_btlaneg_infusion_corr = corr.test.list(cd8_btlaneg_infusion_diversity)
write.table(cd8_btlaneg_infusion_corr, 'results/diversity/correlation/cd8_btlaneg_infusion_corr.txt', sep = '\t', quote = F, row.names = F)


cd8_btlaneg_infusion_diversity = lapply(colnames(cd8_btlaneg_infusion_diversity), wilcox.test.list, df = cd8_btlaneg_infusion_diversity)
cd8_btlaneg_infusion_diversity = do.call(rbind, cd8_btlaneg_infusion_diversity) %>% 
  mutate(p.adj = p.adjust(p.value, method = 'BH')) %>% arrange(p.value) %>% 
  select(p.value, p.adj, name, mean_r, mean_n, up) %>% 
  mutate(p.value = round(p.value, 3), 
         p.adj = round(p.adj, 3),
         mean_r = round(mean_r, 3),
         mean_n = round(mean_n, 3))
cd8_btlaneg_infusion_diversity
write.table(cd8_btlaneg_infusion_diversity, 'results/diversity/wilcox/cd8_btlaneg_infusion_diversity.txt', sep = '\t', quote = F, row.names = F)

