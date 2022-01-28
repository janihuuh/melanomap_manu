
## Fold change for diversity metrics
helsinki_fc         <- preprocessFoldchange(helsinki_diversity_0m, helsinki_diversity_3m, clin_df = helsinki_clin)
tumeh_fc            <- preprocessFoldchange(tumeh_diversity_0m, tumeh_diversity_3m, clin_df = tumeh_clin)
riaz_fc             <- preprocessFoldchange(riaz_diversity_0m, riaz_diversity_3m, clin_df = riaz_clin)
robert_fc           <- preprocessFoldchange(robert_diversity_0m, robert_diversity_3m, clin_df = robert_clin)

yusko_mnc_IpNi_fc   <- preprocessFoldchange(df1 = yusko_mnc_IpNi_0m, df2 = yusko_mnc_IpNi_3m, clin_df = yusko_mnc_clin)
yusko_mnc_NiIp_fc   <- preprocessFoldchange(yusko_mnc_NiIp_0m, yusko_mnc_NiIp_3m, clin_df = yusko_mnc_clin)

yusko_cd8_IpNi_fc   <- preprocessFoldchange(yusko_cd8_IpNi_0m, yusko_cd8_IpNi_3m, clin_df = yusko_cd8_clin)
yusko_cd8_NiIp_fc   <- preprocessFoldchange(yusko_cd8_NiIp_0m, yusko_cd8_NiIp_3m, clin_df = yusko_cd8_clin)


plotFoldchange(helsinki_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/helsinki.pdf', width = 12, height = 12)

plotFoldchange(tumeh_fc)    + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/tumeh.pdf', width = 12, height = 12)

plotFoldchange(riaz_fc)     + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/riaz.pdf', width = 12, height = 12)

plotFoldchange(robert_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/robert.pdf', width = 12, height = 12)

plotFoldchange(yusko_mnc_IpNi_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/yusko_mnc_IpNi.pdf', width = 12, height = 12)

plotFoldchange(yusko_mnc_NiIp_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/yusko_mnc_NiIp.pdf', width = 12, height = 12)

plotFoldchange(yusko_cd8_IpNi_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/yusko_cd8_IpNi.pdf', width = 12, height = 12)

plotFoldchange(yusko_cd8_NiIp_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/yusko_cd8_NiIp.pdf', width = 12, height = 12)


## Combine pd1 
rbind(tumeh_fc, riaz_fc) %>% plotFoldchange + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/pd1.pdf', width = 12, height = 12)

rbind(tumeh_fc, riaz_fc) %>% filter(io_stat == "Prior.IO") %>% plotFoldchange + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/pd1_prior_io.pdf', width = 12, height = 12)

rbind(tumeh_fc, riaz_fc) %>% filter(io_stat == "IO.naive") %>% plotFoldchange + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/pd1_io_naive.pdf', width = 12, height = 12)






## Study the effects of regimens, irrespectively to the response
plotFoldchangeTimepoint(helsinki_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/helsinki_regimen.pdf', width = 12, height = 12)

plotFoldchangeTimepoint(tumeh_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/tumeh_regimen.pdf', width = 12, height = 12)

plotFoldchangeTimepoint(riaz_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/riaz_regimen.pdf', width = 12, height = 12)

plotFoldchangeTimepoint(robert_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/robert_regimen.pdf', width = 12, height = 12)

plotFoldchangeTimepoint(yusko_mnc_IpNi_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/yusko_mnc_IpNi_regimen.pdf', width = 12, height = 12)

plotFoldchangeTimepoint(yusko_mnc_NiIp_fc) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/yusko_mnc_NiIp_regimen.pdf', width = 12, height = 12)



## Pool samples
rbind(tumeh_fc, riaz_fc) %>% plotFoldchangeTimepoint() + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/pd1_regimen.pdf', width = 12, height = 12)



## After the therapy, which regimen produces the largest differences?


## Pool the tumor samples together
rbind(tumeh_fc, riaz_fc, yusko_mnc_IpNi_fc, yusko_mnc_NiIp_fc) %>% 
  plotFoldchangeComparisons(map_signif_level = T) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/tumor_effects_regimen.pdf', width = 12, height = 12)

rbind(tumeh_fc, riaz_fc, yusko_mnc_IpNi_fc, yusko_mnc_NiIp_fc)  %>% filter(overall == "R") %>% 
  plotFoldchangeComparisons + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/tumor_effects_regimen_R.pdf', width = 12, height = 12)

rbind(tumeh_fc, riaz_fc, yusko_mnc_IpNi_fc, yusko_mnc_NiIp_fc)  %>% filter(overall == "N") %>% 
  plotFoldchangeComparisons + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/tumor_effects_regimen_N.pdf', width = 12, height = 12)


## Study the clonality more closely
rbind(tumeh_fc, riaz_fc, yusko_mnc_IpNi_fc, yusko_mnc_NiIp_fc) %>% filter(variable == "normalizedShannonWienerIndex_mean") %>% 
  plotFoldchangeComparisons
ggsave('results/diversity/foldchange/tumor_effects_regimen_clonality.pdf', width = 4, height = 3)

rbind(tumeh_fc, riaz_fc, yusko_mnc_IpNi_fc, yusko_mnc_NiIp_fc) %>% filter(variable == "normalizedShannonWienerIndex_mean") %>% 
  filter(overall == "R") %>% 
  plotFoldchangeComparisons 
ggsave('results/diversity/foldchange/tumor_effects_regimen_clonality.pdf', width = 4, height = 3)
ggsave('results/diversity/foldchange/tumor_effects_regimen_clonality_R.pdf', width = 4, height = 3)

rbind(tumeh_fc, riaz_fc, yusko_mnc_IpNi_fc, yusko_mnc_NiIp_fc) %>% filter(variable == "normalizedShannonWienerIndex_mean") %>% 
  filter(overall == "N") %>% 
  plotFoldchangeComparisons 
ggsave('results/diversity/foldchange/tumor_effects_regimen_clonality.pdf', width = 4, height = 3)
ggsave('results/diversity/foldchange/tumor_effects_regimen_clonality_N.pdf', width = 4, height = 3)





## Pool the blood samples together
rbind(yusko_cd8_IpNi_fc, yusko_cd8_NiIp_fc, robert_fc, helsinki_fc) %>% 
  plotFoldchangeComparisons(calculate_y_pos = T, map_signif_level = T) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/blood_effects_regimen.pdf', width = 12, height = 12)

rbind(yusko_cd8_IpNi_fc, yusko_cd8_NiIp_fc, robert_fc, helsinki_fc)  %>% filter(overall == "R") %>% 
  plotFoldchangeComparisons(calculate_y_pos = T, map_signif_level = T) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/blood_effects_regimen_R.pdf', width = 12, height = 12)

rbind(yusko_cd8_IpNi_fc, yusko_cd8_NiIp_fc, robert_fc, helsinki_fc) %>% filter(overall == "N") %>% 
  plotFoldchangeComparisons(calculate_y_pos = T, map_signif_level = T) + facet_wrap(~variable, scales = 'free') + facets_nice
ggsave('results/diversity/foldchange/blood_effects_regimen_N.pdf', width = 12, height = 12)






## Study the different sample types in Yusko cohort; is there a difference?
rbind(yusko_mnc_IpNi_fc, yusko_cd8_IpNi_fc) %>% 
  plotFoldchangeComparisons(x_axis = 'type', plotSigf = F) + facet_wrap(~variable, scales = 'free') + facets_nice +
  ggsignif::geom_signif(comparisons = list(c("Blood", "Tumor")))
ggsave('results/diversity/foldchange/IpNi_tumor_vs_blood.pdf', width = 12, height = 12)

rbind(yusko_mnc_IpNi_fc, yusko_cd8_IpNi_fc) %>% filter(overall == "R") %>% 
  plotFoldchangeComparisons(x_axis = 'type', plotSigf = F) + facet_wrap(~variable, scales = 'free') + facets_nice +
  ggsignif::geom_signif(comparisons = list(c("Blood", "Tumor")))
ggsave('results/diversity/foldchange/IpNi_tumor_vs_blood_R.pdf', width = 12, height = 12)

rbind(yusko_mnc_IpNi_fc, yusko_cd8_IpNi_fc) %>% filter(overall == "N") %>% 
  plotFoldchangeComparisons(x_axis = 'type', plotSigf = F) + facet_wrap(~variable, scales = 'free') + facets_nice +
  ggsignif::geom_signif(comparisons = list(c("Blood", "Tumor")))
ggsave('results/diversity/foldchange/IpNi_tumor_vs_blood_N.pdf', width = 12, height = 12)

a = rbind(yusko_mnc_IpNi_fc, yusko_cd8_IpNi_fc) %>% filter(variable == "inverseSimpsonIndex_mean") %>% filter(overall == "N")
wilcox.test(filter(a, type == "Blood")$foldchange, filter(a, type == "Tumor")$foldchange)




rbind(yusko_mnc_NiIp_fc, yusko_cd8_NiIp_fc) %>% 
  plotFoldchangeComparisons(x_axis = 'type', plotSigf = F) + facet_wrap(~variable, scales = 'free') + facets_nice +
  ggsignif::geom_signif(comparisons = list(c("Blood", "Tumor")))
ggsave('results/diversity/foldchange/NiIp_tumor_vs_blood.pdf', width = 12, height = 12)

rbind(yusko_mnc_NiIp_fc, yusko_cd8_NiIp_fc) %>% filter(overall == "R") %>% 
  plotFoldchangeComparisons(x_axis = 'type', plotSigf = F) + facet_wrap(~variable, scales = 'free') + facets_nice +
  ggsignif::geom_signif(comparisons = list(c("Blood", "Tumor")))
ggsave('results/diversity/foldchange/NiIp_tumor_vs_blood_R.pdf', width = 12, height = 12)

rbind(yusko_mnc_NiIp_fc, yusko_cd8_NiIp_fc) %>% filter(overall == "N") %>% 
  plotFoldchangeComparisons(x_axis = 'type', plotSigf = F) + facet_wrap(~variable, scales = 'free') + facets_nice +
  ggsignif::geom_signif(comparisons = list(c("Blood", "Tumor")))
ggsave('results/diversity/foldchange/NiIp_tumor_vs_blood_N.pdf', width = 12, height = 12)


## Plot paired analysis
IpNi_fc = merge(yusko_mnc_IpNi_fc, yusko_cd8_IpNi_fc, by = "merge_id")
IpNi_fc %>% mutate(foldchange = foldchange.x / foldchange.y) %>% dplyr::rename(overall = overall.x, variable = variable.x) %>% 
  plotFoldchange + facet_wrap(~variable, scales = 'free') + facets_nice 
ggsave('results/diversity/foldchange/IpNi_tumor_vs_blood_paired.pdf', width = 12, height = 12)

NiIp_fc = merge(yusko_mnc_NiIp_fc, yusko_cd8_NiIp_fc, by = "merge_id")
NiIp_fc %>% mutate(foldchange = foldchange.x / foldchange.y) %>% dplyr::rename(overall = overall.x, variable = variable.x) %>% 
  plotFoldchange + facet_wrap(~variable, scales = 'free') + facets_nice 
ggsave('results/diversity/foldchange/NiIp_tumor_vs_blood_paired.pdf', width = 12, height = 12)




