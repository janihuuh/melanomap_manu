
## Calculate on baseline
helsinki0m_wilcox  <- preprocess_for_wilcox.test(df = helsinki_diversity_0m, clin_df = helsinki_clin)
tumeh0m_wilcox     <- preprocess_for_wilcox.test(df = tumeh_diversity_0m, clin_df = tumeh_clin)
riaz0m_wilcox      <- preprocess_for_wilcox.test(df = riaz_diversity_0m, clin_df = riaz_clin) %>% filter(overall != "NA")
robert0m_wilcox    <- preprocess_for_wilcox.test(df = robert_diversity_0m, clin_df = robert_clin)

yusko_mnc_IpNi0m_wilcox <- preprocess_for_wilcox.test(df = yusko_mnc_IpNi_0m, clin_df = yusko_mnc_clin)
yusko_mnc_NiIp0m_wilcox <- preprocess_for_wilcox.test(df = yusko_mnc_NiIp_0m, clin_df = yusko_mnc_clin)

yusko_cd8_IpNi0m_wilcox <- preprocess_for_wilcox.test(df = yusko_cd8_IpNi_0m, clin_df = yusko_cd8_clin)
yusko_cd8_NiIp0m_wilcox <- preprocess_for_wilcox.test(df = yusko_cd8_NiIp_0m, clin_df = yusko_cd8_clin)

## Calculate
helsinki0m_wilcox_res <- wilcox.test.list(helsinki0m_wilcox) %>% polishDiversityResults
tumeh0m_wilcox_res    <- wilcox.test.list(tumeh0m_wilcox) %>% polishDiversityResults
riaz0m_wilcox_res     <- wilcox.test.list(riaz0m_wilcox) %>% polishDiversityResults
robert0m_wilcox_res   <- wilcox.test.list(robert0m_wilcox) %>% polishDiversityResults

yusko_mnc_IpNi0m_wilcox_res  <- wilcox.test.list(df = yusko_mnc_IpNi0m_wilcox) %>% polishDiversityResults
yusko_mnc_NiIp0m_wilcox_res  <- wilcox.test.list(yusko_mnc_NiIp0m_wilcox) %>% polishDiversityResults

yusko_cd8_IpNi0m_wilcox_res  <- wilcox.test.list(yusko_cd8_IpNi0m_wilcox) %>% polishDiversityResults
yusko_cd8_NiIp0m_wilcox_res  <- wilcox.test.list(yusko_cd8_NiIp0m_wilcox) %>% polishDiversityResults

write.table(helsinki0m_wilcox_res, "results/diversity/wilcox/helsinki_bl.txt", sep = "\t", quote = F, row.names = F)
write.table(tumeh0m_wilcox_res, "results/diversity/wilcox/tumeh_bl.txt", sep = "\t", quote = F, row.names = F)
write.table(riaz0m_wilcox_res, "results/diversity/wilcox/riaz_bl.txt", sep = "\t", quote = F, row.names = F)
write.table(robert0m_wilcox_res, "results/diversity/wilcox/robert_bl.txt", sep = "\t", quote = F, row.names = F)

write.table(yusko_mnc_IpNi0m_wilcox_res, "results/diversity/wilcox/yusko_mnc_IpNi_bl.txt", sep = "\t", quote = F, row.names = F)
write.table(yusko_mnc_NiIp0m_wilcox_res, "results/diversity/wilcox/yusko_mnc_NiIp_bl.txt", sep = "\t", quote = F, row.names = F)

write.table(yusko_cd8_IpNi0m_wilcox_res, "results/diversity/wilcox/yusko_cd8_IpNi_bl.txt", sep = "\t", quote = F, row.names = F)
write.table(yusko_cd8_NiIp0m_wilcox_res, "results/diversity/wilcox/yusko_cd8_NiIp_bl.txt", sep = "\t", quote = F, row.names = F)


## Combine pd1
pd10m_wilcox = rbind(tumeh0m_wilcox, riaz0m_wilcox)
pd10m_wilcox_res = wilcox.test.list(pd10m_wilcox) %>% polishDiversityResults()
write.table(pd10m_wilcox_res, "results/diversity/wilcox/pd1_bl.txt", sep = "\t", quote = F, row.names = F)


## The only significant findings from baseline
fancy_boxplot(yusko_mnc_IpNi, "inverseSimpsonIndex_mean")
ggsave("results/diversity/plots/yusko_mnc_IpNi_inverseSimpsonIndex_mean.pdf", width = 4, height = 3)

fancy_boxplot(yusko_mnc_IpNi, "shannonWienerIndex_mean")
ggsave("results/diversity/plots/yusko_mnc_IpNi_shannonWienerIndex_mean.pdf", width = 4, height = 3)

fancy_boxplot(rbind(tumeh_diversity, filter(riaz_diversity, overall != "NA")), "normalizedShannonWienerIndex_mean") + facet_wrap(~io_stat)
ggsave("results/diversity/plots/pd1_normalizedShannonWienerIndex_mean.pdf", width = 6, height = 3)

              
