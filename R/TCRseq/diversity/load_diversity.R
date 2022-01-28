
## Calculated with vdjtools

## Exact
helsinki_diversity    <- read.delim('results/diversity/calculated/helsinki.diversity.aa.exact.txt') %>% mutate(cohort = "helsinki")
tumeh_diversity       <- read.delim('results/diversity/calculated/tumeh.diversity.aa.exact.txt') %>% mutate(cohort = "tumeh")
riaz_diversity        <- read.delim('results/diversity/calculated/riaz.diversity.aa.exact.txt') %>% mutate(cohort = "riaz")
robert_diversity      <- read.delim('results/diversity/calculated/robert.diversity.aa.exact.txt') %>% mutate(cohort = "robert")
yusko_mnc_diversity   <- read.delim('results/diversity/calculated/yusko_mnc.diversity.aa.exact.txt') %>% mutate(cohort = "yusko_mnc")
yusko_cd8_diversity   <- read.delim('results/diversity/calculated/yusko_cd8.diversity.aa.exact.txt') %>% mutate(cohort = "yusko_cd8")

## Resampled
# helsinki_diversity  <- read.delim('results/diversity/calculated/helsinki.diversity.aa.resampled.txt')
# tumeh_diversity     <- read.delim('results/diversity/calculated/tumeh.diversity.aa.resampled.txt') 
# riaz_diversity      <- read.delim('results/diversity/calculated/riaz.diversity.aa.resampled.txt') 
# robert_diversity    <- read.delim('results/diversity/calculated/robert.diversity.aa.resampled.txt') 
# yusko_mnc_diversity <- read.delim('results/diversity/calculated/yusko_mnc.diversity.aa.resampled.txt') 
# yusko_cd8_diversity <- read.delim('results/diversity/calculated/yusko_cd8.diversity.aa.resampled.txt') 







## Only baseline
helsinki_diversity_0m  <- helsinki_diversity %>% cbind(breakName(helsinki_diversity$sample_id)) %>% filter(type == "Blood" & timepoint == '0m')
tumeh_diversity_0m     <- tumeh_diversity %>% cbind(breakName(tumeh_diversity$sample_id)) %>% filter( timepoint == '0m')
riaz_diversity_0m      <- riaz_diversity %>% cbind(breakName(riaz_diversity$sample_id)) %>% filter(overall != "NA" & timepoint == '0m')
robert_diversity_0m    <- robert_diversity %>% cbind(breakName(robert_diversity$sample_id)) %>% filter(timepoint == '0m')

yusko_mnc_IpNi_0m      <- yusko_mnc_diversity %>% cbind(breakName(yusko_mnc_diversity$sample_id)) %>% filter(overall != "NA" & regimen == "antiCTLA4+antiPD1" & timepoint == "0m")
yusko_mnc_NiIp_0m      <- yusko_mnc_diversity %>% cbind(breakName(yusko_mnc_diversity$sample_id)) %>% filter(overall != "NA" & regimen == "antiPD1+antiCTLA4" & timepoint == "0m")

yusko_cd8_IpNi_0m      <- yusko_cd8_diversity %>% cbind(breakName(yusko_cd8_diversity$sample_id)) %>% filter(overall != "NA" & regimen == "antiCTLA4+antiPD1" & timepoint == "0m")
yusko_cd8_NiIp_0m      <- yusko_cd8_diversity %>% cbind(breakName(yusko_cd8_diversity$sample_id)) %>% filter(overall != "NA" & regimen == "antiPD1+antiCTLA4" & timepoint == "0m")




## 3m 
helsinki_diversity_3m  <- helsinki_diversity %>% cbind(breakName(helsinki_diversity$sample_id)) %>% filter(type == "Blood" & timepoint == '3m')
tumeh_diversity_3m     <- tumeh_diversity %>% cbind(breakName(tumeh_diversity$sample_id)) %>% filter(timepoint == '3m')
riaz_diversity_3m      <- riaz_diversity %>% cbind(breakName(riaz_diversity$sample_id)) %>% filter(overall != "NA" & timepoint == '3m')
robert_diversity_3m    <- robert_diversity %>% cbind(breakName(robert_diversity$sample_id)) %>% filter(timepoint == '3m')

yusko_mnc_IpNi_3m      <- yusko_mnc_diversity %>% cbind(breakName(yusko_mnc_diversity$sample_id)) %>% filter(overall != "NA" & regimen == "antiCTLA4+antiPD1" & timepoint == "3m")
yusko_mnc_NiIp_3m      <- yusko_mnc_diversity %>% cbind(breakName(yusko_mnc_diversity$sample_id)) %>% filter(overall != "NA" & regimen == "antiPD1+antiCTLA4" & timepoint == "3m")

yusko_cd8_IpNi_3m      <- yusko_cd8_diversity %>% cbind(breakName(yusko_cd8_diversity$sample_id)) %>% filter(overall != "NA" & regimen == "antiCTLA4+antiPD1" & timepoint == "3m")
yusko_cd8_NiIp_3m      <- yusko_cd8_diversity %>% cbind(breakName(yusko_cd8_diversity$sample_id)) %>% filter(overall != "NA" & regimen == "antiPD1+antiCTLA4" & timepoint == "3m")


