

## Only the significantly expanded clonotypes
helsinki_expanded       <- fread("results/expansion/expanded/helsinki_sigf_expanded.txt") 
tumeh_expanded          <- fread("results/expansion/expanded/tumeh_sigf_expanded.txt") 
yusko_mnc_expanded      <- fread("results/expansion/expanded/yusko_mnc_sigf_expanded.txt")
yusko_cd8_expanded      <- fread("results/expansion/expanded/yusko_cd8_sigf_expanded.txt")
robert_expanded         <- fread("results/expansion/expanded/robert_sigf_expanded.txt") 
riaz_expanded           <- fread("results/expansion/expanded/riaz_sigf_expanded.txt") 

helsinki_expanded       <- helsinki_expanded  %>% bind_cols(breakName(helsinki_expanded$sample1_name))  %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))
tumeh_expanded          <- tumeh_expanded     %>% bind_cols(breakName(tumeh_expanded$sample1_name))     %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))
yusko_mnc_expanded      <- yusko_mnc_expanded %>% bind_cols(breakName(yusko_mnc_expanded$sample1_name)) %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))
yusko_cd8_expanded      <- yusko_cd8_expanded %>% bind_cols(breakName(yusko_cd8_expanded$sample1_name)) %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))
robert_expanded         <- robert_expanded    %>% bind_cols(breakName(robert_expanded$sample1_name))    %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))
riaz_expanded           <- riaz_expanded      %>% bind_cols(breakName(riaz_expanded$sample1_name))      %>% mutate(clonotypename = paste(name, cdr3aa, sep = "_"))

## All clonotypes; not needed in any scripts at least now
# helsinki_exp            <- fread("results/expansion/helsinki_expansion.txt")
# tumeh_exp               <- fread("results/expansion/tumeh_expansion.txt")
# yusko_mnc_exp           <- fread("results/expansion/yusko_mnc_expansion.txt")
# yusko_cd8_exp           <- fread("results/expansion/yusko_cd8_expansion.txt")
# robert_exp              <- fread("results/expansion/robert_expansion.txt")
# riaz_expanded           <- fread("results/expansion/riaz_expansion.txt")
