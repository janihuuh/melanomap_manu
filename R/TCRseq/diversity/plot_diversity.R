

## On ACT-samples
act_ontherapy = read.delim('results/diversity/act_on_therapy.diversity.aa.exact.txt', stringsAsFactors = F) %>% 
  mutate(name = extractName(sample_id),
         date = extractDate(sample_id)) %>% 
  left_join(act_clin_df) %>% arrange(date)


## Follow patient by patient
act_ontherapy %>% filter(overall != 'Unknown') %>% 
  ggplot(aes(date, 1 - normalizedShannonWienerIndex_mean, group = name, color = overall)) + geom_path(lwd = 0.5) + facet_wrap(~overall, scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = '', y = 'clonality') + 
  response_col + facets_nice 
  # geom_hline(yintercept = 0.2, linetype = 3) + 
  # geom_hline(yintercept = 0.4, linetype = 3)
ggsave('results/diversity/evolution/patient.pdf', width = 8, height = 4)

act_ontherapy %>% filter(overall != 'Unknown') %>% 
  ggplot(aes(date, 1 - normalizedShannonWienerIndex_mean, group=name, color = overall)) + geom_path(lwd = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = '', y = 'clonality') +  
  response_col
ggsave('results/diversity/evolution/all_patient.pdf', width = 7, height = 4)


## Show by box plots
act_ontherapy %>% filter(overall != 'Unknown') %>% 
  ggplot(aes(overall, 1 - normalizedShannonWienerIndex_mean, fill = overall)) + geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~date,nrow=1) + geom_jitter(shape = 21, size = 2) + response_fill + 
  ggsignif::geom_signif(comparisons = list(c('N', 'R'))) + labs(x = '', y = 'clonality') +
  facets_nice + theme(legend.position = 'none')
ggsave('results/diversity/evolution/box.pdf', width = 14, height = 4)

## Show evolution of clonality as fold change
df = act_ontherapy %>% filter(overall != 'Unknown') %>% filter(date %in% c('BL', 'wk001-2'))
cast(df, name+overall~date,value = 'normalizedShannonWienerIndex_mean') %>% 
  mutate(`wk001-2` = 1 - `wk001-2`,
         BL = 1 - BL) %>% 
  mutate(fc = `wk001-2` / BL) %>% 
  
  ggplot(aes(overall,fc,fill=overall)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(size = 2, shape = 21) + ggsignif::geom_signif(comparisons = list(c('N', 'R'))) +
    response_fill + labs(x = '', y = 'fold change on clonality') + theme(legend.position = 'none')
ggsave('results/diversity/evolution/fc_clonality_bl_first.pdf', width = 3, height = 4)


## Median effect and 0.75% CI
ci_df = aggregate(normalizedShannonWienerIndex_mean~date+overall, act_ontherapy, gmodels::ci, confidence = 0.75)
md_df = aggregate(normalizedShannonWienerIndex_mean~date+overall, act_ontherapy, median)
df    = cbind(md_df, ci_df$normalizedShannonWienerIndex_mean)

df %>% filter(overall != 'Unknown') %>% 
  
  ggplot(aes(date, 1 - normalizedShannonWienerIndex_mean, group=overall, color = overall)) + 
    geom_path(lwd = 1) + 
    geom_errorbar(aes(ymin=1-`CI lower`, ymax=1-`CI upper`), width=0.3, alpha = 0.4, lwd = 0.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = '', y = 'clonality') +  
    response_col + ylim(c(0,0.45))
ggsave('results/diversity/evolution/median.pdf', width = 6, height = 4)


aggregate(normalizedShannonWienerIndex_mean~date+overall, act_ontherapy, median) %>% 
  filter(overall != 'Unknown') %>% 
  ggplot(aes(date, 1 - normalizedShannonWienerIndex_mean, group=overall, color = overall)) + geom_path(lwd = 1) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = '', y = 'clonality') +  
    response_col



## === On cehckpoint inhibitor samples

plotFoldchange             <- function(foldchange_file){
    
    limit = foldchange_file$foldchange[is.finite(foldchange_file$foldchange)] %>% abs %>% max %>% round(0)
    
    foldchange_file %>%
    ggplot(aes(overall,foldchange, fill = overall)) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = 3, color = 'black') +
    geom_jitter(size = 0.5) +
    ggsignif::geom_signif(comparisons = list(c('N', 'R')), na.rm = T) +
    response_fill + facets_nice + theme(legend.position = 'none') +
    labs(x = "", y = "log2(fc)") +
    ylim(-limit, limit)
    
}

plotFoldchangeTimepoint    <- function(foldchange_file){
    
    df = melt(foldchange_file, id = c("merge_id", "name", "variable", "foldchange", "response", "type", "regimen", "overall", "timepoint", "celltype", "io_stat"))
    colnames(df)[which(colnames(df) == "variable")[2]] = "timepoints"
    levels(df$timepoints) = c("0m", "3m")
    
    df$timepoint
    limit = df$foldchange[is.finite(df$foldchange)] %>% abs %>% max %>% round(0)
    
    df %>%
    ggplot(aes(timepoints,value, fill = timepoints)) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 0, linetype = 3, color = 'black') +
    geom_jitter(size = 0.5) +
    ggsignif::geom_signif(comparisons = list(c('0m', '3m')), na.rm = T) +
    facets_nice + theme(legend.position = 'none') +
    labs(x = "", y = "value") +
    scale_fill_manual(values = brewer.pal(2, 'Pastel2'))
    
}


plotFoldchangeComparisons      <- function(foldchange_file, x_axis = "regimen", calculate_y_pos = F, map_signif_level = F, plotSigf = T){
    
    ## Set upper and lower limits so that the figures are easier to interpret
    limit = foldchange_file$foldchange[is.finite(foldchange_file$foldchange)] %>% abs %>% max %>% round(0)
    
    ## Calculate different comparisons
    if(plotSigf){
        
        regimens    = combn(unique(foldchange_file$regimen), 2) %>% t %>% as.data.frame()
        comparisons = NULL
        for(i in 1:nrow(regimens)){
            comparisons[[i]] <- regimens[i,] %>% sapply(as.character) %>% as.vector()
        }
        
        ## Where the significance lines should be plotted
        y_positions = c(limit-0.2, limit-1.5,limit-3)
        
        if(calculate_y_pos){
            nComparisons  = length(comparisons)
            y_positions   = seq(1, limit, length.out = nComparisons)
        }
        
    }
    
    
    
    ## Actual plot
    foldchange_file %>%
    ggplot(aes_string(x_axis, 'foldchange', fill = x_axis)) +
    geom_hline(yintercept = 0, linetype = 3, color = 'black') +
    
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.5) +
    ggsignif::geom_signif(comparisons = comparisons, na.rm = T, y_position = y_positions, map_signif_level = map_signif_level) +
    facets_nice + theme(legend.position = 'none') +
    labs(x = "", y = "log2(fc)") +
    scale_fill_manual(values = brewer.pal(9, 'Pastel2')) +
    ylim(-limit, limit) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
}





## Combine all tumor samples
tumor_diversity <- rbind(tumeh_diversity, riaz_diversity, yusko_mnc_diversity)

data = tumor_diversity %>% select(reads:cohort) %>% select(-contains("std")) %>% melt(id = "cohort") 

ggplot(tumor_diversity, aes(diversity)) + geom_histogram() 

sort(tumor_diversity$diversity)
summary(tumor_diversity$diversity)

table(tumor_diversity$diversity > 1e3)
table(tumor_diversity$diversity > 500)

samples_to_use <- tumor_diversity %>% filter(diversity < 5e2) %>% pull(sample_id)
names_to_use   <- extractName(samples_to_use) %>% unique 



p <- ggplot(data, aes(cohort,value)) + geom_boxplot() + facet_wrap(~variable, scales = "free")
ggsave(plot = p, "results/diversity/plots/tumor/meta.pdf", width = 24, height = 12)
