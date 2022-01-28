

# Stylized Boxplot
boxplot =  geom_boxplot(outlier.colour = NULL, aes_string(colour="SampleClass", fill="SampleClass")) + # geom_boxplot(notch=T) to compare groups
  stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })
# Different axis options. If you are using facets and a legend the x-axis is redundant.

#No X Axis
theme = theme_update(axis.text.x=element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), axis.title.x=element_blank())
#No Y Axis
theme = theme_update(axis.text.y=element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(), axis.title.y=element_blank())
#No Y Axis Label + Grey Axis Numbers
theme = theme_update(axis.line.y = element_blank(), axis.title.y=element_blank(), axis.text.y = element_text(colour="grey"), axis.ticks.y= element_line(colour="grey"))
# Allows multiple scales and altered grouping. Try them out.

#Different Scale Per Facet
boxplot = boxplot + facet_wrap(~ Gland, nrow = 1, scales="free")
#Same Scale Per Facet
boxplot = boxplot + facet_grid(facets = ". ~ Gland")



## Tumor turnover - how many of the clonotypes can be seen on the second time point
fetchNonOnverlapping <- function(folder){

  files <- list.files(folder, pattern = "table.collapsed.txt")
  non_overlap <- NULL

  for(i in 1:length(files)){

    file = files[i]
    print(file)

    df <- fread(paste0(folder, file))

    if(nrow(df) == 1){

      df = data.frame("baseline" = 0, "followup" = 0, "filename" = substr(file, 1,17))
      non_overlap[[i]] <- df
    }

    else{
      df[nrow(df),15] = sum(df[-nrow(df),15])
      df[nrow(df),16] = sum(df[-nrow(df),16])

      df =  df %>% filter(cdr3aa == "NonOverlapping")
      df <- df[,15:16]

      colnames(df) <- c("baseline", "followup")
      non_overlap[[i]] <- df %>% mutate(filename = substr(file, 1,17))
      }
  }

  non_overlap <- do.call(rbind, non_overlap)
  return(non_overlap)

}


plotOverlapPath <- function(df){

  melt(select(df, baseline,followup,name,overall), id = c("name","overall")) %>%

    filter(variable %in% c("baseline", "followup")) %>%
    ggplot(aes(variable, value, group=name, color = overall)) + geom_path() + labs(x = "", y = "freq of shared cells") +
    scale_color_manual(values = c("lightgrey", "salmon"))

}



tumeh_overlap <- fetchNonOnverlapping("results/turnover/Tumeh/files/")
tumeh_overlap <- tumeh_overlap %>% bind_cols(breakName(tumeh_overlap$filename))

tumeh_overlap %>% plotOverlapPath
ggsave("results/turnover/plots/path_tumeh.pdf", width = 6, height = 4)

tumeh_overlap %>% plotOverlapBox
ggsave("results/turnover/plots/box_tumeh.pdf", width = 6, height = 4)

robert_overlap <- fetchNonOnverlapping("results/turnover/robert/files/")
robert_overlap <- robert_overlap %>% bind_cols(breakName(robert_overlap$filename))

robert_overlap %>% plotOverlapPath
ggsave("results/turnover/plots/path_robert.pdf", width = 6, height = 4)

robert_overlap %>% plotOverlapBox
ggsave("results/turnover/plots/box_robert.pdf", width = 6, height = 4)



yusko_mnc_overlap <- fetchNonOnverlapping("results/turnover/yusko/mnc/files/")
yusko_mnc_overlap <- yusko_mnc_overlap %>% bind_cols(breakName(yusko_mnc_overlap$filename))

yusko_mnc_overlap %>% filter(regimen == "antiCTLA4+antiPD1") %>% plotOverlapPath
ggsave("results/turnover/plots/path_yusko_mnc_ipi_nivo.pdf", width = 6, height = 4)

yusko_mnc_overlap %>% filter(regimen == "antiPD1+antiCTLA4") %>% plotOverlapPath
ggsave("results/turnover/plots/path_yusko_mnc_nivo_ipi.pdf", width = 6, height = 4)

yusko_mnc_overlap %>% filter(regimen == "antiCTLA4+antiPD1") %>% plotOverlapBox
ggsave("results/turnover/plots/box_yusko_mnc_ipi_nivo.pdf", width = 6, height = 4)

yusko_mnc_overlap %>% filter(regimen == "antiPD1+antiCTLA4") %>% plotOverlapBox
ggsave("results/turnover/plots/box_yusko_mnc_nivo_ipi.pdf", width = 6, height = 4)



yusko_cd8_overlap <- fetchNonOnverlapping("results/turnover/yusko/cd8/files/")
yusko_cd8_overlap <- yusko_cd8_overlap %>% bind_cols(breakName(yusko_cd8_overlap$filename))

yusko_cd8_overlap %>% filter(regimen == "antiCTLA4+antiPD1") %>% plotOverlapPath
ggsave("results/turnover/plots/path_yusko_cd8_ipi_nivo.pdf", width = 6, height = 4)

yusko_cd8_overlap %>% filter(regimen == "antiPD1+antiCTLA4") %>% plotOverlapPath
ggsave("results/turnover/plots/path_yusko_cd8_nivo_ipi.pdf", width = 6, height = 4)

yusko_cd8_overlap %>% filter(regimen == "antiCTLA4+antiPD1") %>% plotOverlapBox
ggsave("results/turnover/plots/box_yusko_cd8_ipi_nivo.pdf", width = 6, height = 4)

yusko_cd8_overlap %>% filter(regimen == "antiPD1+antiCTLA4") %>% plotOverlapBox
ggsave("results/turnover/plots/box_yusko_cd8_nivo_ipi.pdf", width = 6, height = 4)
