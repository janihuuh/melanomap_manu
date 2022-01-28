


## Basic dotplot
plot_dot <- function(x){
  
  ggplot() +
    geom_point(data = x, aes(Sample1_count, Sample2_count, color = BH.sigf), alpha = 0.3, color = "lightgrey") +
    geom_point(data = subset(x, BH.sigf == "Sigf" & abs(log2_FC_count) > 1), aes(Sample1_count, Sample2_count, color = direction), alpha = 0.5) +
    scale_color_manual(values = c("salmon", "dodgerblue")) + 
    scale_x_log10() + scale_y_log10() +
    labs(color = expression("p"["adj"]~"< 0.05"), x = expression("Time point"["1"]~"count"), y = expression("Time point"["2"]~"count"))
  
}



## Plot the amount of expanded clonotypes per patient as a rank plot
plot_rank <- function(x){
  x <- x %>% mutate(name = paste(overall, "_", name))
  ggplot(x, aes(reorder(name, n), n, fill=overall)) + geom_bar(stat = "identity") + response_fill + theme_classic() + 
    labs(x = "", y = "expanded clonotypes") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


## Basic boxplot
plot_box <- function(x){
  
  ggplot(x, aes(overall,n,color=overall,fill=overall)) + geom_boxplot(outlier.shape = NA) + 
    stat_summary(geom = "crossbar", width=0.65, fatten=0, color="white", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
    geom_jitter(size = 0.5, color = "black") + 
    ggsignif::geom_signif(comparisons = list(c("N", "R")), color = "black") +
    labs(x = "", y = "expanded clonotypes") +
    response_col +
    response_fill 

  }

