

# === plot_tcrgpBaseline

## Plot boxplot to see differences between R and N
plot_box <- function(x, target = target, x_axis = 'overall'){
  
  ggplot(x, aes_string(x_axis, 'freq', fill = x_axis)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(size = 0.5) + 
    facet_wrap(~target, scales = "free") + facets_nice +
    
    ggsignif::geom_signif(comparisons = list(c("N", "R")), margin_top = 0.01,  test = "wilcox.test", textsize = 3) +
    theme(legend.position = "none") + 
    response_fill + 
    labs(x = "", y = "fraction of cells") 
  
  }

# === plot_tcrgpEvolution
plot_line <- function(tot_df, med_df){
  
  ggplot() +
    # geom_point(data = tot_df, aes(timepoint.x,freq,group=name,color=overall), alpha = 0.3) + 
    geom_path(data = tot_df, aes(timepoint.x,freq,group=name,color=overall), alpha = 0.3) +
    
    # geom_point(data = med_df, aes(timepoint.x,freq,group=overall,color=overall), lwd = 2) + 
    geom_path(data = med_df, aes(timepoint.x,freq,group=overall,color=overall), lwd = 2) +
    
    facet_wrap(~target,ncol=4,scales="free") + response_col + labs(x = "", y = "fraction of cells")
}


## Analyze the effect of drugs on melanoma recognising antigens, only on responding patients
plot_line_drug <- function(tot_df, med_df){
  
  nCol <- tot_df$regimen %>% unique %>% length
  
  ggplot() +
    geom_path(data = tot_df, aes(timepoint.x,freq,group=name,color=regimen), alpha = 0.3) +
    geom_path(data = med_df, aes(timepoint.x,freq,group=regimen,color=regimen), lwd = 2, alpha = 0.8) +
    
    facet_wrap(~target,ncol=4,scales="free") + labs(x = "", y = "fraction of cells") +
    scale_color_manual(values = brewer.pal(nCol, "Set1"))
}



## === plot_tcrgpExpanded
## How many of the expanded clonotypes from pooled responders and non-responders are against which target / species / epitope
plot_bar_prop <- function(x){
  
  ggplot(x, aes(overall,prop,fill=target)) + 
    geom_bar(stat = "identity") + 
    scale_fill_manual(values = c(brewer.pal(3, "Set1"), "grey")) + 
    coord_flip() + labs(x = "", y = "frequency of cells")
  
  }

bar_prop <- function(x){
  
  x %>% distinct(clonotypename, .keep_all = TRUE) %>% 
    group_by(overall,target) %>% summarise(n = n()) %>% mutate(prop = n/sum(n)) 
  
}

massage_pie <- function(x,y){
  
  x %>% 
    distinct(clonotypename, .keep_all = TRUE) %>% 
    group_by(name,target) %>% summarise(n = n()) %>% mutate(prop = n/sum(n), total = sum(n)) %>% 
    left_join(y, by = "name")
  
}

plot_pie <- function(x,nCol = 6){
  
  # x <- x %>% mutate(name2 = paste(name, "\nExpanded:", total), name2 = reorder(name2, total))
  # ggplot(x, aes(x = "",prop,fill=target)) + geom_bar(stat = "identity") + scale_fill_manual(values = c(brewer.pal(3, "Set1"), "grey")) + coord_flip() + labs(x = "", y = "frequency of cells") +
  # coord_polar("y", start=0) +
  # facet_wrap(reorder(name2, desc(total))~overall,ncol=nCol) + theme_void() + guides(colour = guide_legend(override.aes = list(size = 10)))
  
  # x = helsinki_pred_expanded %>% filter(type == "Blood") %>% massage_pie(y = helsinki_clin)
  x <- x %>% mutate(name2 = paste(name, io_stat, "\nExpanded:", total), name2 = reorder(name2, total)) 
  
  x %>%
    
    ggplot(aes(x = "", y = prop, fill = target)) +
      geom_bar(width = 1, size = 1, color = "white", stat = "identity") +
      coord_polar("y") +
      ggrepel::geom_text_repel(aes(label = paste0(target, ": ", round(prop, 3)*100, "%")), 
                               position = position_stack(vjust = 0.5)) +
      labs(x = "", y = "frequency of cells")  +  
      scale_fill_manual(values = c(brewer.pal(3, "Set1"), "grey")) +
      # scale_fill_manual(values = c("salmon", "#bcbcbc", "#ffd700")) +
      facet_wrap(reorder(name2, desc(total))~overall,ncol=nCol) + theme_void() + theme(legend.position = "none")  
    
}


## How much each patients expansion is explained by melanoma targeting clonotypes
massage_rank <- function(x,y,target_temp="anti-melanoma"){
  
  x %>% distinct(clonotypename, .keep_all = TRUE) %>% 
    mutate(name = as.character(name)) %>% 
    group_by(name,target, .drop = FALSE) %>% 
    summarise(n = n()) %>% mutate(prop = n/sum(n)) %>% 
    left_join(y, by = "name") %>% 
    filter(target == target_temp)  
}

plot_rank <- function(x){ggplot(x, aes(reorder(name,prop),prop,fill=overall)) + geom_bar(stat = "identity") + response_fill + labs(x = "", y = "frequency of expanding clonotypes") + theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))}


## How much each patients expansion is explained by melanoma targeting clonotypes by epitope
plot_rank2       <- function(x){
  
  plot_rank1 <- function(x){ ggplot(x, aes(reorder(name,prop),prop,fill=overall)) + geom_bar(stat = "identity") + response_fill + labs(x = "", y = "frequency of expanding clonotypes", title = unique(x$pred_epitope)) + theme_classic() + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))}
  
  x <- x %>% mutate(name = paste0(overall, "_", name))
  dfs <- x %>% split(x$pred_epitope)
  
  plot_list <- lapply(dfs, plot_rank1)
  do.call(gridExtra::grid.arrange, c(plot_list, nrow = 1))
  
}

epitope_expanded <- function(x, y){
  
  z <- x %>% 
    distinct(clonotypename, .keep_all = TRUE) %>% 
    mutate(name = as.character(name)) %>% 
    group_by(name, target, pred_epitope) %>% 
    summarise(n = n())
  
  z <- aggregate(n ~ name + pred_epitope, z, FUN = sum, drop = F) %>% 
    tidyr::replace_na(list(n = 0)) %>% 
    mutate(prop = n/sum(n)) %>% 
    filter(pred_epitope %in% melanoma_epitopes) %>% 
    left_join(y, by = "name")
  
  return(z)
  
}
