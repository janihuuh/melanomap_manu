
library(dplyr)
library(ggplot2)
library(reshape)
library(gridExtra)
library(RColorBrewer)
library(data.table)

me=system("whoami", intern = T)
setwd(paste0("/Users/", me, "/Dropbox/MelanoMAP/"))

## Plotting
theme_set(theme_bw())
response_col  <- scale_color_manual(values = c("lightgrey", "salmon"))
response_fill <- scale_fill_manual(values = c("lightgrey", "salmon")) 
iostat_fill   <- scale_fill_manual(values = brewer.pal(2, 'Set1'))
iostat_col    <- scale_color_manual(values = brewer.pal(2, 'Set1'))

facets_nice <- theme(strip.background = element_rect(fill="grey96"), strip.text = element_text(colour = 'black')) 

add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))
getPalette6 <- colorRampPalette(brewer.pal(8, "Dark2"))
getPalette7 <- colorRampPalette(brewer.pal(8, "Spectral"))
getPalette8 <- colorRampPalette(brewer.pal(8, "Accent"))

nice_colors <- function(n){scale_color_manual(values = brewer.pal(n, "Set1"))}
nice_fill   <- function(n){scale_fill_manual(values  = brewer.pal(n, "Set1"))}
corner      <- function(df){df[1:5,1:5]}

## Current melanoma epitopes
# melanoma_epitopes         <- c("melana_cdr3b", "meloe1_cdr3b", "melana_cdr3b, meloe1_cdr3b", "meloe1_cdr3b, melana_cdr3b")
  melanoma_epitopes         <- c("melana_cdr3b", "meloe1_cdr3b", 'ELAGIGILTV_cdr3b_comb', 'mart1_cdr3b', 'AMFWSVPTV_cdr3b', 'FLYNLLTRV_cdr3b')

interesting_species       <- c("CMV", "EBV", "InfluenzaA", "Melanoma", "MELANOMA", "HomoSapiens")
interesting_species_tcgrp <- c("CMV", "EBV", "INFa", "MELANOMA", 'ELAGIGILTV', 'mart1')
interesting_species_gliph <- c("CMV", "EBV", "InfluenzaA", "HomoSapiens", "MELANOMA")

melanoma_antigens_vdjdb   <- c("MLANA", "TKT", "SEC24A", "AKAP13", "CTAG1B", "EXOC8","PABPC1", "NDC1", "WT1", "MelanA", "Meloe1", "Mart1")


## Run function-files
# detach("package:org.Hs.eg.db", unload=TRUE)

for(file in list.files("src/jani/R/tcrb/", recursive = T, full.names = T, pattern = "fun")){
  message(file)
  source(file)
}

for(file in list.files("src/jani/R/scrnaseq/", recursive = T, full.names = T, pattern = "fun")){
  message(file)
  source(file)
}

## Load clinical data
source("src/jani/R/tcrb/general/load_clinical.R")

## Open multiple files
# files <- list.files("src/jani/R/tcrb/gliph/", full.names = TRUE)
# file.edit(files)



add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

getPalette  <- colorRampPalette(brewer.pal(9, "Set1"))
getPalette2 <- colorRampPalette(brewer.pal(8, "Set2"))
getPalette3 <- colorRampPalette(brewer.pal(9, "Set3"))
getPalette4 <- colorRampPalette(brewer.pal(9, "Pastel1"))
getPalette5 <- colorRampPalette(brewer.pal(8, "Pastel2"))
getPaletteRdBu <- colorRampPalette(brewer.pal(8, "RdBu"))
getPalette_uchi <- colorRampPalette(c(ggsci::pal_uchicago(palette = "light")(9)))
getPalette_uchi2 <- colorRampPalette(c(ggsci::pal_uchicago(palette = "default")(9)))
getPalette_jco   <- colorRampPalette(c(ggsci::pal_jco()(10)))


facets_nice2 <- theme(strip.text.x = element_text(size=15, angle=0, hjust = 0), strip.background = element_rect(colour="white", fill="white"))
facets_nice3 <- theme(strip.text.x = element_text(size=8, angle=0, hjust = 0), strip.background = element_rect(colour="white", fill="white"))
facets_nice4 <- theme(strip.text.x = element_text(size=11, angle=0, hjust = 0), strip.background = element_rect(colour="white", fill="white"))

add_block <- theme(panel.border = element_rect(colour = "black", fill=NA, size=1), axis.line=element_blank())
