#!/bin/bash
cd /Users/hru/Dropbox/MelanoMAP/

RNKFILE=results/scrnaseq/li/exhausted/anti_melanoma_expression_fc.rnk
GMT=applications/gsea/h.all.v6.2.symbols.gmt # path to gmt file
OUTDIR=results/scrnaseq/li/exhausted/gsea_antimelanoma
LABEL=$RNKFILE

    java -cp /Users/hru/Documents/Laaketieteen_tohtori/Applications/gsea-3.0.jar \
        -Xmx5g xtools.gsea.GseaPreranked \
        -rpt_label $LABEL \
        -rnk $RNKFILE \
        -gmx $GMT \
        -out $OUTDIR \
        -plot_top_x 250 \
        -collapse false \
        -mode Max_probe \
        -norm meandiv \
        -scoring_scheme weighted \
        -include_only_symbols true \
        -make_sets true \
        -rnd_seed 149 \
        -zip_report false \
        -gui false \
        -nperm 1000 \
        -set_min 5 \
        -set_max 500
