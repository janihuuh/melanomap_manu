## Query against vdjdb. We use our own vdjdb, where we have included our epitope-specific data as well
## Use confidence score 0, then filter in downstream analyses if need to

## Initialize
files_path=/Users/hru/Dropbox/MelanoMAP
output_path=/Users/hru/Dropbox/MelanoMAP/results/vdjdb

# vdjdb_path=/Users/hru/Documents/Laaketieteen_tohtori/Applications/vdjdb-1.1.5/vdjdb-1.1.5.jar
vdjdb_path=/Users/hru/Dropbox/MelanoMAP/applications/vdjdb-1.1.5/vdjdb-1.1.5.jar
vdjdb_new=/Users/hru/Dropbox/MelanoMAP/data/selected_TCRseq/vdjdb_new

cd $files_path/data/unselected_TCRseq/Helsinki/;
helsinki=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Tumeh/;
tumeh=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Yusko/MNC_corrected/;
yusko_mnc_baseline=$(ls -d "$PWD"/*0m*);

cd $files_path/data/unselected_TCRseq/Yusko/CD8/;
yusko_cd8=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Robert/;
robert=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Riaz/;
riaz=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Healthy/Helsinki_CD8;
healthy_cd8=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Healthy/Emerson_subsample;
healthy=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/TCGA/vdjtools/beta/;
tcga_beta=$(ls -d "$PWD"/*);



## ACT
cd  $files_path/data/unselected_TCRseq/ACT/On_therapy/PBMC/baseline/
act_baseline=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/ACT/Infusion/CD8_infusion
act_infusion_cd8=$(ls -d "$PWD"/*);


## Run vdjdb

## Unselected TCRb
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $helsinki $output_path/Helsinki/
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $tumeh $output_path/Tumeh/
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $yusko_mnc $output_path/Yusko/MNC/
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $yusko_cd8 $output_path/Yusko/CD8/
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $robert $output_path/Robert/
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $riaz $output_path/Riaz/
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $healthy_cd8 $output_path/healthy_cd8/
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $healthy $output_path/healthy/

java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $tcga_beta $output_path/tcga/beta/


java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $act_baseline $output_path/act_baseline/
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 $act_infusion_cd8 $output_path/act_infusion_cd8/


## Expanded TCRb
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 /Users/hru/Dropbox/MelanoMAP/results/expansion/expanded/helsinki_vdjdb.txt $output_path/expanded/Helsinki_expanded
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 /Users/hru/Dropbox/MelanoMAP/results/expansion/expanded/tumeh_vdjdb.txt $output_path/expanded/Tumeh_expanded
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 /Users/hru/Dropbox/MelanoMAP/results/expansion/expanded/yusko_mnc_vdjdb.txt $output_path/expanded/Yusko_MNC_expanded
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 /Users/hru/Dropbox/MelanoMAP/results/expansion/expanded/yusko_cd8_vdjdb.txt $output_path/expanded/Yusko_CD8_expanded
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0 /Users/hru/Dropbox/MelanoMAP/results/expansion/expanded/robert_vdjdb.txt $output_path/expanded/Robert_expanded

## scRNAseq
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0  /Users/hru/Dropbox/MelanoMAP/data/scRNAseq+TCRseq/sf_vdjdb.txt $output_path/sf_tcrb
java -Xmx4G -jar $vdjdb_path --database $vdjdb_new -R TRB -S human --vdjdb-conf 0  /Users/hru/Dropbox/MelanoMAP/data/scRNAseq+TCRseq/li_vdjdb.txt $output_path/li_tcrb
