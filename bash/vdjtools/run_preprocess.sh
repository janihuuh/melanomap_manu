
## Convert
cd $files_path/data/unselected_TCRseq/new_batch/;
java -Xmx4G -jar $vdj Convert -S ImmunoSeqV2 $new_batch new/


## Filter non/functional
cd $files_path/data/unselected_TCRseq/Helsinki/;
java -Xmx4G -jar $vdj FilterNonFunctional $helsinki functional/

cd $files_path/data/unselected_TCRseq/Tumeh/;
java -Xmx4G -jar $vdj FilterNonFunctional $tumeh functional/

cd $files_path/data/unselected_TCRseq/Yusko/MNC_corrected/;
java -Xmx4G -jar $vdj FilterNonFunctional $yusko_mnc functional/

cd $files_path/data/unselected_TCRseq/Yusko/CD8/;
java -Xmx4G -jar $vdj FilterNonFunctional $yusko_cd8 functional/

cd $files_path/data/unselected_TCRseq/Robert/;
java -Xmx4G -jar $vdj FilterNonFunctional $robert functional/

cd $files_path/data/unselected_TCRseq/Riaz/;
java -Xmx4G -jar $vdj FilterNonFunctional $riaz functional/

cd $files_path/data/unselected_TCRseq/Pruessmann/;
java -Xmx4G -jar $vdj FilterNonFunctional $pruessmann functional/

cd $files_path/data/unselected_TCRseq/new_batch/;
java -Xmx4G -jar $vdj FilterNonFunctional $new_batch functional/



cd $files_path/data/unselected_TCRseq/Healthy/Helsinki_CD8/;
java -Xmx4G -jar $vdj FilterNonFunctional $healthy_cd8 functional/

cd $files_path/data/unselected_TCRseq/Healthy/Emerson_subsample;
java -Xmx4G -jar $vdj FilterNonFunctional $healthy functional/

cd $files_path/data/unselected_TCRseq/Healthy/Emerson_subsample;
java -Xmx4G -jar $vdj FilterNonFunctional $healthy functional/

cd /Users/hru/Dropbox/ACT_Helsinki/On_therapy/PBMC/baseline/
java -Xmx4G -jar $vdj FilterNonFunctional $healthy functional/



## Select top n clonotypes
java -Xmx4G -jar $vdj SelectTop -x 20 $act_baseline $output_path/top_clones/act_baseline
