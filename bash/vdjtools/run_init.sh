
me=$(whoami)

##########
files_path=/Users/$me/Dropbox/MelanoMAP
vdj=/Users/$me/Dropbox/MelanoMAP/applications/vdjtools-1.2.1/vdjtools-1.2.1.jar
output_path=/Users/$me/Dropbox/MelanoMAP/results/
##########

cd $files_path/data/unselected_TCRseq/Helsinki/;
helsinki=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/new_batch/;
new_batch=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Tumeh/;
tumeh=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Yusko/MNC_corrected/;
yusko_mnc=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Yusko/CD8/;
yusko_cd8=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Robert/;
robert=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Riaz/;
riaz=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Pruessmann/;
pruessmann=$(ls -d "$PWD"/*);



cd $files_path/data/unselected_TCRseq/Healthy/Helsinki_CD8/;
healthy_cd8=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Healthy/Emerson_subsample;
healthy=$(ls -d "$PWD"/*);

## ACT
cd $files_path/data/unselected_TCRseq/ACT/On_therapy/PBMC/wk0
act_baseline=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/Healthy/Helsinki_CD8/;
healthy_cd8=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/ACT/Infusion/CD8_infusion/;
cd8_infusion=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/ACT/Infusion/CD8_btlapos_infusion/;
CD8_btlapos_infusion=$(ls -d "$PWD"/*);

cd $files_path/data/unselected_TCRseq/ACT/Infusion/CD8_btlaneg_infusion/;
CD8_btlaneg_infusion=$(ls -d "$PWD"/*);

cd /Users/janihuuh/Dropbox/Emerson/data/vdjt/;
emerson=$(ls -d "$PWD"/*);

###############
cd $files_path
###############
clear
