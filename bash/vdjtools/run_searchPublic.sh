## Look for public clonotypes, shared by matching aa


## Initialize
files_path=/Users/$me/Dropbox/MelanoMAP
output_path=/Users/$me/Dropbox/MelanoMAP/results/public
vdj=/Users/hru/Documents/Laaketieteen_tohtori/Applications/vdjtools-1.1.4/vdjtools-1.1.4.jar

cd $files_path/data/unselected_TCRseq/Helsinki/;
helsinki_baseline=$(ls -d "$PWD"/*PB*0m*);

cd $files_path/data/unselected_TCRseq/Tumeh/;
tumeh_baseline=$(ls -d "$PWD"/*0m*);

cd $files_path/data/unselected_TCRseq/Riaz/;
riaz_baseline=$(ls -d "$PWD"/*0m*);

cd $files_path/data/unselected_TCRseq/Yusko/MNC_corrected/;
yusko_mnc_baseline=$(ls -d "$PWD"/*0m*);

cd $files_path/data/unselected_TCRseq/Yusko/CD8/;
yusko_cd8_baseline=$(ls -d "$PWD"/*0m*);

cd $files_path/data/unselected_TCRseq/Robert/;
robert_baseline=$(ls -d "$PWD"/*0m*);

cd $files_path/data/unselected_TCRseq/Healthy/Helsinki_CD8;
healthy_cd8=$(ls -d "$PWD"/);

baseline_all=$(echo $helsinki_baseline $tumeh_baseline $riaz_baseline $yusko_mnc_baseline $yusko_cd8_baseline $robert_baseline)



## Search for public clonotypes at baseline in vdj-tools
vdjtools JoinSamples $helsinki_baseline $output_path/Helsinki/helsinki_baseline
vdjtools JoinSamples $tumeh_baseline $output_path/Tumeh/tumeh_baseline
vdjtools JoinSamples $riaz_baseline $output_path/Riaz/riaz_baseline
vdjtools JoinSamples $yusko_mnc_baseline $output_path/Yusko/MNC/yusko_mnc_baseline
vdjtools JoinSamples $yusko_cd8_baseline $output_path/Yusko/CD8/yusko_cd8_baseline
vdjtools JoinSamples $robert_baseline $output_path/Robert/roebrt_baseline

## Pool all baseline
java -Xmx16G -jar $vdj JoinSamples $baseline_all $output_path/baseline_all
