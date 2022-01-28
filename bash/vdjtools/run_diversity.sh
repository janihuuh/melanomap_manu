## Run diversity calculations with vdjtools

## Calc diversity stats for each sample
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $helsinki $output_path/diversity/calculated/helsinki
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $tumeh $output_path/diversity/calculated/tumeh
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $riaz $output_path/diversity/calculated/riaz
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $yusko_mnc $output_path/diversity/calculated/yusko_mnc
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $yusko_cd8 $output_path/diversity/calculated/yusko_cd8
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $robert $output_path/diversity/calculated/robert
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $healthy_cd8 $output_path/diversity/calculated/healthy_cd8
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $pruessmann $output_path/diversity/calculated/pruessmann




## Calculate diversity on ACT samples
cd /Users/hru/Dropbox/MelanoMAP/data/unselected_TCRseq/ACT/On_therapy/PBMC
files=$(find "$PWD" -type f) | grep .txt | uniq)
# files=$(find "$PWD" -type f) | egrep "txt|tsv" | uniq


java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $files $output_path/act_on_therapy

java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $cd8_infusion $output_path/cd8_infusion
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $CD8_btlapos_infusion $output_path/CD8_btlapos_infusion
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa $CD8_btlaneg_infusion $output_path/CD8_btlaneg_infusion


## Calc rarefraction plots
java -Xmx4G -jar $vdj RarefactionPlot --wide-plot $helsinki $output_path/helsinki
java -Xmx4G -jar $vdj RarefactionPlot --wide-plot $tumeh $output_path/tumeh
java -Xmx4G -jar $vdj RarefactionPlot --wide-plot $yusko_mnc $output_path/yusko_mnc
java -Xmx4G -jar $vdj RarefactionPlot --wide-plot $yusko_cd8 $output_path/yusko_cd8
java -Xmx4G -jar $vdj RarefactionPlot --wide-plot $robert $output_path/robert
java -Xmx4G -jar $vdj RarefactionPlot --wide-plot $healthy_cd8 $output_path/healthy_cd8
java -Xmx4G -jar $vdj RarefactionPlot --wide-plot $riaz $output_path/riaz


## Calc diversity on resampled samples;
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 10000 $helsinki $tumeh $yusko_mnc $yusko_cd8 $robert $riaz $healthy_cd8 $emerson  $output_path/calculated/total_10k
java -Xmx4G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 20000 $helsinki $tumeh $yusko_mnc $yusko_cd8 $robert $riaz $healthy_cd8 $emerson  $output_path/calculated/total_20k
java -Xmx32G -jar $vdj CalcDiversityStats --intersect-type aa --downsample-to 40000 $helsinki $tumeh $yusko_mnc $yusko_cd8 $robert $riaz $healthy_cd8 $emerson  $output_path/calculated/total_40k

###############
cd $files_path
###############
clear
