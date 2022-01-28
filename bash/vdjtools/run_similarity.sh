## Run diversity calculations with vdjtools
output_path=/Users/hru/Dropbox/MelanoMAP/results/similarity
mkdir $output_path

## Calc diversity stats for each sample
java -Xmx4G -jar $vdj CalcPairwiseDistances --intersect-type aa $helsinki $output_path/calculated/helsinki
java -Xmx4G -jar $vdj CalcPairwiseDistances --intersect-type aa $tumeh $output_path/calculated/tumeh
java -Xmx4G -jar $vdj CalcPairwiseDistances --intersect-type aa $riaz $output_path/calculated/riaz
java -Xmx4G -jar $vdj CalcPairwiseDistances --intersect-type aa $yusko_mnc $output_path/calculated/yusko_mnc
java -Xmx4G -jar $vdj CalcPairwiseDistances --intersect-type aa $yusko_cd8 $output_path/calculated/yusko_cd8
java -Xmx4G -jar $vdj CalcPairwiseDistances --intersect-type aa $robert $output_path/calculated/robert
java -Xmx4G -jar $vdj CalcPairwiseDistances --intersect-type aa $healthy_cd8 $output_path/calculated/healthy_cd8

## Cluster based on the similarity measures
java -Xmx8G -jar $vdj ClusterSamples --plot $output_path/calculated/helsinki $output_path/calculated/clustered_helsinki
java -Xmx8G -jar $vdj ClusterSamples --plot $output_path/calculated/tumeh $output_path/calculated/clustered_tumeh
java -Xmx8G -jar $vdj ClusterSamples --plot $output_path/calculated/riaz $output_path/calculated/clustered_riaz
java -Xmx8G -jar $vdj ClusterSamples --plot $output_path/calculated/yusko_mnc $output_path/calculated/clustered_yusko_mnc
java -Xmx8G -jar $vdj ClusterSamples --plot $output_path/calculated/yusko_cd8 $output_path/calculated/clustered_yusko_cd8
java -Xmx8G -jar $vdj ClusterSamples --plot $output_path/calculated/robert $output_path/calculated/clustered_robert
