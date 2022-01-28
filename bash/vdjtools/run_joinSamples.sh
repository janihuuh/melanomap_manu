## Combine all the samples together, like RNAseq data

java -Xmx16G -jar $vdj JoinSamples --intersect-type strict --times-detected 1 -compress $helsinki $output_path/pooled/helsinki
java -Xmx4G -jar  $vdj JoinSamples --intersect-type strict --times-detected 1 -compress $tumeh $output_path/pooled/tumeh
java -Xmx4G -jar  $vdj JoinSamples --intersect-type strict --times-detected 1 -compress $riaz $output_path/pooled/riaz
java -Xmx16G -jar $vdj JoinSamples --intersect-type strict --times-detected 1 -compress $yusko_mnc $output_path/pooled/yusko_mnc
java -Xmx16G -jar $vdj JoinSamples --intersect-type strict --times-detected 1 -compress $yusko_cd8 $output_path/pooled/yusko_cd8
java -Xmx16G -jar $vdj JoinSamples --intersect-type strict --times-detected 1 -compress $robert $output_path/pooled/robert
