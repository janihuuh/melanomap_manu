
## Downsample the TIL MNC cells to 10k
java -Xmx4G -jar $vdj DownSample --size 10000 $yusko_mnc            data/unselected_tcrb/resampled/Yusko_mnc/
java -Xmx4G -jar $vdj DownSample --size 10000 $riaz         data/unselected_tcrb/resampled/Riaz
java -Xmx4G -jar $vdj DownSample --size 10000 $tumeh         data/unselected_tcrb/resampled/Tumeh/
