
## Track the persisting clonotypes, find in all the time points

## without baseline sample

cd $files_path/data/unselected_TCRseq/ACT/On_therapy/PBMC/

id_166=$(ls -ld $(find .))
id_166=$(echo $id_166 | tr ' ' '\n' | grep 166_);
id_166=$(echo $id_166 | tr ' ' '\n' | grep -v _BL);

id_188=$(ls -ld $(find .))
id_188=$(echo $id_188 | tr ' ' '\n' | grep 188_);
id_188=$(echo $id_188 | tr ' ' '\n' | grep -v _BL);

id_237=$(ls -ld $(find .))
id_237=$(echo $id_237 | tr ' ' '\n' | grep 237_);
id_237=$(echo $id_237 | tr ' ' '\n' | grep -v _BL);

id_335=$(ls -ld $(find .))
id_335=$(echo $id_335 | tr ' ' '\n' | grep 335_);
id_335=$(echo $id_335 | tr ' ' '\n' | grep -v _BL);

id_341=$(ls -ld $(find .))
id_341=$(echo $id_341 | tr ' ' '\n' | grep 341_);
id_341=$(echo $id_341 | tr ' ' '\n' | grep -v _BL);

id_408=$(ls -ld $(find .))
id_408=$(echo $id_408 | tr ' ' '\n' | grep 408_);
id_408=$(echo $id_408 | tr ' ' '\n' | grep -v _BL);

id_420=$(ls -ld $(find .))
id_420=$(echo $id_420 | tr ' ' '\n' | grep 420_);
id_420=$(echo $id_420 | tr ' ' '\n' | grep -v _BL);

id_441=$(ls -ld $(find .))
id_441=$(echo $id_441 | tr ' ' '\n' | grep 441_);
id_441=$(echo $id_441 | tr ' ' '\n' | grep -v _BL);

id_444=$(ls -ld $(find .))
id_444=$(echo $id_444 | tr ' ' '\n' | grep 444_);
id_444=$(echo $id_444 | tr ' ' '\n' | grep -v _BL);

id_518=$(ls -ld $(find .))
id_518=$(echo $id_518 | tr ' ' '\n' | grep 518_);
id_518=$(echo $id_518 | tr ' ' '\n' | grep -v _BL);

id_520=$(ls -ld $(find .))
id_520=$(echo $id_520 | tr ' ' '\n' | grep 520_);
id_520=$(echo $id_520 | tr ' ' '\n' | grep -v _BL);

id_538=$(ls -ld $(find .))
id_538=$(echo $id_538 | tr ' ' '\n' | grep 538_);
id_538=$(echo $id_538 | tr ' ' '\n' | grep -v _BL);


vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_166 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/166
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_188 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/188
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_237 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/237
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_335 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/335
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_408 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/408
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_420 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/420
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_441 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/441
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_444 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/444
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_518 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/518
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_520 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/520
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_538 $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/538



cd $files_path/results/act_tracked_clonotypes/raw/withoutBaseline/
mkdir $files_path/results/act_tracked_clonotypes/plots/withoutBaseline/
var=$(ls *pdf)
mv $var $files_path/results/act_tracked_clonotypes/plots/withoutBaseline/








## With baseline sample

cd $files_path/data/unselected_TCRseq/ACT/On_therapy/PBMC/

id_166=$(ls -ld $(find .))
id_166=$(echo $id_166 | tr ' ' '\n' | grep 166_);

id_188=$(ls -ld $(find .))
id_188=$(echo $id_188 | tr ' ' '\n' | grep 188_);

id_237=$(ls -ld $(find .))
id_237=$(echo $id_237 | tr ' ' '\n' | grep 237_);

id_335=$(ls -ld $(find .))
id_335=$(echo $id_335 | tr ' ' '\n' | grep 335_);

id_341=$(ls -ld $(find .))
id_341=$(echo $id_341 | tr ' ' '\n' | grep 341_);

id_408=$(ls -ld $(find .))
id_408=$(echo $id_408 | tr ' ' '\n' | grep 408_);

id_420=$(ls -ld $(find .))
id_420=$(echo $id_420 | tr ' ' '\n' | grep 420_);

id_441=$(ls -ld $(find .))
id_441=$(echo $id_441 | tr ' ' '\n' | grep 441_);

id_444=$(ls -ld $(find .))
id_444=$(echo $id_444 | tr ' ' '\n' | grep 444_);

id_518=$(ls -ld $(find .))
id_518=$(echo $id_518 | tr ' ' '\n' | grep 518_);

id_520=$(ls -ld $(find .))
id_520=$(echo $id_520 | tr ' ' '\n' | grep 520_);

id_538=$(ls -ld $(find .))
id_538=$(echo $id_538 | tr ' ' '\n' | grep 538_);


vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_166 $files_path/results/act_tracked_clonotypes/raw/withBaseline/166
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_188 $files_path/results/act_tracked_clonotypes/raw/withBaseline/188
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_237 $files_path/results/act_tracked_clonotypes/raw/withBaseline/237
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_335 $files_path/results/act_tracked_clonotypes/raw/withBaseline/335
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_408 $files_path/results/act_tracked_clonotypes/raw/withBaseline/408
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_420 $files_path/results/act_tracked_clonotypes/raw/withBaseline/420
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_441 $files_path/results/act_tracked_clonotypes/raw/withBaseline/441
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_444 $files_path/results/act_tracked_clonotypes/raw/withBaseline/444
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_518 $files_path/results/act_tracked_clonotypes/raw/withBaseline/518
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_520 $files_path/results/act_tracked_clonotypes/raw/withBaseline/520
vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $id_538 $files_path/results/act_tracked_clonotypes/raw/withBaseline/538



cd $files_path/results/act_tracked_clonotypes/raw/withBaseline/
mkdir $files_path/results/act_tracked_clonotypes/plots/withBaseline/
var=$(ls *pdf)
mv $var $files_path/results/act_tracked_clonotypes/plots/withBaseline/
