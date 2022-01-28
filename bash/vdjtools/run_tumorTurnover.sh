#!/bin/bash

## Search for TIL "turnover" following therapy
cd /Users/hru/Dropbox/MelanoMAP/src
output_folder=/Users/hru/Dropbox/MelanoMAP/results/turnover
melanomap=/Users/hru/Dropbox/MelanoMAP
vdjtools_path=/Users/hru/Dropbox/MelanoMAP/applications/vdjtools-1.2.1/vdjtools-1.2.1.jar

### More than pair 

## Helsinki
helsinki_folder=$melanomap/data/unselected_TCRseq/Helsinki
helsinki_file=$melanomap/results/expansion/helsinki_sample_matrix.txt

while read line; do
   
    sample1=$(echo $line | awk '{print $1;}')
    sample2=$(echo $line | awk '{print $2;}')
    sample3=$(echo $line | awk '{print $3;}')

    output=$(echo $sample1 | cut -c1-18)

    sample1=$helsinki_folder/$sample1
    sample2=$helsinki_folder/$sample2
    sample3=$helsinki_folder/$sample3

    vdjtools TrackClonotypes --track-sample 0 --plot --top 20 --intersect-type aa $sample1 $sample2 $sample3 $output_folder/helsinki/$output

done < $helsinki_file

## Organize
cd $output_folder/Helsinki/
mkdir plots/
mkdir files

var=$(ls *.pdf)
mv $var plots/

var=$(ls *.txt)
mv $var files/

cd /Users/hru/Dropbox/MelanoMAP/src














### Paired data

## Tumeh
tumeh_folder=$melanomap/data/unselected_TCRseq/Tumeh
tumeh_file=$melanomap/results/expansion/tumeh_sample_pairs.txt

while read line; do
   
    sample1=$(echo $line | awk '{print $1;}')
    sample2=$(echo $line | awk '{print $2;}')
    output=$(echo $sample1 | cut -c1-18)

    sample1=$tumeh_folder/$sample1
    sample2=$tumeh_folder/$sample2

    vdjtools OverlapPair -p $sample1 $sample2 $output_folder/Tumeh/$output

done < $tumeh_file

## Organize
cd $output_folder/Tumeh/
mkdir plots/
mkdir files

var=$(ls *.pdf)
mv $var plots/

var=$(ls *.txt)
mv $var files/

cd /Users/hru/Dropbox/MelanoMAP/src




## yusko_mnc
yusko_mnc_folder=$melanomap/data/unselected_TCRseq/Yusko/MNC_corrected
yusko_mnc_file=$melanomap/results/expansion/yusko_mnc_sample_pairs.txt

while read line; do
   
    sample1=$(echo $line | awk '{print $1;}')
    sample2=$(echo $line | awk '{print $2;}')
    output=$(echo $sample1 | cut -c1-18)

    sample1=$yusko_mnc_folder/$sample1
    sample2=$yusko_mnc_folder/$sample2

    vdjtools OverlapPair -p $sample1 $sample2 $output_folder/Yusko/MNC/$output

done < $yusko_mnc_file

## Organize
cd $output_folder/Yusko/MNC/
mkdir plots/
mkdir files

var=$(ls *.pdf)
mv $var plots/

var=$(ls *.txt)
mv $var files/

cd /Users/hru/Dropbox/MelanoMAP/src




## yusko_cd8
yusko_cd8_folder=$melanomap/data/unselected_TCRseq/Yusko/CD8
yusko_cd8_file=$melanomap/results/expansion/yusko_cd8_sample_pairs.txt

while read line; do
   
    sample1=$(echo $line | awk '{print $1;}')
    sample2=$(echo $line | awk '{print $2;}')
    output=$(echo $sample1 | cut -c1-18)

    sample1=$yusko_cd8_folder/$sample1
    sample2=$yusko_cd8_folder/$sample2

    vdjtools OverlapPair -p $sample1 $sample2 $output_folder/Yusko/CD8/$output

done < $yusko_cd8_file

## Organize
cd $output_folder/Yusko/cd8/
mkdir plots/
mkdir files

var=$(ls *.pdf)
mv $var plots/

var=$(ls *.txt)
mv $var files/

cd /Users/hru/Dropbox/MelanoMAP/src



## Robert
robert_folder=$melanomap/data/unselected_TCRseq/robert
robert_file=$melanomap/results/expansion/robert_sample_pairs.txt

while read line; do
   
    sample1=$(echo $line | awk '{print $1;}')
    sample2=$(echo $line | awk '{print $2;}')
    output=$(echo $sample1 | cut -c1-18)

    sample1=$robert_folder/$sample1
    sample2=$robert_folder/$sample2

    vdjtools OverlapPair -p $sample1 $sample2 $output_folder/robert/$output

done < $robert_file

## Organize
cd $output_folder/robert
mkdir plots/
mkdir files

var=$(ls *.pdf)
mv $var plots/

var=$(ls *.txt)
mv $var files/

cd /Users/hru/Dropbox/MelanoMAP/src

