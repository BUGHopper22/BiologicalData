#!/bin/bash
echo "Start"
echo "insert folder"

read folderName

inputDataset="../init_input/kinaseDataset.fasta"
inputFolder="../results/R05_hmmerTest/$folderName"
outputFolder="../results/R06_hmmerTestCompareDataset/$folderName/"

mkdir -p $outputFolder

count=0

for file in "$inputFolder"/*
do
	outputTitleTemp=${file%.hmm}
	outputTitle=${outputTitleTemp##*/}
	outputTblout="$outputFolder$outputTitle.tblout"
	outputAli="$outputFolder$outputTitle.ali"

	hmmsearch "--tblout" $outputTblout $file $inputDataset > $outputAli
	let "count+=1"
	echo "file analyzed: $count"
done

echo "finish"

