#!/bin/bash
echo "Start"
echo "inserisci il nome della cartella dentro ./results/R05hmmerTest/folderName in cui sono presenti i files da analizzare, esempio di input: 2019-05-28_001"
read folderName;

inputPathFolderName="../results/R05_hmmerTest/$folderName"
echo $inputPathFolderName

outputPathName="../results/R06_hmmerTestCompareDataset/$folderName/"
echo $outputPathName

mkdir -p $outputPathName

for file in "$inputPathFolderName"/*
do
	outputTitleTemp=${file%.hmm}
	outputTitle=${outputTitleTemp##*/}

    #tblout output
	outputFullPath1="$outputPathName$outputTitle"

    #.ali output
	outputFullPath2="$outputPathName$outputTitle"

    dataset="../init_input/KinaseDataset.fasta"
	echo $file
	echo $outputFullPath
	hmmsearch "–tblout" $outputFullPath1 $file $dataset > $outputFullPath2
done
echo "Press any key to close"
