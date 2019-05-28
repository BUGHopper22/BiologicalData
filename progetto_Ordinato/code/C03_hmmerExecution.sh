#!/bin/bash
echo "Start"
echo "inserisci il nome della cartella dentro ./result/folderName in cui sono presenti i files da analizzare"
read folderName

inputPathFolderName="../results/R04_msaCleanedTest/$folderName"
echo $inputPathFolderName

outputPathName="../results/R05_hmmerTest/$folderName/"
echo $outputPathName

mkdir -p $outputPathName

for file in "$inputPathFolderName"/*
do
	
	outputTitle=${file%.clw}
	outputTitle2=${outputTitle##*/}
	outputFullPath="$outputPathName$outputTitle2.hmm"
	echo $file
	echo $outputFullPath
	hmmbuild $outputFullPath $file	
done
echo "Press any key to close"
