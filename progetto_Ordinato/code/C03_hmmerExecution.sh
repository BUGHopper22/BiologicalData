#!/bin/bash
echo "Start"
echo "inserisci il nome della cartella dentro ./result/folderName in cui sono presenti i files da analizzare"
read folderName

inputPathFolderName="../results/R04_msaCleanedTest/$folderName"

outputPathName="../results/R05_hmmerTest/$folderName/"

mkdir -p $outputPathName

count=0

for file in "$inputPathFolderName"/*
do
	outputTitle=${file%.clw}
	outputTitle2=${outputTitle##*/}
	outputFullPath="$outputPathName$outputTitle2.hmm"
	hmmbuild $outputFullPath $file
	count+=1
	echo "file analyzed: $count"
done
echo "Press any key to close"
