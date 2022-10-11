#!/bin/bash
outputDir="genbankDir"
echo $#
if [ $# -ge 1 ]
then
	outputDir=$1
fi
mkdir -p $outputDir
zcat bact_GenBank.list.gz | while read -r line ; 
do
    fname=$(echo $line | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
    #echo "$line/$fname" ;
    wget -c "$line/$fname" ;
		mv "$fname" $outputDir
done
