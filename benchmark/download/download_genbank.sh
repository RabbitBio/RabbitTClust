#!/bin/bash

wget https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
awk -F '\t' 'NR>2 {print $20}' assembly_summary.txt >ftp.list

outputDir="genbankDir"
echo $#
if [ $# -ge 1 ]
then
	outputDir=$1
fi
mkdir -p $outputDir
cat ftp.list | while read line
do
    fname=$(echo $line | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/')
    #echo "$line/$fname"
    wget -c "$line/$fname" ;
		mv "$fname" $outputDir
done
