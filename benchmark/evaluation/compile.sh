#!/bin/bash
if [ $# -ge 1 ] && [ $1 == "clean" ];
then
	set -x
	rm mapGenome calGenome calResult2 calNMI calF1
else
	set -x
	g++ ./src/mapGenome.cpp -lz -o mapGenome
	g++ ./src/calGenome.cpp -lz -o calGenome
	#g++ ./src/calResult1.cpp -lz -o calResult1
	g++ ./src/calResult2.cpp -lz -o calResult2
	g++ ./src/calNMI.cpp -o calNMI
	g++ ./src/calF1.cpp -o calF1
fi


