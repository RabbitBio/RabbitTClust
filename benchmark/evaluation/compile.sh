#!/bin/bash
if [ $# -ge 1 ] && [ $1 == "clean" ];
then
	set -x
	rm mapGenome calGenome calLabel calViral
else
	set -x
	g++ ./src/mapGenome.cpp -lz -o mapGenome
	g++ ./src/calGenome.cpp -lz -o calGenome
	g++ ./src/calLabel.cpp -o calLabel
	g++ ./src/calViral.cpp -o calViral
fi


