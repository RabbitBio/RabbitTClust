all: mapGenome calGenome calLabel calViral calPurity getRepresentativeList

mapGenome: src/mapGenome.cpp
	g++ -O3 ./src/mapGenome.cpp -lz -o mapGenome
calGenome: src/calGenome.cpp
	g++ -O3 ./src/calGenome.cpp -lz -o calGenome
calLabel: src/calLabel.cpp
	g++ -O3 ./src/calLabel.cpp -g -o calLabel
calViral: src/calViral.cpp
	g++ -O3 ./src/calViral.cpp -o calViral
calPurity: src/calViral.cpp
	g++ -O3 ./src/calPurity.cpp -o calPurity
getRepresentativeList: src/getRepresentativeList.cpp
	g++ -O3 ./src/getRepresentativeList.cpp -o getRepresentativeList

clean:
	rm mapGenome calGenome calLabel calViral calPurity getRepresentativeList

