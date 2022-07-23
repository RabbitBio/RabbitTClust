all: mapGenome calGenome calLabel calPurity getRepresentativeList

mapGenome: src/mapGenome.cpp
	g++ -O3 ./src/mapGenome.cpp -lz -o mapGenome
calGenome: src/calGenome.cpp
	g++ -O3 ./src/calGenome.cpp -lz -o calGenome
calLabel: src/calLabel.cpp src/groundTruth.cpp
	g++ -O3 ./src/calLabel.cpp ./src/groundTruth.cpp -o calLabel
calPurity: src/calPurity.cpp src/groundTruth.cpp
	g++ -O3 ./src/calPurity.cpp ./src/groundTruth.cpp -o calPurity
getRepresentativeList: src/getRepresentativeList.cpp
	g++ -O3 ./src/getRepresentativeList.cpp -o getRepresentativeList

clean:
	rm mapGenome calGenome calLabel calPurity getRepresentativeList

