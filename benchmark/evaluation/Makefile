all: calLabel calPurity getRepresentativeList analysisPurity checkTaxonomyStatus

calLabel: src/calLabel.cpp src/groundTruth.cpp
	g++ -O3 ./src/calLabel.cpp ./src/groundTruth.cpp -o calLabel
calPurity: src/calPurity.cpp src/groundTruth.cpp
	g++ -O3 ./src/calPurity.cpp ./src/groundTruth.cpp -o calPurity
getRepresentativeList: src/getRepresentativeList.cpp
	g++ -O3 ./src/getRepresentativeList.cpp -o getRepresentativeList
analysisPurity: src/analysisPurity.cpp
	g++ -O3 ./src/analysisPurity.cpp -o analysisPurity
checkTaxonomyStatus: src/checkTaxonomyStatus.cpp
	g++ -O3 ./src/checkTaxonomyStatus.cpp -o checkTaxonomyStatus

clean:
	rm calLabel calPurity getRepresentativeList analysisPurity checkTaxonomyStatus

