all: simulate-longSequence create-containment-bacteria
simulate-longSequence: src/simulate_longSequence.cpp
	g++ -O3 src/simulate_longSequence.cpp -o simulate-longSequence
create-containment-bacteria: src/create_containment_bacteria.cpp
	g++ -O3 src/create_containment_bacteria.cpp -o create-containment-bacteria -lz

clean:
	rm simulate-longSequence create-containment-bacteria
	
