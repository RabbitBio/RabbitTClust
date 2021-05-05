TARGET=clust

CXX=g++
CXXFLAGS=-g -O3 -std=c++14 -fopenmp -D THREADPOOL_MINHASH

INCLUDE=./RabbitSketch/build/include
LIBS=./RabbitSketch/build/lib
OBJECTS=main.o MST.o SketchInfo.o
#OBJECTS=demo.o MST.o
#OBJECTS=test.o MST.o

$(TARGET):$(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -I$(INCLUDE) -lz -lrabbitsketch -L$(LIBS)


#$(OBJECTS):main.cpp MST.cpp
#	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INCLUDE) -lz

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -I$(INCLUDE) -lz

clean:
	rm *.o $(TARGET)
