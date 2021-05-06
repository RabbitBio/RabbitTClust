TARGET=clust

CXX=g++
#CXXFLAGS=-g -O3 -std=c++14 -fopenmp -D THREADPOOL_MINHASH
CXXFLAGS=-g -O3 -std=c++14 -fopenmp -D RABBIT_IO

INCLUDE=./RabbitSketch/build/include
INCLUDE1=./RabbitIO/build/include
LIBS=./RabbitSketch/build/lib
LIBS1=./RabbitIO/build/io
OBJECTS=main.o MST.o SketchInfo.o
#OBJECTS=demo.o MST.o
#OBJECTS=test.o MST.o

$(TARGET):$(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -I$(INCLUDE) -I$(INCLUDE1) -lz -lrabbitsketch -lio_lib -L$(LIBS) -L$(LIBS1)


#$(OBJECTS):main.cpp MST.cpp
#	$(CXX) $(CXXFLAGS) -c $< -o $@ -I$(INCLUDE) -lz

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -I$(INCLUDE) -I$(INCLUDE1) -lz

clean:
	rm *.o $(TARGET)
