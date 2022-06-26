![RabbitTClust](rabbittclust.png)

RabbitTClust is an efficient sketch-based genome clustering application 
using minimum spanning trees. It enables processing of large-scale 
datasets by combining dimensionality reduction techniques with streaming 
and parallelization on modern multi-core platforms. 

## Installation
RabbitTClust version 1.0 can only support 64-bit Linux Systems.

### Dependancy
* cmake v.3.0 or later
* c++14
* [zlib](https://zlib.net/)

### Compile and install automatically
```bash
git clone --recursive https://github.com/RabbitBio/RabbitTClust.git
cd RabbitTClust
./install.sh

```

### Compile and install manually 
```bash
git clone --recursive https://github.com/RabbitBio/RabbitTClust.git
cd RabbitTClust

#make rabbitSketch library
cd RabbitSketch
mkdir build && cd build 
cmake -DCXXAPI=ON -DCMAKE_INSTALL_PREFIX=. ..
make -j8 && make install
cd ../../

#make rabbitFX library
cd RabbitFX
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make -j8 && make install
cd ../../

mkdir build && cd build
cmake -DUSE_RABBITFX=ON ..
make -j8 && make install
cd ../

```

## Usage
```bash
usage: clust-mst [-h] [-l] [-t] <int> [-d] <double> -F <string> [-o] <string> -i <string>
usage: clust-mst [-h] [-f] [-E] [-d] <double> -i <string> <string> -o <string>
usage: clust-greedy [-h] [-l] [-t] <int> [-d] <double> -F <string> [-o] <string> -i <string>
usage: clust-greedy [-h] [-f] [-d] <double> -i <string> <string> -o <string>
-h         : this help message
-k <int>   : set kmer size, default 21, for both clust-greedy and clust-mst
-s <int>   : set sketch size, default 1000, for both clust-greedy and clust-mst
-l         : cluster for genomes(not sequences), list input. Lines in each <input> specify paths to genome files, one per line. for both clust-greedy and clust-mst
-c <int>   : compute the containment of genomes, set proportion sketchSize = genomeSize/compress, ATTENTION with MinHash function, for both clust-greedy and clust-mst
-t <int>   : set the thread number, default 1, for both clust-greedy and clust-mst
-d <double>: set the threshold of the clusters from the Minimum Spanning Tree and greedy cluster threshold, default 0.05(0.01) for clust-mst(clust-greedy)
-f         : input files are genomeInfo and MST contents(sketch contents) for clust-mst(clust-greedy)
-E         : input files are genomeInfo and sketch contents contents for clust-mst
-F <string>: sketch function, includes MinHash, WMH, OMH, HLL, default MinHash, for both clust-greedy and clust-mst
-o <string>: path of output file, for both clust-greedy and clust-mst
-i <string>: path of input file, ATTENTION with -f and -E option

```

## Example:
```bash
#for genomes clustering, input is a genome file list:
./clust-mst -l -t 48 -i bacteriaList -o bacteria.clust
./clust-greedy -l -t 48 -i bacteriaList -o bacteria.clust

#for genomes clustering, input is a single genome file:
./clust-mst -d 0.05 -t 48 -i bacteria.fna -o bacteria.clust
./clust-greedy -d 0.05 -t 48 -i bacteria.fna -o bacteria.clust

#for redundancy detection with containment, input is a genome file list:
./clust-mst -l -c 10000 -t 48 -i bacteriaList -o bacteria.out
./clust-greedy -l -c 10000 -t 48 -i bacteriaList -o bacteria.out

#for redundancy detection with containment, input is a single genome file:
./clust-mst -c 10000 -t 48 -i bacteria.fna -o bacteria.out
./clust-greedy -c 10000 -t 48 -i bacteria.fna -o bacteria.out

#for generator cluster from exist MST:
./clust-mst -f -d 0.05 -t 48 -i bacteriaList.MinHashGenomeInfo bacteriaList.MinHashMSTInfo -o bacteria.clust
ATTENTION: the -f must in front of the -i option

#for generator cluster from exist sketches:
./clust-mst -E -d 0.05 -t 48 -i bacteriaList.MinHashGenomeInfo bacteriaList.MinHashSketchInfo -o bacteria.clust
./clust-greedy -f -d 0.05 -t 48 -i bacteriaList.MinHashGenomeInfo bacteriaList.MinHashSketchInfo -o bacteria.clust
#ATTENTION: the -E and -f must in front of the -i option

```

## Bug Report
All bug reports, comments and suggestions are welcome.

## Cite
RabbitTClust paper is under review now.
