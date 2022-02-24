![RabbitTClust](rabbittclust.png)

RabbitTClust is an efficient sketch-based genome clustering application 
using minimum spanning trees. It enables processing of large-scale 
datasets by combining dimensionality reduction techniques with streaming 
and parallelization on modern multi-core platforms. 

## Installation
RabbitTClust version 1.0

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
make && make install
cd ../../

#make rabbitFX library
cd RabbitFX
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make && make install
cd ../../

mkdir build && cd build
cmake -DUSE_RABBITFX=ON ..
make && make install
cd ../

```

## Usage
```bash
Usage: clust [-h] [-l] [-t] <int> [-d] <double> [-F] <string> [-i] <string> [-o] <string> 
Usage: clust [-h] [-f] [-d] <double> [-i] <string> <string> [-o] <string>
-h          : this help message
-k <int>    : set kmer size, default <21>
-s <int>    : set sketch size, default <1000>
-l          : list input. Lines in each <input> specify paths to genome files, one per line.
-c <int>    : compute the containment of genomes, set proportion sketchSize = genomeSize/compress, ATTENTION with MinHash function. 
-d <double> : set the threshold cluster from the Minimum Spanning Tree, default 0.05 (0.3 for containment)
-f          : generate cluster from the existing genomeInfo and MST content,
-F <string> : set the sketch function, including <MinHash>, <WMH>, <OMH>, <HLL>, default <MinHash>
-o <string> : path of result file
-i <strings>: path of input files. 

```

## Example:
```bash

#The bacteria.fna is a single files for multi-genomes .
./clust -t 48 -i data/bacteria.fna -o bacteria.clust

#The refList is the list path of the RefSeq genome files.
./clust -l -t 48 -i data/refList -o ref.clust

#For redundancy detection, run with containment:
#input is a single file:
./clust -c 1000 -t 48 -i data/data.fna -o data.out

#For redundancy detection, run with containment:
#input is a file list:
./clust -l -c 1000 -t 48 -i data/fileList -o file.out

#get the clustering result by inputing MST info.
#ATTENTION the -f option must in front of the -i option
./clust -f -d 0.05 -i refListMinHashGenomeInfo refListMinHashMSTInfo -o result.clust

#get more help info.
./clust -h

```

## Bug Report
All bug reports, comments and suggestions are welcome.

## Cite
RabbitTClust paper is under review now.
