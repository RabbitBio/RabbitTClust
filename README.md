# RabbitTClust

RabbitTClust is a clustering method based on Minimum Spanning Tree. The distance computing is based on Sketch.

## Installation

### Dependancy
* c++14
* [zlib](https://zlib.net/)

### Starting
```bash
git clone --recursive git@github.com:RabbitBio/RabbitTClust.git
#make rabbitSketch library
cd RabbitSketch
mkdir build && cd build 
cmake -DCXXAPI=ON -DCMAKE_INSTALL_PREFIX=. ..
make && make install
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
cd ../../

#make rabbitFX library
cd RabbitFX
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make && make install
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
cd ../../

mkdir build && cd build
cmake -DUSE_RABBITFX=ON ..
make && make install
cd ../

```

## Usage
```bash
Usage: clust [-h] [-l] [-t] <int> [-d] <double> [-F] <string> [-o] <string> [-i] <string>
Usage: clust [-h] [-f] [-d] <double> [-i] <string> <string>
-h          : this help message
-k <int>    : set kmer size, default <21>
-s <int>    : set sketch size, default <1000>
-l          : list input. Lines in each <input> specify paths to genome files, one per line.
-c <int>    : compute the containment of genomes, set proportion sketchSize = genomeSize/compress, ATTENTION with MinHash function. 
-d <double> : set the threshold cluster from the Minimum Spanning Tree, default 0.05 (0.3 for containment)
-f          : generate cluster from the existing genomeInfo and MST content
-F <string> : set the sketch function, including <MinHash>, <WMH>, <OMH>, <HLL>, default <MinHash>
-o <string> : path of result file
-i <string> : path of input file. 
                (example 1) Default, path of input file, cluster for sequences within one file.
                (example 2) With the cooperation of '-l' option, list input. Lines in each <inputFile> specify paths to genome files, one per line.
                (example 3) With the cooperation of '-f' option, two input file, the former genomeInfo, the latter MST content.

```

## Example:
```bash

#The refList is the list path of the RefSeq genome files.
./clust -l -t 48 -i data/refList -o ref.clust

#The bacteria.fna is a single files for multi-genomes .
./clust -t 48 -i data/bacteria.fna -o bacteria.clust

#For redundancy detection, run with containment:
#input is a file list:
./clust -l -c 1000 -t 48 -i data/fileList -o file.out

#For redundancy detection, run with containment:
#input is a single file:
./clust -c 1000 -t 48 -i data/data.fna -o data.out

#get the clustering result by inputing MST info.
#ATTENTION the -f option must in front of the -i option
./clust -f -d 0.05 -i refListMinHashGenomeInfo refListMinHashMSTInfo -o result.clust

#get more help info.
./clust -h

```

## Bug Report
All bug reports, comments and suggestions are welcome.
