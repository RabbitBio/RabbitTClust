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

#make rabbitIO library
cd RabbitIO
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=. ..
make && make install
export LD_LIBRARY_PATH=`pwd`/lib:$LD_LIBRARY_PATH
cd ../../

mkdir build && cd build
cmake -DUSE_RABBITIO=ON ..
make && make install
cd ../

#The refList is the list path of the RefSeq genome files.
./clust -l -t 48 -i refList -o ref.out

#get the clustering result by inputing MST info.
./clust -f -d 0.01 -i refListMinHashGenomeInfo refListMinHashMSTInfo -o result.out

#get more help info.
./clust -h

```

## Usage
```bash
Usage: clust [-h] [-l] [-t] <int> [-d] <double> [-F] <string> [-o] <string> [-i] <string>
Usage: clust [-h] [-f] [-d] <double> [-i] <string> <string>
-h          : this help message
-k <int>    : set kmer size, default <21>
-s <int>    : set sketch size, default <1000>
-l          : cluster for genomes(not sequences). 
-c          : compute the containment of genomes(sequences), cooperate with sketch function MinHash
-d <double> : set the threshold cluster from the Minimum Spanning Tree, default 0.05 (0.3 for containment)
-f          : generate cluster from the existing genomeInfo and MST content
-F <string> : set the sketch function, including <MinHash>, <WMH>, <OMH>, <HLL>, default <MinHash>
-o <string> : path of result file
-i <string> : path of input file. 
                (example 1) Default, path of input file, cluster for sequences within one file.
                (example 2) With the cooperation of '-l' option, list input. Lines in each <inputFile> specify paths to genome files, one per line.
                (example 3) With the cooperation of '-f' option, two input file, the former genomeInfo, the latter MST content.

Example as follows:

    #example 1: cluster for sequences in ref.fna
    ./clust -d 0.05 -t 48 -F MinHash -o ref.out -i ref.fna

    #example 2: cluster for genomes, 'refList' is list input.
    ./clust -l -d 0.05 -t 48 -F MinHash -o ref.out -i refList

    #example 3: cluster from the existing MST.
    #refListMinHashGenomeInfo and refMSTInfo are generated from former cluster as example 1 or example 2.
    ./clust -f -d 0.15 -i refListMinHashGenomeInfo refListMinHashMSTInfo -o ref.out

```

## Bug Report
All bug reports, comments and suggestions are welcome.