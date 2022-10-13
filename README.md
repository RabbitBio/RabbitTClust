![RabbitTClust](rabbittclust.png)

RabbitTClust is a fast and memory-efficient genome clustering tool based on sketch-based distance estimations.
It enables processing of large-scale datasets by combining dimensionality reduction techniques with streaming and parallelization on modern multi-core platforms.
RabbitTClust supports classical single-linkage hierarchical (clust-mst) and greedy incremental clustering (clust-greedy) algorithms for different scenarios. 

## Installation
RabbitTClust version 2.0 can only support 64-bit Linux Systems.

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
cd RabbitSketch &&
mkdir -p build && cd build &&
cmake -DCXXAPI=ON -DCMAKE_INSTALL_PREFIX=. .. &&
make -j8 && make install &&
cd ../../ &&

#make rabbitFX library
cd RabbitFX && 
mkdir -p build && cd build &&
cmake -DCMAKE_INSTALL_PREFIX=. .. &&
make -j8 && make install && 
cd ../../ &&

#compile the clust-greedy
mkdir -p build && cd build &&
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=ON .. && 
make -j8 && make install &&
cd ../ &&

#compile the clust-mst
cd build &&
cmake -DUSE_RABBITFX=ON -DUSE_GREEDY=OFF .. &&
make -j8 && make install &&
cd ../ 
```

## Usage
```bash
usage: clust-mst [-h] [-l] [-t] <int> [-d] <double> [-F] <string> [-i] <string> [-o] <string>
usage: clust-mst [-h] [-f] [-E] [-d] <double> [-i] <string> <string> [-o] <string>
usage: clust-greedy [-h] [-l] [-t] <int> [-d] <double> [-F] <string> [-i] <string> [-o] <string>
usage: clust-greedy [-h] [-f] [-d] <double> [-i] <string> <string> [-o] <string>
-h         : this help message
-k <int>   : set kmer size, default 21, for both clust-mst and clust-greedy
-s <int>   : set sketch size, default 1000, for both clust-mst and clust-greedy
-c <int>   : set sampling ratio to compute variable sketchSize, sketchSize = genomeSize/samplingRatio, only support with MinHash sketch function, for clust-greedy
-d <double>: set the distance threshold, default 0.05 for both clust-mst and clust-greedy
-t <int>   : set the thread number, default take full usage of platform cores number, for both clust-mst and clust-greedy
-l         : input is a file list, not a single genome file. Lines in the input file list specify paths to genome files, one per line, for both clust-mst and clust-greedy
-i <string>: path of input file. One file list or single genome file. Two input file with -f and -E option
-f         : two input files, genomeInfo and MSTInfo files for clust-mst; genomeInfo and sketchInfo files for clust-greedy 
-E         : two input files, genomeInfo and sketchInfo for clust-mst
-F <string>: set the sketch function, including MinHash and KSSD, default MinHash, for both clust-mst and clust-greedy
-o <string>: path of output file, for both clust-mst and clust-greedy
-e         : not save the intermediate file generated from the origin genome file, such as the GenomeInfo, MSTInfo, and SketchInfo files, for both clust-mst and clust-greedy

```

## Example:
```bash
#input is a file list, one genome path per line:
./clust-mst -l -i bact_refseq.list -o bact_refseq.mst.clust
./clust-greedy -l -i bact_genbank.list -o bact_genbank.greedy.clust

#input is a single genome file in FASTA format, one genome as a sequence:
./clust-mst -i bacteria.fna -o bacteria.mst.clust
./clust-greedy -i bacteria.fna -o bacteria.greedy.clust

#the sketch size (reciprocal of sampling proportion), kmer size, and distance threshold can be specified by -s (-c), -k, and -d options.
./clust-mst -l -k 21 -s 1000 -d 0.05 -i bact_refseq.list -o bact_refseq.mst.clust
./clust-greedy -l -k 21 -c 1000 -d 0.05 -i bact_genbank.list -o bact_genbank.greedy.clust


#for redundancy detection with clust-greedy, input is a genome file list:
#use -d to specify the distance threshold corresponding to various degrees of redundancy.
./clust-greedy -d 0.001 -l -i bacteriaList -o bacteria.out

#for generator cluster from exist MST with a distance threshold of 0.045:
#ATTENTION: the -f must in front of the -i option
./clust-mst -d 0.05 -f -i bact_refseq.list.MinHashGenomeInfo bact_refseq.list.MinHashMSTInfo -o bact_refseq.mst.d.045.clust

#for generator cluster from exist sketches of clust-greedy with a distance threshold of 0.001:
#ATTENTION: the -f must in front of the -i option
./clust-greedy -d 0.001 -f -i bact_genbank.list.MinHashGenomeInfo bact_genbank.list.MinHashSketchInfo -o bact_genbank.greedy.d.001.clust

```

## Bug Report
All bug reports, comments and suggestions are welcome.

## Cite
RabbitTClust paper is under review now.
