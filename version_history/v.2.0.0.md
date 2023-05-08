# `v.2.0.0`
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
## Output
The output file is in a CD-HIT output format and is slightly different when running with varying input options (*-l* and *-i*).  
Option *-l* means input as a FASTA file list, one file per genome, and *-i* means input as a single FASTA file, one sequence per genome.

#### Output format for a FASTA file list input
With *-l* option, the tab-delimited values in the lines beginning with tab delimiters are:
* local index in a cluster
* global index of the genome
* genome length
* genome file name (including genome assembly accession number)
* sequence name (first sequence in the genome file)
* sequence comment (remaining part of the line)

**Example:**
```txt
the cluster 0 is:
    0   0   14782125nt  bacteria/GCF_000418325.1_ASM41832v1_genomic.fna     NC_021658.1     Sorangium cellulosum So0157-2, complete sequence
    1   1   14598830nt  bacteria/GCF_004135755.1_ASM413575v1_genomic.fna    NZ_CP012672.1   Sorangium cellulosum strain So ce836 chromosome, complete genome

the cluster 1 is:
    0   2   14557589nt  bacteria/GCF_002950945.1_ASM295094v1_genomic.fna    NZ_CP012673.1   Sorangium cellulosum strain So ce26 chromosome, complete genome

the cluster 2 is:
    0   3   13673866nt  bacteria/GCF_019396345.1_ASM1939634v1_genomic.fna   NZ_JAHKRM010000001.1    Nonomuraea guangzhouensis strain CGMCC 4.7101 NODE_1, whole genome shotgun sequence

......
```

#### Output format for a single FASTA file input
With *-i* option, the tab-delimited values in the lines beginning with tab delimiters are:
* local index in a cluster
* global index of the genome
* genome length
* sequence name 
* sequence comment (remaining part of this line)

**Example:**
```txt
the cluster 0 is:
    0   0   11030030nt  NZ_GG657755.1   Streptomyces  himastatinicus ATCC 53653 supercont1.2, whole genome shotgun sequence
    1   1   11008137nt  NZ_RIBZ01000339.1   Streptomyces  sp. NEAU-LD23 C2041, whole genome shotgun sequence

the cluster 1 is:
    0   2   11006208nt  NZ_KL647031.1   Nonomuraea  candida strain NRRL B-24552 Doro1_scaffold1, whole genome shotgun sequence
    
the cluster 2 is:
    0   3   10940472nt  NZ_VTHK01000001.1   Amycolatopsis anabasis strain EGI 650086 RDPYD18112716_A.Scaf1, whole genome shotgun sequence

......
```


# Bug Report
All bug reports, comments and suggestions are welcome.

## Cite
[Xu, X. et al. (2022). RabbitTClust: enabling fast clustering analysis of
millions bacteria genomes with minhash sketches. bioRxiv.](https://doi.org/10.1101/2022.10.13.512052)
