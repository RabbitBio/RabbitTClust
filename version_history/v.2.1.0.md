# `RabbitTClust v.2.1.0`
RabbitTClust is a fast and memory-efficient genome clustering tool based on sketch-based distance estimations.
It enables processing of large-scale datasets by combining dimensionality reduction techniques with streaming and parallelization on modern multi-core platforms.
RabbitTClust supports classical single-linkage hierarchical (clust-mst) and greedy incremental clustering (clust-greedy) algorithms for different scenarios. 

## Installation
RabbitTClust version 2.1.0 can only support 64-bit Linux Systems.

The detailed update information for this version, as well as the version history, can be found in the [`version_history`](./history.md) document.

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
# clust-mst, minimum-spanning-tree-based module for RabbitTClust
Usage: ./clust-mst [OPTIONS]
Options:
  -h,--help                   Print this help message and exit
  -t,--threads INT            set the thread number, default all CPUs of the platform
  -m,--min-length UINT        set the filter minimum length (minLen), genome length less than minLen will be ignore, default 10,000
  -c,--containment INT        use AAF distance with containment coefficient, set the containCompress, the sketch size is in proportion with 1/containCompress
  -k,--kmer-size INT          set the kmer size
  -s,--sketch-size INT        set the sketch size for Jaccard Index and Mash distance, default 1000
  -l,--inputlist              input is genome list, one genome per line
  -e,--no-save                not save the intermediate files, such as sketches or MST
  -d,--threshold FLOAT        set the distance threshold for clustering
  -F,--function TEXT          set the sketch function, such as MinHash, KSSD, default MinHash
  -o,--output TEXT REQUIRED   set the output name of cluster result
  -i,--input TEXT             set the input file
  --presketched TEXT          clustering by the pre-generated sketch files rather than genomes
  --premsted TEXT             clustering by the pre-generated mst files rather than genomes for clust-mst

# clust-greedy, greedy incremental clustering module for RabbitTClust
Usage: ./clust-greedy [OPTIONS]
Options:
  -h,--help                   Print this help message and exit
  -t,--threads INT            set the thread number, default all CPUs of the platform
  -m,--min-length UINT        set the filter minimum length (minLen), genome length less than minLen will be ignore, default 10,000
  -c,--containment INT        use AAF distance with containment coefficient, set the containCompress, the sketch size is in proportion with 1/containCompress
  -k,--kmer-size INT          set the kmer size
  -s,--sketch-size INT        set the sketch size for Jaccard Index and Mash distance, default 1000
  -l,--inputlist              input is genome list, one genome per line
  -e,--no-save                not save the intermediate files, such as sketches or MST
  -d,--threshold FLOAT        set the distance threshold for clustering
  -F,--function TEXT          set the sketch function, such as MinHash, KSSD, default MinHash
  -o,--output TEXT REQUIRED   set the output name of cluster result
  -i,--input TEXT             set the input file
  --presketched TEXT          clustering by the pre-generated sketch files rather than genomes
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
./clust-greedy -d 0.001 -l -i bacteria.list -o bacteria.out

#for last running of clust-mst, it generated a folder name in year_month_day_hour-minute-second format, such as 2023_05_06_08-49-15.
#this folder contains the sketch, mst files.
#for generator cluster from exist MST with a distance threshold of 0.045:
./clust-mst -d 0.045 --premsted 2023_05_06_08-49-15 -o bacteria.mst.d.045.clust
#for generator cluster from exist sketches files of clust-mst with a distance threshold of 0.045:
./clust-mst -d 0.045 --presketched 2023_05_06_08-49-15 -o bacteria.mst.d.045.clust

#for generator cluster from exist sketches of clust-greedy with a distance threshold of 0.001:
# folder 2023_05_06_08-49-15 contains the sketch files.
./clust-greedy -d 0.001 --presketched 2023_05_06_08-49-15 -o bact_genbank.greedy.d.001.clust
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
We highly appreciate all bug reports, comments, and suggestions from our users.  
Please feel free to raise any concerns or feedback with us without hesitation by `issue`. 

## Cite
[Xu, X. et al. (2022). RabbitTClust: enabling fast clustering analysis of
millions bacteria genomes with minhash sketches. bioRxiv.](https://doi.org/10.1101/2022.10.13.512052)
