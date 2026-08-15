[![install with conda](
https://anaconda.org/bioconda/rabbittclust/badges/version.svg)](https://anaconda.org/bioconda/rabbittclust)
[![install with conda](
https://anaconda.org/bioconda/rabbittclust/badges/latest_release_date.svg)](https://anaconda.org/bioconda/rabbittclust)
[![install with conda](
https://anaconda.org/bioconda/rabbittclust/badges/platforms.svg)](https://anaconda.org/bioconda/rabbittclust)
[![install with conda](
https://anaconda.org/bioconda/rabbittclust/badges/downloads.svg)](https://anaconda.org/bioconda/rabbittclust)

![RabbitTClust](rabbittclust.png)

# `RabbitTClust v2.4.0`
RabbitTClust is a fast and memory-efficient genome clustering tool based on sketch-based distance estimations.
It enables processing of large-scale datasets by combining dimensionality reduction techniques with streaming and parallelization on modern multi-core platforms.
RabbitTClust supports classical single-linkage hierarchical (clust-mst), greedy incremental clustering (clust-greedy), and graph-based clustering (clust-leiden) algorithms for different scenarios. 

## Benchmark data
The large-scale bacterial benchmark in this study uses the bacterial collection from [NCBI RefSeq Release 233](https://ncbiinsights.ncbi.nlm.nih.gov/2026/02/19/refseq-release-233/), not the continuously updated RefSeq bacteria directory. The version-pinned, release-wide accession catalog is available in the [NCBI RefSeq archive](https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/archive/RefSeq-release233.catalog.gz); bacterial records are identified by `bacteria` in the RefSeq release-directory field. Do not substitute the live `assembly_summary.txt`, because its contents change as RefSeq is updated.

## Installation
`RabbitTClust v2.4.0` supports 64-bit Linux systems.

The detailed update information for this version, as well as the version history, can be found in the [`version_history`](version_history/history.md) document.

### Install from bioconda 
RabbitTClust is available from [Bioconda](https://anaconda.org/bioconda/rabbittclust).

Ensure that your machine supports at least AVX2 instructions.


### Install from source code
#### Dependencies
* cmake v.3.0 or later
* C++17
* [zlib](https://zlib.net/)
* [igraph](https://igraph.org/) (optional, required for clust-leiden)

#### Compile and install
```bash
git clone --recursive https://github.com/RabbitBio/RabbitTClust.git
cd RabbitTClust
./install.sh
```

This will compile `clust-mst` and `clust-greedy` by default. If igraph is detected, `clust-leiden` will also be compiled.

#### Optional: Install igraph for clust-leiden
The `clust-leiden` module requires the [igraph](https://igraph.org/) library. If igraph is not found during installation, you will see a warning message, but `clust-mst` and `clust-greedy` will still be available.

**Option 1: Install via package manager** (if available)
```bash
# Ubuntu/Debian
sudo apt-get install libigraph-dev

# macOS
brew install igraph
```

**Option 2: Compile from source** (recommended for CentOS/RHEL)
```bash
cd ~
wget https://github.com/igraph/igraph/releases/download/0.10.10/igraph-0.10.10.tar.gz
tar xzf igraph-0.10.10.tar.gz
cd igraph-0.10.10
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/local
make -j8 && make install
```

After installing igraph, return to the RabbitTClust directory and run `./install.sh` again to compile `clust-leiden`.
## Usage

### Sketch selection

Use `--sketch` to select the sketch implementation. The default remains RabbitTClust's original `minhash` mode.

| `--sketch` value | Implementation and similarity estimator |
| --- | --- |
| `minhash` | Original RabbitTClust MinHash (default) |
| `kssd` | RabbitSketch KSSD; tune sampling with `--drlevel 0..3` |
| `fastkmv` | FastKMV bottom-k sketch |
| `probminhash` | Frequency-weighted ProbMinHash; estimates probability/weighted Jaccard |
| `bindash` | BinDash b-bit sketch |
| `hll` | HyperLogLog joint estimator |
| `setsketch` | SetSketch joint estimator |

All seven choices work with fresh FASTA or genome-list input in `clust-greedy`, `clust-mst`, and `clust-dbscan`. `--fast` remains available as a deprecated alias for `--sketch kssd`. For DBSCAN, KSSD retains its inverted-index accelerator; the other sketches use an exact all-pairs neighborhood search.

`-s/--sketch-size` is the requested sketch/register budget. BinDash rounds it up to a multiple of 64, while HLL and SetSketch round it up to a power of two with a minimum of 1024 registers. Distances are reported on the Mash-distance scale; ProbMinHash first estimates frequency-weighted Jaccard. Containment mode is limited to MinHash and KSSD.

The common CLI defaults are `-k 19`, `-s 1000`, `-d 0.05`, `-m 10000`, and all detected CPUs. The effective sketch defaults are:

| Sketch | Effective default configuration |
| --- | --- |
| MinHash | 1,000 hashes; the legacy automatic tuner may change MinHash parameters when neither `-k` nor `-s` is supplied |
| KSSD | CLI `k=19` is represented as an even `k=20`; `--drlevel 3` (sampling 1/4096); `-s` does not apply |
| FastKMV | 1,000 bottom-k values; seed 42 |
| ProbMinHash | 1,000 registers; frequency-weighted ProbMinHash4; seed 42; unlimited register updates (`max_L=0`); canonical k-mer occurrences sampled at 1/8 before frequency aggregation |
| BinDash | 1,024 bins (16 groups of 64); 16 bits per bin; seed 42 |
| HLL | 1,024 registers (`p=10`); Ertl joint-MLE Jaccard estimator |
| SetSketch | 1,024 registers (`p=10`); base 2, `a=5`; joint-MLE Jaccard estimator |

Persistent `--presketched`, `--append`, RepDB, and `--save-rep` workflows currently remain limited to MinHash and KSSD. The five additional implementations run in memory.

```bash
# MinHash is the default; these two commands are equivalent.
./clust-greedy -l -i bacteria.list -o minhash.cluster
./clust-greedy --sketch minhash -l -i bacteria.list -o minhash.cluster

# Select another RabbitSketch implementation.
./clust-mst --sketch probminhash -s 1024 -l -i bacteria.list -o probminhash.cluster
./clust-dbscan --sketch setsketch -s 1024 --eps 0.05 --minpts 5 -l -i bacteria.list -o setsketch.cluster
```

```bash
# clust-mst, minimum-spanning-tree-based module for RabbitTClust
Usage: ./clust-mst [OPTIONS]
Options:
  -h,--help                   Print this help message and exit
  -t,--threads INT            set the thread number, default all CPUs of the platform
  -m,--min-length UINT        set the filter minimum length (minLen), genome length less than minLen will be ignore, default 10,000
  -c,--containment INT        use AAF distance with containment coefficient, set the containCompress, the sketch size is in proportion with 1/containCompress  -k,--kmer-size INT          set the kmer size
  -s,--sketch-size INT        set the sketch size for Jaccard Index and Mash distance, default 1000
  -l,--list                   input is genome list, one genome per line
  -e,--no-save                not save the intermediate files, such as sketches or MST
  -d,--threshold FLOAT        set the distance threshold for clustering
  -o,--output TEXT REQUIRED   set the output name of cluster result
  -i,--input TEXT Excludes: --append
                              set the input file, single FASTA genome file (without -l option) or genome list file (with -l option)
  --presketched TEXT          clustering by the pre-generated sketch files rather than genomes
  --sketch TEXT               sketch algorithm: minhash, kssd, fastkmv, probminhash, bindash, hll, or setsketch (default: minhash)
  --premsted TEXT             clustering by the pre-generated mst files rather than genomes for clust-mst
  --newick-tree               output the newick tree format file for clust-mst
  --fast                      deprecated alias for --sketch kssd
  --dense                     optional: enable density/ANI stats and MST noise-removal pass (high memory; default is off)
  --dedup-dist FLOAT          within each cluster, collapse near-duplicate nodes connected by forest edges with dist <= dedup-dist; output to <output>.dedup
  --reps-per-cluster INT      select up to k representatives per cluster (after optional dedup); output to <output>.reps
  --append TEXT Excludes: --input
                              append genome file or file list with the pre-generated sketch or MST files

# clust-greedy, greedy incremental clustering module for RabbitTClust
Usage: ./clust-greedy [OPTIONS]
Options:
  -h,--help                   Print this help message and exit
  -t,--threads INT            set the thread number, default all CPUs of the platform
  -m,--min-length UINT        set the filter minimum length (minLen), genome length less than minLen will be ignore, default 10,000
  -c,--containment INT        use AAF distance with containment coefficient, set the containCompress, the sketch size is in proportion with 1/containCompress  -k,--kmer-size INT          set the kmer size
  -s,--sketch-size INT        set the sketch size for Jaccard Index and Mash distance, default 1000
  -l,--list                   input is genome list, one genome per line
  -e,--no-save                not save the intermediate files, such as sketches or MST
  -d,--threshold FLOAT        set the distance threshold for clustering
  -o,--output TEXT REQUIRED   set the output name of cluster result
  -i,--input TEXT Excludes: --append
                              set the input file, single FASTA genome file (without -l option) or genome list file (with -l option)
  --presketched TEXT          clustering by the pre-generated sketch files rather than genomes
  --sketch TEXT               sketch algorithm: minhash, kssd, fastkmv, probminhash, bindash, hll, or setsketch (default: minhash)
  --append TEXT Excludes: --input
                              append genome file or file list with the pre-generated sketch or MST files
  --save-rep                  save representative inverted index for incremental clustering (note: may slightly affect performance)
  --dense                     optional: enable density/ANI-related MST post-processing (high memory; default is off; KSSD presketched path)

# clust-leiden, graph-based clustering module for RabbitTClust (requires igraph)
Usage: ./clust-leiden [OPTIONS]
Options:
  -h,--help                   Print this help message and exit
  -t,--threads INT            set the thread number, default all CPUs of the platform
  -m,--min-length UINT        set the filter minimum length (minLen), genome length less than minLen will be ignore, default 10,000
  -k,--kmer-size INT          set the kmer size
  -l,--list                   input is genome list, one genome per line
  -e,--no-save                not save the intermediate files, such as sketches
  -d,--threshold FLOAT        set the distance threshold for graph edge construction
  -o,--output TEXT REQUIRED   set the output name of cluster result
  -i,--input TEXT Excludes: --presketched
                              set the input file, single FASTA genome file (without -l option) or genome list file (with -l option)
  --presketched TEXT          clustering by the pre-generated sketch files rather than genomes
  --pregraph TEXT             clustering from pre-built graph (fast resolution adjustment without rebuilding graph)
  --fast                      use the kssd algorithm for sketching and distance computing (required)
  --resolution FLOAT          resolution parameter for clustering (higher = more clusters, default 1.0)
  --louvain                   use Louvain algorithm instead of Leiden (default: Leiden)
  --drlevel INT               set the dimension reduction level for Kssd sketches, default 3 with a dimension reduction of 1/4096
```

## Example:
```bash
# input is a file list, one genome path per line:
./clust-mst -l -i bact_refseq.list -o bact_refseq.mst.clust
./clust-greedy -l -i bact_genbank.list -o bact_genbank.greedy.clust

# input is a single genome file in FASTA format, one genome as a sequence:
./clust-mst -i bacteria.fna -o bacteria.mst.clust
./clust-greedy -i bacteria.fna -o bacteria.greedy.clust

# the sketch size (reciprocal of sampling proportion), kmer size, and distance threshold can be specified by -s (-c), -k, and -d options.
./clust-mst -l -k 21 -s 1000 -d 0.05 -i bact_refseq.list -o bact_refseq.mst.clust
./clust-greedy -l -k 21 -c 1000 -d 0.05 -i bact_genbank.list -o bact_genbank.greedy.clust


# for redundancy detection with clust-greedy, input is a genome file list:
# use -d to specify the distance threshold corresponding to various degrees of redundancy.
./clust-greedy -d 0.001 -l -i bacteria.list -o bacteria.out

# v.2.1.0 or later
# for last running of clust-mst, it generated a folder name in year_month_day_hour-minute-second format, such as 2023_05_06_08-49-15.
# this folder contains the sketch, mst files.
# for generator cluster from exist MST with a distance threshold of 0.045:
./clust-mst -d 0.045 --premsted 2023_05_06_08-49-15/ -o bact_refseq.mst.d.045.clust
# for generator cluster from exist sketches files of clust-mst with a distance threshold of 0.045:
./clust-mst -d 0.045 --presketched 2023_05_06_08-49-15/ -o bact_refseq.mst.d.045.clust

# for generator cluster from exist sketches of clust-greedy with a distance threshold of 0.001:
# folder 2023_05_06_08-49-15 contains the sketch files.
./clust-greedy -d 0.001 --presketched 2023_05_06_09-37-23/ -o bact_genbank.greedy.d.001.clust

# v.2.2.0 or later
# for generator cluster from exist part sketches (presketch_A_dir) and append genome set (genome_B.list) to incrementally clustering 
./clust-mst --presketched 2023_05_06_08-49-15/ -l --append genome_B.list -o append_refseq.mst.clust
./clust-mst --presketched 2023_05_06_09-37-23/ -l --append genome_B.list -o append_genbank.greedy.clust

# v.2.2.1 or later
# output the newick tree format for clust-mst, use the --newick-tree flag.
./clust-mst -l -i bacteria.list --newick-tree -o bacteria.mst.clust 

# v.2.3.0 or later
# use the efficient Kssd sketch strategy for clust-mst, use the --fast flag.
./clust-mst --fast -l -i bacteria.list -o bacteria.fast.mst.clust

# v.2.4.0 or later
# clust-greedy also supports the efficient Kssd sketch strategy with the --fast flag.
./clust-greedy --fast -l -i bacteria.list -o bacteria.fast.greedy.clust

# clust-greedy supports incremental clustering with --presketched and --append
./clust-greedy --fast --presketched 2023_05_06_09-37-23/ -l --append genome_B.list -o append.greedy.clust

# clust-leiden: graph-based clustering with Leiden algorithm (default)
# does not require a strict distance threshold, automatically finds community structure
./clust-leiden --fast -l -i bacteria.list -o bacteria.leiden.clust

# clust-leiden with Louvain algorithm (faster, simpler alternative)
./clust-leiden --fast -l -i bacteria.list -o bacteria.louvain.clust --louvain

# clust-leiden with pre-sketched data and pre-graph
./clust-leiden --fast --presketched 2023_05_06_09-37-23/ -o bacteria.leiden.clust
./clust-leiden --pregraph 2023_05_06_09-37-23 -o bacteria_r20.leiden.clust --resolution 2.0
```


## Output
The output file is in a CD-HIT output format and is slightly different when running with or without `-l` input option.  
When using the `-l` option, the input is expected to be a FASTA file list, with each file representing a genome. Without the `-l` option, the input should be a single FASTA file, with each sequence representing a genome.

#### Output format for a FASTA file list input
With `-l*` option, the tab-delimited values in the lines beginning with tab delimiters are:
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
Without `-l` option, the tab-delimited values in the lines beginning with tab delimiters are:
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

#### Output the newick tree format (v.2.2.1 or latter)
When the `--newick-tree` option is used, an additional output file will be generated in the Newick tree format with a suffix name of ".newick.tree".


# Bug Report
We highly appreciate all bug reports, comments, and suggestions from our users.  
Please feel free to raise any concerns or feedback with us without hesitation by `issue`. 

## Cite

- Zhang T, Xu X, Yin Z, et al. RabbitTClust2: Fast, Scalable, and Versatile Clustering for Massive Genomic Datasets. In: *2025 IEEE International Conference on Bioinformatics and Biomedicine (BIBM)*. IEEE; 2025:1495-1502. https://doi.org/10.1109/BIBM66473.2025.11356063

- Xu X, Yin Z, Yan L, et al. RabbitTClust: enabling fast clustering analysis of millions of bacteria genomes with MinHash sketches. *Genome Biology*. 2023;24:121. https://doi.org/10.1186/s13059-023-02961-6
