# download the genomes from NCBI RefSeq and GenBank
RabbitTClust supports an input list of genomes in original FASTA format and gzips format.
We recommend using the decompressed genome files as input to filter out the broken download compressed files.
If you have to use the input list of the compressed files, you must check the md5 value.

For the input of a single FASTA file (each sequence means a genome), RabbitTClust only supports decompressed FASTA format.

## RefSeq version used in the study

The large-scale bacterial benchmark is fixed to [NCBI RefSeq Release 233](https://ncbiinsights.ncbi.nlm.nih.gov/2026/02/19/refseq-release-233/). Its versioned, release-wide accession catalog is archived as [`RefSeq-release233.catalog.gz`](https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/archive/RefSeq-release233.catalog.gz); bacterial records are identified by `bacteria` in the RefSeq release-directory field. Do not replace this snapshot with the live RefSeq bacteria `assembly_summary.txt`, which changes over time.

## Download genomes from RefSeq
The download script comes from [Bonsai](https://github.com/dnbaker/bonsai/tree/ac6f8c7ee1b2ae1128970a8f6dc01ddad19fdb37).

For general-purpose downloads from the current RefSeq directory, use `download_refseq.py` as follows:

* `python3 download_refseq.py bacteria`
* `python3 download_refseq.py -h` for additional options.

This convenience command accesses the live RefSeq directory and is not a replacement for the fixed Release 233 benchmark snapshot.

## Download genomes from GenBank

The latest release of GenBank bacterial genomes can be downloaded by `download_genbank.sh` as follows:
* `./download_genbank.sh`
