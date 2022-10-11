# download the genomes from NCBI RefSeq and GenBank
RabbitTClust supports an input list of genomes in original FASTA format and gzips format.
We recommend using the decompressed genome files as input to filter out the broken download compressed files.
If you have to use the input list of the compressed files, you must check the md5 value.

For the input of a single FASTA file (each sequence means a genome), RabbitTClust only supports decompressed FASTA format.

## download genomes from RefSeq
The download script comes from [Bonsai](https://github.com/dnbaker/bonsai/tree/ac6f8c7ee1b2ae1128970a8f6dc01ddad19fdb37).

RefSeq bacterial genomes can be downloaded by `download_refseq.py` as follows:

* `python3 download_genomes.py bacteria`  
* `python3 download_genomes.py -h` more details of help infos.

## download genomes from GenBank
The FTP paths of the GenBank assembled bacterial genomes are listed in `bact_GenBank.list.gz`, which is generated from [assembly_summary_genbank.txt](https://ftp.ncbi.nlm.nih.gov/genomes/genbank/).

GenBank bacterial genomes can be downloaded by `download_genbank.sh` as follows:
* `./download_genbank.sh`
