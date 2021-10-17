# benchmark of cluster evaluation and simulate databese

* sub-directory evaluation is used for cluster evaluation.
* sub-directory simulateData is used for simulating genome databases with very different size for containment evaluation.

# download genomes from RefSeq
The download script is from [Bonsai](https://github.com/dnbaker/bonsai/tree/ac6f8c7ee1b2ae1128970a8f6dc01ddad19fdb37).

RefSeq genomes can be downloaded by download_genomes.py as follows:

`python3 download_genomes.py all`

Parameter `all` includes ***archaea***, ***bacteria***, ***fungi***, ***viral***, ***plant***, ***protozoa***, ***human***, ***vertebrate_mammalian*** and ***vertebrate_other***.

 Use `python3 download_genomes.py -h` for more details.
	
