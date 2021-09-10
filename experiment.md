# data sets 

* viral dataset from RefSeq release 99
    * 9,033 genomes(11,562 sequences)
    * 318MB

* bacteria dataset from RefSeq release 99
    * 100,321 genomes(18,408,846 sequences)
    * 89,854 genomes(over 17 million sequences)

# platform 
* gold
* platinum
* amd2(default platform)

# software 
new thought:
* the fastANI says: Distribution of ANI values with each comparison labeled by the ***nomenclature of genomes*** being compared(in Fig.3c). What is about ***nomenclature of genomes***?
* the result of fastANI(based on MashMap to estimate the Average Nucleotide Identity) may be different with RabbitTClust in the method aspect but has some identity with the final result.(this is just a new thought, need to verify)
* the distance computing methods based on Mash is sensitive to the length of two sequences(genomes). Whether the fastANI has consider this question? (need to verify) 


## rabbitTClust
time on viral dataset on amd2 with parameters as follows:(other parameter is default)  
* run `./clust -l -t 48 -d 0.05 viralList >out`  
time is: **2.361s**


* run `./clust -l -t 48 -d 0.05 -F MinHash bacteriaList >out`  
time is: **44m57.519s** //sketchSize = 10000  
time is: **8m56.104s** //sketchSize = 1000


## fastANI
Time on viral dataset on amd2 with parameters as follows:(other parameter is default)

* run `./fastANI --rl viralList --ql viralList -k 21 -t 48 -o viral.out`  
time is: **1m6.081s**

* run `./fastANI --rl bacteriaList --ql bacteriaList -k 21 -t 48 -o bacteria.out`  
the memory footprint is over 1.5T which cannot run on the platinum workstation.
