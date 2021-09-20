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
* node1 3090 

# software 
new thought:
* the fastANI says: Distribution of ANI values with each comparison labeled by the ***nomenclature of genomes*** being compared(in Fig.3c). What is about ***nomenclature of genomes***?
* the result of fastANI(based on MashMap to estimate the Average Nucleotide Identity) may be different with RabbitTClust in the method aspect but has some identity with the final result.(this is just a new thought, need to verify)
* the distance computing methods based on Mash is sensitive to the length of two sequences(genomes). Whether the fastANI has consider this question? (need to verify) 


## rabbitTClust
Run time and parameters of bacteria database on different platform, clust genomes(with parameter -l)  
note: 
* **k:** kmer size
* **Containment:** use containment of minHash
* **march:** SIMD instruction level
* **T_Sketch:** time of sketch generating
* **T_Dist_MST:** time of distance matrix computing and MST generating

|database   |platform|k  |sketchSize   |threads|Containment|    march |time   |T_Sketch|T_Dist_MST|
|:-:        |:-:     |:-:|:-:          |:-:    |:-:        |    :-:   |:-:    |:-:     |:-:       |
|bacteria   |amd2    |21 |length/1000  |48     |true       |  -O2  none  |160m48s|551s    |9070s     |
|bacteria   |amd2    |21 |length/10000 |48     |true       |  -O2  none  | 20m43s|279s    |932s      |
|bacteria   |amd2    |21 |1000         |48     |false      |  -O2  avx2  |  9m50s|322s    |242s      |
|bacteria   |amd2    |21 |1000         |48     |false      |  -O3  avx2  |  9m46s|311s    |251s      |
|bacteria   |amd2    |21 |1000         |48     |false      |  -O2  none  | 23m25s|296s    |1084s     |
|bacteria   |amd2    |21 |1000         |48     |false      |  -O3  none  |  24m3s|287s    |1130s     |
|bacteria   |platinum|21 |1000         |48     |false      |  -O2  avx512| 11m40s|350s    |326s      |
|bacteria   |platinum|21 |1000         |48     |false      |  -O3  avx512|  9m55s|353s    |216s      |
|bacteria   |platinum|21 |1000         |48     |false      |  -O3  avx2  | 10m32s|367s    |240s      |
|bacteria   |platinum|21 |1000         |48     |false      |  -O2  avx2  | 11m26s|362s    |298s      |
|bacteria   |platinum|21 |1000         |48     |false      |  -O3  none  | 33m33s|345s    |1643s     |
|bacteria   |platinum|21 |1000         |48     |false      |  -O2  none  | 29m57s|341s    |1431s     |
|bacteria   |gold    |21 |1000         |40     |false      |  -O2  avx512| 12m11s|462s    |242s      |
|bacteria   |gold    |21 |1000         |40     |false      |  -O3  avx512| 14m26s|449s    |388s      |
|bacteria   |node1   |21 |1000         |10     |false      |  -O2  avx2  |220m35s|12534s  |680s      |


interprate for the result:
* for resemblance:
  * time of sketch on node1 is much more than that on amd2 since the 335G bacteria database is saved in SSD on amd2 while in hard disk on node1. 
  * Time of distance computing without vectorization(march=none) are much larger since the vectorization implemented by SIMD instructions can avoid branch misprediction efficiently when computing intersection.
* for containment:
  * In the table, the sketch Size of row 1 is 10 times of row 2, the distance computing is on proportion of sketchSize , so the T_Dist_MST is  10 times of line 2 as well.
  * compared with row 1, row 2 has less sketch time since the overhead of maintaining the MinHashHeap with smaller size is less.(the sketchSize of row 1 is distributed on the range of 1000 to 10000 while the sketch size range of row 2 is between 100 to 1000 with the same distribution.)

questions for the experient result for the table:  
* [ ] For RabbitSketch, the -O2 is faster than -O3 compile option, especially in distance computing.

## fastANI
Time on viral dataset on amd2 with parameters as follows:(other parameter is default)

* run `./fastANI --rl viralList --ql viralList -k 21 -t 48 -o viral.out`  
time is: **1m6.081s**

* run `./fastANI --rl bacteriaList --ql bacteriaList -k 21 -t 48 -o bacteria.out`  
the memory footprint is over 1.5T which cannot run on the platinum workstation.
