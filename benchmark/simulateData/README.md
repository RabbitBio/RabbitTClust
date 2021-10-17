# several source files are in different usage
the compile.sh is used to compile the src file.
* redundancyRemove.cpp: 
  * It is used as a filter to choose the refGenomes for generator subGenomes with very different genome size.
  * Run with `./redunRemoveExec bacteriaList`, the refGenome fileList is saved in nonRedundantList.
  * You should copy the genome files into the noRedundant directory according to the nonRedundantList.
* simulate.cpp: 
  * It is used to generate the subGenomes of refGenomes with very different length.
  * Run with `./simulateExec noRedundant`,(the noRedundant is the fileList including the file path of refGenomes). 
  * The subGenomes is saved in the same directory of refGenomes.
  * You can run RabbitTClust with containment option with the database including refGenomes and subGenomes.
* createBacteria.cpp:
	* It is used for comparing the performance of RabbitTClust and gclust since gclust cannot finish the cluster of the whole bacteria genomes.
	* Run with `./createBacteria bacteriaList`, (the bacteriaList is the fileList including the file path for each line). 
	* You can run RabbitTClust and gclust on subBacteria.fna database to compare the performance.
