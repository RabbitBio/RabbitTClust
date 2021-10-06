# several source files are in different usage
the compile.sh is used to compile the src file.
* redundancyRemove.cpp: 
  * it is used as a filter to choose the refGenomes for generator subGenomes with very different genome size.
  * run with `./redunRemoveExec bacteriaList`, the refGenome fileList is saved in nonRedundantList.
  * you should copy the genome files into the noRedundant directory according to the nonRedundantList.
* simulate.cpp: 
  * it is used to generate the subGenomes of refGenomes with very different length.
  * run with `./simulateExec noRedundant`,(the noRedundant is the fileList including the file path of refGenomes). 
  * The subGenomes is saved in the same directory of refGenomes.
  * you can run rabbitTClust with containment option with the database including refGenomes and subGenomes.
