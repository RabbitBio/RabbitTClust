# several source files are in different usage 
the compile.sh is used to compile the src file.

* mapGenome.cpp: 
  * it is used to check whether there are different namenclature type in the same genome.  
  * run with `./mapGenome bacteriaList`, the result is saved in mapType.out.
* calGenome.cpp: 
  * it is used to calculate the number of each namenclature type of first two name in comment.  
  * run with `./calGenome2 bacteriaList`, the result is saved in genomeType2.out
* calResult.cpp: 
  * it is used to calculate the evaluation cluster result. It will output the purity, precision, recall, rand index and F-score of the result.  
  * run with `./calResult2 bacteria.out`, the bacteria.out is the output result file of RabbitTClust.






