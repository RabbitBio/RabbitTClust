# several source files are in different usage 
The compile.sh is used to compile the src file.

* mapGenome.cpp: 
  * It is used to check whether there are different namenclature type in the same genome.  
  * Run with `./mapGenome bacteriaList`, the result is saved in mapType.out.
* calGenome.cpp: 
  * It is used to calculate the number of each namenclature type of first two name in comment.  
  * Run with `./calGenome2 bacteriaList`, the result is saved in genomeType2.out
* calResult.cpp: 
  * It is used to calculate the evaluation cluster result. It will output the purity, precision, recall, rand index and F1-score of the result.  
  * Run with `./calResult1 bacteria.out` or `./calResult2 bacteria.out` for differnt label levels. The bacteria.out is the output result file of RabbitTClust.






