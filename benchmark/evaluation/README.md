# several source files are in different usage 
The compile.sh is used to compile the src file.

* mapGenome.cpp: 
  * It is used to check whether there are different namenclature type in the same genome.  
  * Run with `./mapGenome bacteriaList`, the result is saved in mapType.out.

* calGenome.cpp: 
  * It is used to calculate the number of each namenclature type of first two name in comment.  
  * Run with `./calGenome2 bacteriaList`, the result is saved in genomeType2.out

* calLabel.cpp: 
  * It is the first step to calculate the precision, recall, F1-score, and NMI of cluster result.
* getResult.py:
	* It is the second step to caculate the precision, recall, F1-score, and NMI of the cluster result.








