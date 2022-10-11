# evaluation script  
Run `make` to compile the script.

## calLabel && getNMI.py
The `calLabel` is used for generating the label file from the clustering result of RabbitTClust. 

**==============================================================================**  
Example:  
`./calLabel bacteria.groundTruth -l bacteria.mst.clust bacteria.mst.label`  
It will generate the label file (`bacteria.mst.label`) from the clustering result (`bacteria.mst.clust`) of RabbitTClust.  
Subsequently, run `python3 getNMI.py bacteria.mst.label` to get the NMI score.

**==============================================================================**  
* run as: `./calLabel groundTruth sketchOption clustFile labelFile`
  * The 0 parameter(`./calLabel`) is the application name
  * The 1 parameter(`groundTruth`) is input file, groundTruth file, `<accession, taxid, organismName>` per line (first line is header)
  * The 2 parameter(`sketchOption`) is input option, sketch options, `-l` means sketchByFile (input as a genome list), `-i` means sketchBySequence (input as a single genome file)
  * The 3 parameter(`clustFile`) is input file, cluster result file generated from RabbitTClust
  * The 4 parameter(`labelFile`) is output file, label file according to the groundTruth 

## calPurity 
The `calPurity` is used for computing the purity of the clustering result of RabbitTClust.

**==============================================================================**  
Example:  
`./calPurity -l bacteria.groundTruth bacteria.clust bacteria.purity`  
It will compute the total purity and the coverage of the clustering result and generate three information files: `bacteria.purity`, `bacteria.purity.accession.purity`, and `bacteria.purity.accession.unpurity`.
* `bacteria.purity` is the detail purity for each cluster.
* `bacteria.purity.accession.purity` is the list of first genome in each purity cluster.
* `bacteria.purity.accession.unpurity` is the list of first dominant genome and the unpurity genomes for each cluster.

**==============================================================================**  
* run as: `./calPurity options(-l, -i) groundTruth clustFile bacteria.purity`  
  * The 0 parameter(`./calPurity`) is the application name
  * The 1 parameter(`options(-l, -i)`) is input option, sketch option for clust, -l or -i
  * The 2 parameter(`groundTruth`) is input file, the groundTruth file, `<assembly_accession, species_taxid, genomeName>` per line
  * The 3 parameter(`clustFile`) is input file, cluster result file from RabbitTClust
  * The 4 parameter(`bacteria.purity`) is output file, output purity info file, including total result file(`bacteria.purity`) and accession file(`<accession, taxid>` per line)

## getRepresentativeList
The `getRepresentativeList` is used for generating the representative genome list from the clustering result.
**==============================================================================**  
Example:  
`./getRepresentativeList -l bacteria.greedy.clust bacteria_representative.list`  
It will choose the representative genome for each cluster and generate a list of these representative genomes.
**==============================================================================**  
run as: `./getRepresentativeList -i/-l clustFile representative_list`  
  * The 0 parameter(`./getRepresentativeList`) is the application name
  * The 1 parameter(`-i/-l`) is input parameter, sketch parameter for the cluster file, -i means sketchBySequence, -l means sketchByFile
  * The 2 parameter(`clustFile`) is input file, the cluster result from RabbitTClust
  * The 3 parameter(`representative_list`) is output file, the representative list of genomes file or sequences name.








