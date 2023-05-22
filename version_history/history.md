# Latest version: `v.2.2.1` 
* add `--newick-tree` option to output the Newick tree format for `clust-mst`.

## [`v.2.2.0`](v.2.2.0.md)
* support incrementally clustering by option `--append` accompanied with `--presketched` or `--premsted` options.

Note:  
* When considering the clustering of the genome set `A+B` using a pre-generated sketch `A_sketch` and an appending genome set `B`, it is important to note that the sketch parameter for the pre-generated sketch `A_sketch` and the appending set `B` may differ from that of the whole genome set `A+B`. However, the impact of changes in the genome lengths of set `B` on the automatically generated parameters will be minimal if they are not significant.

    * This is because the sketch parameters, including the $k$-mer size, sketch size, and containment compress ratio, for the appending genome set `B` are the same as those of the pre-generated sketch `A_sketch`. Additionally, the automatic parameter generation method, which is carried out using the `tune_parameters()` function, depends on whole genome information such as minimum, maximum, and mean genome length.
    Therefore, the changes in the genome lengths of the appending set `B` are unlikely to have a significant effect on the automatically generated parameters if they are not substantial.

* In the context of genome clustering, the sketches are sorted by unstable sort in a decreasing order of their genome length. Consequently, the order of sketches may undergo slight changes if there are genomes with identical lengths. However, this does not significantly affect the outcome of the clustering process.

## [`v.2.1.0`](v.2.1.0.md)
* change the parameter parsing by [CLI11](https://github.com/CLIUtils/CLI11).
* save the intermediate files (sketch, mst files) in binary format.
* abrogate the `-f` option for loading pre-generated sketch or MST file, replaced by `--presketched` and `--premsted` option.

More details by `clust-mst --help` or `clust-greedy --help`.

## `v.2.0.3`
* add the parameter `-m` to set the minimum genome length (*minLen*), genomes with lengths less than *minLen* will be ignored.

## `v.2.0.2`
* update the `calSize` of gz files for automatically generating $k$-mer size .

## `v.2.0.1`
* Update the latest version of [robin-hood-hashing](https://github.com/martinus/robin-hood-hashing) to solve the compile error with `g++ 12.0+`.

## [`v.2.0.0`](v.2.0.0.md)
* Add the `clust-greedy` module for greedy incremental clustering. 
* Last MST-based clustering module is `clust-mst` module.


## `v.1.0.0`
* First version of RabbitTClust, large-scaled genome clustering tool based on sketch technique and Minimum Spanning Tree (MST).