# Latest version: `v.2.1.0` 
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