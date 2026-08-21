# FracMinHash containment-SCC on RefSeq Release 233

Runs were performed on 2026-08-15 and 2026-08-16 (Asia/Shanghai) using
`233.list`: 208,753 bacterial assemblies and 945,020,243,458 input bases. All
208,753 accessions mapped to `/home/bigssd/zt/ref/233_species_groundtruth`.

## Final operating point

The final conservative assembly preset is `k=19`, `scaled=10000`, `d=0.02`
(directed containment cutoff 0.6838614092), 128 candidate hashes, 256 exact
candidates, and 64 retained outgoing neighbors. Candidate and edge budgets use
a deterministic 75% high-score / 25% diversity allocation to avoid trapping a
large species inside one densely sampled lineage.

```bash
/usr/bin/time -v -o containment_final_end2end.time \
  /tmp/rabbittclust-containment-build/clust-containment \
  -l -i 233.list -o containment_final_end2end.cluster \
  -t 64 -m 10000 -e
```

This command used only CLI defaults other than threads, input, output, and
`--no-save`. It read every FASTA again rather than loading a sketch cache.

| Result | Value |
| --- | ---: |
| End-to-end wall time | 138.61 s |
| Peak RSS | 3.127 GiB |
| Sketching/input stage | 112.14 s |
| Containment graph stage | 24.34 s |
| SCC/DAG stage | 0.71 s |
| Directed edges | 8,658,292 |
| Strict SCCs | 37,784 |
| Condensation-DAG relations | 3,700 |
| Precision | 0.8798 |
| Recall | 0.9765 |
| Pairwise F1 | 0.9256 |
| NMI | 0.947360980 |

Loading the 729 MiB native FracMinHash cache and regenerating the same graph
and outputs took 25.74 s wall time with 2.599 GiB peak RSS. The fresh and cached
strict cluster files were byte-identical.

For comparison, the existing fixed-budget Greedy benchmark on the same input
reported BinDash at 104.86 s / 1.565 GiB / F1 0.9314 and FastKMV at 150.83 s /
3.937 GiB / F1 0.9313. Thus containment-SCC is in the same runtime and memory
regime; its strict F1 is 0.0058 below BinDash-Greedy. The permissive weak view
reached F1 0.9264, but it is reported as a diagnostic upper bound rather than
substituted for the primary SCC result.

Doubling sketch density to `scaled=5000` did not improve accuracy: the strict
result was P/R/F1 = 0.8796/0.9777/0.9260, while end-to-end time rose to 254.31 s
and peak RSS to 4.096 GiB. The default therefore remains `scaled=10000`.

## What the directed graph adds

The strict flat metric does not capture the main new output. The final graph
contains:

- 3,700 SCC-to-SCC containment relations: 2,350 within one species taxid and
  1,350 across taxids.
- 1,134 strongly asymmetric relations with forward-minus-reverse containment
  at least 0.1.
- 330 short-to-long relations whose source/target assembly-length ratio is at
  most 0.8.
- 37,276 of 37,784 SCCs (98.7%) reaching one unique container sink, covering
  197,197 genomes (94.5%). The remaining 508 SCCs, containing 11,556 genomes,
  reach multiple sinks and form a focused taxonomy/structure-review set.
- Among the 2,335 non-sink SCCs, 1,827 (22,459 genomes) have a unique maximal
  container. Of the 508 ambiguous SCCs, 461 branch directly toward multiple
  containers, rather than becoming ambiguous only through a long path.

Three directly inspectable examples show information lost by ordinary flat
clustering:

| Source -> container | Length (Mb) | FASTA records | Forward / reverse containment | Taxonomy and Greedy result |
| --- | ---: | ---: | ---: | --- |
| GCF_001586935.1 -> GCF_900174485.1 | 0.952 -> 3.103 | 3 -> 4 | 1.000 / 0.290 | Both *Sterolibacterium denitrificans*; BinDash-Greedy puts both in cluster 19383 |
| GCF_003177275.1 -> GCF_003568665.1 | 2.664 -> 6.130 | 1 -> 15 | 1.000 / 0.445 | Both *Chitinophaga alhagiae*; BinDash-Greedy puts both in cluster 6354 |
| GCF_018265195.1 -> GCF_000209935.1 | 1.703 -> 3.345 | 198 -> 1 | 0.854 / 0.410 | An *Agathobacter* MAG (taxid 2021311) points to *A. rectalis* (taxid 39491); BinDash-Greedy puts both in cluster 14632 |

The last source has 42 outgoing condensation-DAG relations, all to distinct
*A. rectalis* SCCs and all crossing from taxid 2021311 to taxid 39491. This is a
concrete completeness/taxonomy-review candidate: the directed evidence says
which assembly is the partial object and which references contain it. It is not
by itself proof that the NCBI taxonomic assignment is wrong.

### A genuinely non-tree result: host-associated Tisiphia genome reduction

`GCF_964264975.1`, the complete 1.311 Mb *Candidatus Tisiphia* endosymbiont of
*Neophilaenus lineatus*, has five direct outgoing containment relations to
complete Tisiphia genomes from five other hosts. Forward containment is
0.719--0.833, whereas reverse containment is only 0.462--0.530. Those branches
terminate in three mutually distinct sink taxids. BinDash-Greedy places all six
assemblies in cluster 22069, so the reduction direction and alternative
containers disappear from its result.

An MST may retain one or several of these edges if they help its global tree,
but it necessarily discards non-tree alternatives and supplies no intrinsic
containment orientation. Recovering all five directed parents and the three
incomparable terminal sinks requires returning to the non-tree asymmetric
graph. Biologically, this is a testable candidate for a host-associated genome
reduction/core-retention pattern; gene-level alignment is required before
claiming specific gene losses.

Greedy clustering records only that each example shares a cluster. A symmetric
MST can expose threshold connectivity and length differences, but the directed
containment DAG directly encodes the partial-to-container relation, unique sink,
and multi-sink ambiguity. That is the distinct use case of this module on 233.

## Artifacts

- `containment_final_end2end.cluster`: strict SCC clustering from the no-cache
  terminal run.
- `containment_final_end2end.cluster.{scc,members,dag,assignment,summary}.tsv`:
  component, membership, direction, sink assignment, and timing data.
- `containment_final_s10000_d020.{taxonomy_summary,species_splits,scc_impurity,directed_sources,ambiguous_sinks,dag_annotated}.tsv`:
  taxid-aware analysis. Species is grouped by `species_taxid`, not by the
  strain-bearing `organism_name` string.
- `containment_final_s10000_d020.label`: the two-row label file used for the
  strict clustering metrics.
- `containment_final_end2end.time` and
  `containment_final_s10000_d020_graph.time`: GNU time resource records.
