# RabbitTClust v2.4.0

RabbitTClust v2.4.0 expands the original MST and greedy workflows with the
functionality used in the JBHI study. The fixed paper snapshot is available as
the [`v2.4-JBHI`](https://github.com/RabbitBio/RabbitTClust/releases/tag/v2.4-JBHI)
release.

## Highlights

- Accelerated KSSD and MinHash clustering paths based on inverted indexes.
- Fast and incremental workflows for MST and greedy clustering, including
  pre-sketched input and append operations.
- DBSCAN clustering and representative/MST database workflows.
- Additional tree and matrix outputs, including Newick, PHYLIP, NEXUS, and
  linkage-matrix formats.
- Cluster post-processing for near-duplicate collapse and representative
  selection.
- Benchmark download, simulation, evaluation, and analysis scripts under
  [`benchmark/`](../benchmark/), with study result tables under
  [`discovery/`](../discovery/).
