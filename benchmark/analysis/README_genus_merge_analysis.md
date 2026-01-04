# Complete Workflow for Genus Merge Analysis

## I. Overall Workflow Overview

```
Cluster file (230.003.cluster)
    ↓
[Step 1] Analyze all genus relationships (analyze_genus_species_relationships.py)
    ↓
Generate: top_genus_pairs.tsv, boundary_conflicts.tsv, minority_outliers.tsv
    ↓
[Step 2] Extract cluster distribution for specific genus pairs (analyze_genus_pair_clusters.py)
    ↓
Generate: {genus1}_{genus2}_cluster_distribution.tsv
    ↓
[Step 3] Generate visualization plots (plot_genus_pair_visualization.py)
    ↓
Generate: {genus1}_{genus2}_visualization.png
```

## II. Detailed Steps

### Step 1: Analyze All Genus Relationships (Discover Candidate Genus Pairs)

**Script**: `scripts/analyze_genus_species_relationships.py`

**Functionality**:
- Parse cluster file and extract all accessions and cluster IDs
- Load genus and species groundtruth files
- Identify mixed-genus/mixed-species clusters
- Analyze genus pair co-occurrence patterns
- Classify as "boundary conflict" and "minority outlier"
- Generate top genus pairs list

**Usage**:
```bash
python3 scripts/analyze_genus_species_relationships.py \
    --cluster-file 230.003.cluster \
    --species-groundtruth <species_groundtruth_file> \
    --genus-groundtruth <genus_groundtruth_file> \
    --output-dir genus_species_analysis_v2
```

**Input File Formats**:
- **Species groundtruth file**: TSV format with columns `assembly_accession`, `species_taxid`, `organism_name`
- **Genus groundtruth file**: TSV format with columns `assembly_accession`, `genus_id`, `organism_name`

**Output Files**:
- `top_genus_pairs.tsv`: List of high-frequency genus pairs
- `boundary_conflicts.tsv`: Boundary conflict type genus pairs
- `minority_outliers.tsv`: Minority outlier type genus pairs
- `cluster_summary.tsv`: Statistical information for all clusters

### Step 2: Analyze Cluster Distribution for Specific Genus Pairs

**Script**: `scripts/analyze_genus_pair_clusters.py`

**Functionality**:
- Analyze the distribution of two specific genera across clusters
- Count the number and ratio of each genus in each cluster
- Identify merged clusters (both genera co-occur)
- Calculate merge type (Balanced merge vs Minority merge)
- Count other_count (presence of other genera)

**Usage**:
```bash
python3 scripts/analyze_genus_pair_clusters.py \
    --cluster-file 230.003.cluster \
    --genus-groundtruth <genus_groundtruth_file> \
    --g1-id 1827 \
    --g2-id 3259750 \
    --g1-name Rhodococcus \
    --g2-name Rhodococcoides \
    --output-dir genus_species_analysis_v2
```

**Input File Format**:
- **Genus groundtruth file**: TSV format with columns `assembly_accession`, `genus_id`, `organism_name`

**Output Files**:
- `{genus1}_{genus2}_cluster_distribution.tsv`: Detailed cluster distribution table
- `{genus1}_{genus2}_cluster_distribution_summary.tsv`: Summary statistics

**Key Fields**:
- `cluster_id`: Cluster ID
- `total_genomes`: Total number of genomes
- `{genus1}_count`, `{genus2}_count`: Count of each genus
- `other_count`: Count of other genera (boundary clarity indicator)
- `is_mixed`: Whether merged (both genera co-occur)
- `merge_type`: Merge type (Balanced merge / Minority merge)

### Step 3: Generate Visualization Plots

**Script**: `scripts/plot_genus_pair_visualization.py`

**Functionality**:
- Read cluster distribution table
- Generate 4 subplots:
  1. Stacked bar chart for merged clusters
  2. Overall genome distribution pie chart
  3. Cluster size distribution histogram
  4. Summary statistics text

**Usage**:
```bash
python3 scripts/plot_genus_pair_visualization.py \
    --input genus_species_analysis_v2/rhodococcus_rhodococcoides_cluster_distribution.tsv \
    --output genus_species_analysis_v2/rhodococcus_rhodococcoides_visualization.png \
    --g1-name Rhodococcus \
    --g2-name Rhodococcoides
```

**Output Files**:
- `{genus1}_{genus2}_visualization.png`: Visualization plot

## III. Batch Processing Multiple Genus Pairs

**Script**: `scripts/batch_analyze_genus_pairs.py`

**Functionality**:
- Batch process multiple genus pairs
- Automatically execute Step 2 and Step 3
- Generate statistical tables and visualization plots for each genus pair

**Usage**:
```bash
python3 scripts/batch_analyze_genus_pairs.py
```

**Configuration**: Modify the `GENUS_PAIRS` list in the script to add genus pairs for analysis

## IV. Complete Example: Analyzing Rhodococcus ↔ Rhodococcoides

```bash
# Step 1: Analyze all genus relationships (if not done yet)
python3 scripts/analyze_genus_species_relationships.py \
    --cluster-file 230.003.cluster \
    --species-groundtruth <species_groundtruth_file> \
    --genus-groundtruth <genus_groundtruth_file> \
    --output-dir genus_species_analysis_v2

# Step 2: Analyze specific genus pair
python3 scripts/analyze_genus_pair_clusters.py \
    --cluster-file 230.003.cluster \
    --genus-groundtruth <genus_groundtruth_file> \
    --g1-id 1827 \
    --g2-id 3259750 \
    --g1-name Rhodococcus \
    --g2-name Rhodococcoides \
    --output-dir genus_species_analysis_v2

# Step 3: Generate visualization
python3 scripts/plot_genus_pair_visualization.py \
    --input genus_species_analysis_v2/rhodococcus_rhodococcoides_cluster_distribution.tsv \
    --output genus_species_analysis_v2/rhodococcus_rhodococcoides_visualization.png \
    --g1-name Rhodococcus \
    --g2-name Rhodococcoides
```

## V. Key Discovery Metrics

### 1. Number of Merged Clusters
- More is better, indicating that two genera co-occur in multiple clusters

### 2. Number of Merged Genomes
- Total number of genomes in merged clusters
- Higher ratio indicates higher degree of merging

### 3. Balanced Merge vs Minority Merge
- **Balanced merge**: Both genera have comparable proportions in the cluster (30%-70%)
- **Minority merge**: One genus dominates (>70%)
- Balanced merge is stronger evidence

### 4. other_count (Boundary Clarity)
- **other_count = 0**: Indicates that these two genera, as a whole, have clear boundaries with other genera
- This is an important new observation angle

### 5. Largest Merged Cluster
- The genome count and ratio of the largest merged cluster
- Especially large balanced merge clusters are strong evidence

## VI. File Dependencies

```
230.003.cluster (input)
    ↓
analyze_genus_species_relationships.py
    ↓
top_genus_pairs.tsv (discover candidate genus pairs)
    ↓
analyze_genus_pair_clusters.py (requires genus_id)
    ↓
{genus1}_{genus2}_cluster_distribution.tsv
    ↓
plot_genus_pair_visualization.py
    ↓
{genus1}_{genus2}_visualization.png (final output)
```

## VII. Method to Find genus_id

```bash
# Find genus_id for a genus from the genus groundtruth file
grep -i "^[^[:space:]]*[[:space:]]*[0-9]*[[:space:]]*GenusName_" \
    <genus_groundtruth_file> | \
    head -3 | cut -f1-3
```

This will output: `assembly_accession`, `genus_id`, `organism_name` for the first few matches.

## VIII. Input File Formats

### Cluster File
- Format: RabbitTClust output format
- Contains: Cluster definitions with assembly accessions (e.g., `GCF_000010105.1`)

### Species Groundtruth File
- Format: TSV (tab-separated values)
- Required columns:
  - `assembly_accession`: Assembly accession (e.g., `GCF_000010105.1`)
  - `species_taxid`: Species taxonomy ID (integer)
  - `organism_name`: Organism name (e.g., `Rhodococcus_erythropolis_PR4`)

**Example**:
```
assembly_accession	species_taxid	organism_name
GCF_000001215.4	7227	Drosophila_melanogaster
GCF_000001405.40	9606	Homo_sapiens
GCF_000001635.27	10090	Mus_musculus
GCF_000001735.4	3702	Arabidopsis_thaliana
```

### Genus Groundtruth File
- Format: TSV (tab-separated values)
- Required columns:
  - `assembly_accession`: Assembly accession (e.g., `GCF_000010105.1`)
  - `genus_id`: Genus taxonomy ID (integer)
  - `organism_name`: Organism name (e.g., `Rhodococcus_erythropolis_PR4`)

**Example**:
```
assembly_accession	genus_id	organism_name
GCF_000001215.4	7215	Drosophila_melanogaster
GCF_000001405.40	9605	Homo_sapiens
GCF_000001635.27	10088	Mus_musculus
GCF_000001735.4	3701	Arabidopsis_thaliana
```

## IX. Notes

1. **Cluster File Format**: Ensure the cluster file format is correct and contains accession information
2. **Groundtruth Files**: Ensure groundtruth file paths are correct and in TSV format with required columns
3. **Genus ID**: Use correct genus_id, can be found via grep (see Section VII)
4. **Output Directory**: Ensure output directory exists or the script will create it automatically
5. **Memory**: Large-scale analysis may require significant memory (195,631 genomes)

## X. Frequently Asked Questions

### Q: How to find interesting genus pairs?
A: Run Step 1 first, then check `top_genus_pairs.tsv` and `boundary_conflicts.tsv`

### Q: Why is other_count 0?
A: This is a good sign! It indicates that these two genera, as a whole, have clear boundaries with other genera

### Q: How to determine if genera should be merged?
A: Look at the number of merged clusters, merge ratio, balanced merge situation, and other_count

### Q: Can I analyze species-level merging?
A: Yes, modify the script to use species_taxid instead of genus_id
