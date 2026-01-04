#!/usr/bin/env python3
"""
Analyze genus and species relationships in RabbitTClust clusters.

This script:
1. Maps assembly accessions to NCBI ground truth (genus_id, species_taxid)
2. Identifies mixed clusters (multiple genera/species)
3. Calculates purity, nuniq, majority labels
4. Marks suspects (members inconsistent with majority)
5. Analyzes genus/species co-occurrence patterns
6. Classifies as "boundary conflict" vs "minority outlier" types
7. Outputs key clusters and candidate genomes
"""
import argparse
import re
import collections
import csv
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional


def parse_cluster_file(cluster_file: str) -> List[Tuple[str, int]]:
    """Parse RabbitTClust cluster file and return [(accession, cluster_id), ...]"""
    re_cluster = re.compile(r'^the cluster\s+(\d+)\s+is:', re.I)
    re_acc = re.compile(r'(GC[AF]_\d+\.\d+)')
    
    cur_cluster = None
    acc_cluster = []
    with open(cluster_file, 'r', errors='ignore') as f:
        for line in f:
            s = line.strip()
            m = re_cluster.match(s)
            if m:
                cur_cluster = int(m.group(1))
                continue
            if cur_cluster is None or not s:
                continue
            ma = re_acc.search(line)
            if ma:
                acc_cluster.append((ma.group(1), cur_cluster))
    return acc_cluster


def load_groundtruth(species_file: str, genus_file: str = None) -> Tuple[Dict[str, int], Dict[str, str], Dict[str, int], Dict[int, int], Dict[int, str], Dict[int, str]]:
    """
    Load groundtruth files.
    Args:
        species_file: File with assembly_accession, species_taxid, organism_name
        genus_file: File with assembly_accession, genus_id, organism_name (optional)
    Returns:
        acc_to_species: {accession -> species_taxid}
        acc_to_orgname: {accession -> organism_name}
        acc_to_genus: {accession -> genus_id}
        species_to_genus: {species_taxid -> genus_id} (derived from accession mapping)
        genus_id_to_name: {genus_id -> genus_name}
        species_id_to_name: {species_taxid -> species_name}
    """
    acc_to_species = {}
    acc_to_orgname = {}
    acc_to_genus = {}
    genus_id_to_name = {}
    species_id_to_name = {}
    
    # Load species file
    with open(species_file, 'r', errors='ignore') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            acc = row.get('assembly_accession', '').strip()
            species_taxid = row.get('species_taxid', '').strip()
            orgname = row.get('organism_name', '').strip()
            
            if not acc or not species_taxid:
                continue
            
            try:
                species_taxid_int = int(species_taxid)
                acc_to_species[acc] = species_taxid_int
                acc_to_orgname[acc] = orgname
                
                # Extract species name from organism_name (usually first two words)
                if orgname:
                    parts = orgname.replace('_', ' ').split()
                    if len(parts) >= 2:
                        species_name = ' '.join(parts[:2])  # Genus species
                    else:
                        species_name = parts[0] if parts else orgname
                    species_id_to_name[species_taxid_int] = species_name
            except ValueError:
                continue
    
    # Load genus file if provided
    if genus_file and Path(genus_file).exists():
        with open(genus_file, 'r', errors='ignore') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                acc = row.get('assembly_accession', '').strip()
                genus_id = row.get('genus_id', '').strip()
                orgname = row.get('organism_name', '').strip()
                
                if not acc or not genus_id:
                    continue
                
                try:
                    genus_id_int = int(genus_id)
                    acc_to_genus[acc] = genus_id_int
                    
                    # Extract genus name from organism_name (first word)
                    if orgname:
                        parts = orgname.replace('_', ' ').split()
                        if parts:
                            genus_name = parts[0]
                            genus_id_to_name[genus_id_int] = genus_name
                except ValueError:
                    continue
    
    # Build species_to_genus mapping (from accession-level mapping)
    species_to_genus = {}
    for acc, species_taxid in acc_to_species.items():
        if acc in acc_to_genus:
            species_to_genus[species_taxid] = acc_to_genus[acc]
    
    return acc_to_species, acc_to_orgname, acc_to_genus, species_to_genus, genus_id_to_name, species_id_to_name


def analyze_clusters(acc_cluster: List[Tuple[str, int]], 
                    acc_to_species: Dict[str, int],
                    acc_to_orgname: Dict[str, str],
                    acc_to_genus: Dict[str, int]) -> Dict:
    """
    Analyze clusters and return statistics.
    Returns dict with cluster-level and genome-level information.
    """
    # Build cluster membership
    cluster_members = collections.defaultdict(list)
    for acc, cid in acc_cluster:
        if acc in acc_to_species:
            cluster_members[cid].append(acc)
    
    cluster_stats = {}
    all_genomes = []
    
    for cid, members in cluster_members.items():
        # Count genera and species
        genus_counts = collections.Counter()
        species_counts = collections.Counter()
        genus_species_map = collections.defaultdict(set)
        
        for acc in members:
            species_taxid = acc_to_species.get(acc)
            if not species_taxid:
                continue
            genus_id = acc_to_genus.get(acc)
            if genus_id:
                genus_counts[genus_id] += 1
                genus_species_map[genus_id].add(species_taxid)
            species_counts[species_taxid] += 1
        
        # Calculate metrics
        cluster_size = len(members)
        genus_nuniq = len(genus_counts)
        species_nuniq = len(species_counts)
        
        # Majority labels
        majority_genus = genus_counts.most_common(1)[0][0] if genus_counts else None
        majority_species = species_counts.most_common(1)[0][0] if species_counts else None
        
        # Purity (fraction of majority label)
        genus_purity = (genus_counts[majority_genus] / cluster_size) if majority_genus else 0.0
        species_purity = (species_counts[majority_species] / cluster_size) if majority_species else 0.0
        
        # Identify mixed clusters
        is_mixed_genus = genus_nuniq > 1
        is_mixed_species = species_nuniq > 1
        
        # Mark suspects (inconsistent with majority)
        suspects = []
        for acc in members:
            species_taxid = acc_to_species.get(acc)
            if not species_taxid:
                continue
            genus_id = acc_to_genus.get(acc)
            
            is_suspect = False
            if majority_genus and genus_id != majority_genus:
                is_suspect = True
            if majority_species and species_taxid != majority_species:
                is_suspect = True
            
            if is_suspect:
                suspects.append({
                    'accession': acc,
                    'genus_id': genus_id,
                    'species_taxid': species_taxid,
                    'organism_name': acc_to_orgname.get(acc, ''),
                    'cluster_id': cid
                })
        
        cluster_stats[cid] = {
            'cluster_id': cid,
            'cluster_size': cluster_size,
            'genus_nuniq': genus_nuniq,
            'species_nuniq': species_nuniq,
            'genus_counts': dict(genus_counts),
            'species_counts': dict(species_counts),
            'majority_genus': majority_genus,
            'majority_species': majority_species,
            'genus_purity': genus_purity,
            'species_purity': species_purity,
            'is_mixed_genus': is_mixed_genus,
            'is_mixed_species': is_mixed_species,
            'suspects': suspects
        }
        
        # Store genome-level info
        for acc in members:
            all_genomes.append({
                'accession': acc,
                'cluster_id': cid,
                'genus_id': acc_to_genus.get(acc),
                'species_taxid': acc_to_species.get(acc),
                'organism_name': acc_to_orgname.get(acc, ''),
                'is_suspect': acc in [s['accession'] for s in suspects]
            })
    
    return {
        'cluster_stats': cluster_stats,
        'genomes': all_genomes
    }


def analyze_cooccurrence(cluster_stats: Dict) -> Dict[Tuple[int, int], List[Dict]]:
    """
    Analyze genus/species co-occurrence patterns.
    Returns: {(g1, g2): [cluster_info, ...]}
    """
    cooccurrence = collections.defaultdict(list)
    
    for cid, stats in cluster_stats.items():
        if not stats['is_mixed_genus']:
            continue
        
        # Get all genus pairs in this cluster
        genus_list = list(stats['genus_counts'].keys())
        for i, g1 in enumerate(genus_list):
            for g2 in genus_list[i+1:]:
                pair = tuple(sorted([g1, g2]))
                g1_count = stats['genus_counts'][g1]
                g2_count = stats['genus_counts'][g2]
                g1_ratio = g1_count / stats['cluster_size']
                g2_ratio = g2_count / stats['cluster_size']
                
                cooccurrence[pair].append({
                    'cluster_id': cid,
                    'cluster_size': stats['cluster_size'],
                    'g1_count': g1_count,
                    'g2_count': g2_count,
                    'g1_ratio': g1_ratio,
                    'g2_ratio': g2_ratio,
                    'species_nuniq': stats['species_nuniq'],
                    'genus_nuniq': stats['genus_nuniq'],
                    'genus_purity': stats['genus_purity'],
                    'species_purity': stats['species_purity']
                })
    
    return cooccurrence


def classify_cooccurrence(cooccurrence: Dict[Tuple[int, int], List[Dict]], 
                         threshold_balanced: float = 0.3,
                         threshold_clean: float = 0.7) -> Dict:
    """
    Classify co-occurrence patterns into:
    - "boundary_conflict": both genera have high ratio and cluster is relatively clean
    - "minority_outlier": one genus dominates, the other is minority
    """
    classified = {
        'boundary_conflict': [],
        'minority_outlier': []
    }
    
    for (g1, g2), clusters in cooccurrence.items():
        for cluster_info in clusters:
            g1_ratio = cluster_info['g1_ratio']
            g2_ratio = cluster_info['g2_ratio']
            min_ratio = min(g1_ratio, g2_ratio)
            max_ratio = max(g1_ratio, g2_ratio)
            
            # Boundary conflict: both have significant presence and cluster is clean
            if (min_ratio >= threshold_balanced and 
                cluster_info['genus_purity'] < threshold_clean and
                cluster_info['cluster_size'] >= 10):
                classified['boundary_conflict'].append({
                    'g1': g1,
                    'g2': g2,
                    **cluster_info
                })
            # Minority outlier: one dominates
            elif max_ratio >= 0.7 and min_ratio < 0.3:
                classified['minority_outlier'].append({
                    'g1': g1,
                    'g2': g2,
                    **cluster_info
                })
    
    return classified


def main():
    parser = argparse.ArgumentParser(description='Analyze genus/species relationships in RabbitTClust clusters')
    parser.add_argument('--cluster', required=True, help='RabbitTClust cluster file')
    parser.add_argument('--species-groundtruth', required=True, help='NCBI species groundtruth file (assembly_accession, species_taxid, organism_name)')
    parser.add_argument('--genus-groundtruth', required=True, help='NCBI genus groundtruth file (assembly_accession, genus_id, organism_name)')
    parser.add_argument('--top-k', type=int, default=20, help='Top-K genus pairs to analyze')
    parser.add_argument('--output-dir', default='.', help='Output directory for results')
    args = parser.parse_args()
    
    print("Loading data...")
    acc_cluster = parse_cluster_file(args.cluster)
    acc_to_species, acc_to_orgname, acc_to_genus, species_to_genus, genus_id_to_name, species_id_to_name = load_groundtruth(
        args.species_groundtruth, args.genus_groundtruth
    )
    
    print(f"Loaded {len(acc_cluster)} accessions from {args.cluster}")
    print(f"Loaded {len(acc_to_species)} accessions with species groundtruth")
    print(f"Loaded {len(acc_to_genus)} accessions with genus groundtruth")
    print(f"Loaded {len(genus_id_to_name)} genus name mappings")
    print(f"Loaded {len(species_id_to_name)} species name mappings")
    
    print("Analyzing clusters...")
    results = analyze_clusters(acc_cluster, acc_to_species, acc_to_orgname, acc_to_genus)
    
    # Count mixed clusters
    mixed_genus = sum(1 for s in results['cluster_stats'].values() if s['is_mixed_genus'])
    mixed_species = sum(1 for s in results['cluster_stats'].values() if s['is_mixed_species'])
    print(f"Found {mixed_genus} mixed-genus clusters and {mixed_species} mixed-species clusters")
    
    print("Analyzing co-occurrence patterns...")
    cooccurrence = analyze_cooccurrence(results['cluster_stats'])
    print(f"Found {len(cooccurrence)} unique genus pairs")
    
    # Get top-K by frequency
    top_pairs = sorted(cooccurrence.items(), key=lambda x: len(x[1]), reverse=True)[:args.top_k]
    
    print("Classifying patterns...")
    classified = classify_cooccurrence(cooccurrence)
    
    # Output results
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Top-K genus pairs summary
    with open(output_dir / 'top_genus_pairs.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['g1', 'g1_name', 'g2', 'g2_name', 'cluster_id', 'cluster_size', 'g1_count', 'g2_count', 
                        'g1_ratio', 'g2_ratio', 'species_nuniq', 'score'])
        for (g1, g2), clusters in top_pairs:
            g1_name = genus_id_to_name.get(g1, f'genus_{g1}')
            g2_name = genus_id_to_name.get(g2, f'genus_{g2}')
            for cluster_info in clusters:
                score = min(cluster_info['g1_count'], cluster_info['g2_count'])
                writer.writerow([
                    g1, g1_name, g2, g2_name, cluster_info['cluster_id'], cluster_info['cluster_size'],
                    cluster_info['g1_count'], cluster_info['g2_count'],
                    f"{cluster_info['g1_ratio']:.3f}", f"{cluster_info['g2_ratio']:.3f}",
                    cluster_info['species_nuniq'], score
                ])
    
    # 2. Classified patterns
    with open(output_dir / 'boundary_conflicts.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['g1', 'g1_name', 'g2', 'g2_name', 'cluster_id', 'cluster_size', 'g1_count', 'g2_count',
                        'g1_ratio', 'g2_ratio', 'species_nuniq', 'genus_purity'])
        for item in classified['boundary_conflict']:
            g1_name = genus_id_to_name.get(item['g1'], f'genus_{item["g1"]}')
            g2_name = genus_id_to_name.get(item['g2'], f'genus_{item["g2"]}')
            writer.writerow([
                item['g1'], g1_name, item['g2'], g2_name, item['cluster_id'], item['cluster_size'],
                item['g1_count'], item['g2_count'],
                f"{item['g1_ratio']:.3f}", f"{item['g2_ratio']:.3f}",
                item['species_nuniq'], f"{item['genus_purity']:.3f}"
            ])
    
    with open(output_dir / 'minority_outliers.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['g1', 'g1_name', 'g2', 'g2_name', 'cluster_id', 'cluster_size', 'g1_count', 'g2_count',
                        'g1_ratio', 'g2_ratio', 'species_nuniq'])
        for item in classified['minority_outlier']:
            g1_name = genus_id_to_name.get(item['g1'], f'genus_{item["g1"]}')
            g2_name = genus_id_to_name.get(item['g2'], f'genus_{item["g2"]}')
            writer.writerow([
                item['g1'], g1_name, item['g2'], g2_name, item['cluster_id'], item['cluster_size'],
                item['g1_count'], item['g2_count'],
                f"{item['g1_ratio']:.3f}", f"{item['g2_ratio']:.3f}",
                item['species_nuniq']
            ])
    
    # 3. Suspects (potential mislabels/contamination)
    all_suspects = []
    for stats in results['cluster_stats'].values():
        all_suspects.extend(stats['suspects'])
    
    with open(output_dir / 'suspects.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['accession', 'cluster_id', 'genus_id', 'genus_name', 'species_taxid', 'species_name', 'organism_name'])
        for suspect in all_suspects:
            genus_name = genus_id_to_name.get(suspect['genus_id'], f'genus_{suspect["genus_id"]}') if suspect['genus_id'] else 'Unknown'
            species_name = species_id_to_name.get(suspect['species_taxid'], f'species_{suspect["species_taxid"]}') if suspect['species_taxid'] else 'Unknown'
            writer.writerow([
                suspect['accession'], suspect['cluster_id'],
                suspect['genus_id'], genus_name, suspect['species_taxid'], species_name, suspect['organism_name']
            ])
    
    # 4. Cluster summary
    with open(output_dir / 'cluster_summary.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['cluster_id', 'cluster_size', 'genus_nuniq', 'species_nuniq',
                        'majority_genus', 'majority_genus_name', 'majority_species', 'majority_species_name',
                        'genus_purity', 'species_purity', 'is_mixed_genus', 'is_mixed_species', 'n_suspects'])
        for stats in sorted(results['cluster_stats'].values(), key=lambda x: x['cluster_id']):
            majority_genus_name = genus_id_to_name.get(stats['majority_genus'], f'genus_{stats["majority_genus"]}') if stats['majority_genus'] else 'Unknown'
            majority_species_name = species_id_to_name.get(stats['majority_species'], f'species_{stats["majority_species"]}') if stats['majority_species'] else 'Unknown'
            writer.writerow([
                stats['cluster_id'], stats['cluster_size'],
                stats['genus_nuniq'], stats['species_nuniq'],
                stats['majority_genus'], majority_genus_name, stats['majority_species'], majority_species_name,
                f"{stats['genus_purity']:.3f}", f"{stats['species_purity']:.3f}",
                stats['is_mixed_genus'], stats['is_mixed_species'],
                len(stats['suspects'])
            ])
    
    print(f"\nResults written to {output_dir}/")
    print(f"  - top_genus_pairs.tsv: Top-{args.top_k} genus pairs")
    print(f"  - boundary_conflicts.tsv: {len(classified['boundary_conflict'])} boundary conflict cases")
    print(f"  - minority_outliers.tsv: {len(classified['minority_outlier'])} minority outlier cases")
    print(f"  - suspects.tsv: {len(all_suspects)} suspect genomes")
    print(f"  - cluster_summary.tsv: Summary of all clusters")


if __name__ == '__main__':
    main()

