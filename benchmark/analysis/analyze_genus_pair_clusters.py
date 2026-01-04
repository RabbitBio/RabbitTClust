#!/usr/bin/env python3
"""
Analyze cluster distribution for a pair of genera.
Similar to the Rhodococcus/Rhodococcoides analysis.
"""
import argparse
import re
import csv
from collections import defaultdict
from pathlib import Path


def parse_cluster_file(cluster_file: str) -> dict:
    """Parse cluster file and return {accession -> cluster_id}."""
    re_cluster = re.compile(r'^the cluster\s+(\d+)\s+is:', re.I)
    re_acc = re.compile(r'(GC[AF]_\d+\.\d+)')
    
    acc_to_cluster = {}
    cur_cluster = None
    
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
                acc = ma.group(1)
                acc_to_cluster[acc] = cur_cluster
    
    return acc_to_cluster


def load_genus_groundtruth(genus_file: str, target_genus_ids: set) -> dict:
    """
    Load genus groundtruth.
    Returns:
        acc_to_genus: {accession -> genus_id} for ALL accessions (to count 'other')
        genus_id_to_name: {genus_id -> genus_name} for target genera
    """
    acc_to_genus = {}  # Load ALL accessions, not just target genera
    genus_id_to_name = {}
    
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
                # Load ALL accessions, not just target genera
                acc_to_genus[acc] = genus_id_int
                
                # Only store names for target genera
                if genus_id_int in target_genus_ids and orgname:
                    parts = orgname.replace('_', ' ').split()
                    if parts:
                        genus_name = parts[0]
                        if genus_id_int not in genus_id_to_name:
                            genus_id_to_name[genus_id_int] = genus_name
            except ValueError:
                continue
    
    return acc_to_genus, genus_id_to_name


def analyze_cluster_distribution(acc_to_cluster, acc_to_genus, genus_id_to_name, g1_id, g2_id):
    """Analyze cluster distribution for two genera."""
    # First, find clusters that contain target genera
    target_clusters = set()
    for acc, cluster_id in acc_to_cluster.items():
        if acc in acc_to_genus and acc_to_genus[acc] in {g1_id, g2_id}:
            target_clusters.add(cluster_id)
    
    # Now, get ALL accessions in these clusters (not just target genera)
    cluster_to_accs = defaultdict(list)
    for acc, cluster_id in acc_to_cluster.items():
        if cluster_id in target_clusters:
            cluster_to_accs[cluster_id].append(acc)
    
    # Analyze each cluster
    results = []
    for cluster_id, accs in cluster_to_accs.items():
        g1_count = sum(1 for acc in accs if acc_to_genus.get(acc) == g1_id)
        g2_count = sum(1 for acc in accs if acc_to_genus.get(acc) == g2_id)
        # Other = accessions that are in groundtruth but not in target genera, 
        # OR accessions not in groundtruth at all
        other_count = sum(1 for acc in accs 
                         if acc not in acc_to_genus or acc_to_genus.get(acc) not in {g1_id, g2_id})
        
        total = len(accs)
        if total == 0:
            continue
        
        g1_ratio = g1_count / total
        g2_ratio = g2_count / total
        other_ratio = other_count / total
        
        # Count unique species (simplified - just count unique accessions for now)
        g1_species = set()
        g2_species = set()
        other_species = set()
        
        for acc in accs:
            genus_id = acc_to_genus.get(acc)
            if genus_id == g1_id:
                g1_species.add(acc)  # Simplified - could use species_taxid if available
            elif genus_id == g2_id:
                g2_species.add(acc)
            else:
                other_species.add(acc)
        
        is_mixed = (g1_count > 0 and g2_count > 0)
        
        # Determine merge type
        if is_mixed:
            if g1_ratio >= 0.3 and g2_ratio >= 0.3:
                merge_type = "Balanced merge"
            else:
                merge_type = "Minority merge"
        elif g1_count > 0:
            merge_type = f"{genus_id_to_name.get(g1_id, 'G1')} only"
        elif g2_count > 0:
            merge_type = f"{genus_id_to_name.get(g2_id, 'G2')} only"
        else:
            merge_type = "Other only"
        
        results.append({
            'cluster_id': cluster_id,
            'total_genomes': total,
            'g1_count': g1_count,
            'g2_count': g2_count,
            'other_count': other_count,
            'g1_ratio': g1_ratio,
            'g2_ratio': g2_ratio,
            'other_ratio': other_ratio,
            'g1_species_nuniq': len(g1_species),
            'g2_species_nuniq': len(g2_species),
            'other_species_nuniq': len(other_species),
            'total_species_nuniq': len(g1_species) + len(g2_species) + len(other_species),
            'is_mixed': is_mixed,
            'merge_type': merge_type
        })
    
    return results


def main():
    parser = argparse.ArgumentParser(description='Analyze cluster distribution for genus pair')
    parser.add_argument('--cluster-file', required=True, help='Cluster file (e.g., 230.003.cluster)')
    parser.add_argument('--genus-groundtruth', required=True, help='Genus groundtruth file')
    parser.add_argument('--g1-id', type=int, required=True, help='Genus 1 ID')
    parser.add_argument('--g2-id', type=int, required=True, help='Genus 2 ID')
    parser.add_argument('--g1-name', required=True, help='Genus 1 name')
    parser.add_argument('--g2-name', required=True, help='Genus 2 name')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"Loading cluster file: {args.cluster_file}")
    acc_to_cluster = parse_cluster_file(args.cluster_file)
    print(f"Loaded {len(acc_to_cluster)} accessions from clusters")
    
    print(f"Loading genus groundtruth...")
    target_genus_ids = {args.g1_id, args.g2_id}
    acc_to_genus, genus_id_to_name = load_genus_groundtruth(args.genus_groundtruth, target_genus_ids)
    print(f"Found {len(acc_to_genus)} accessions for target genera")
    
    # Analyze
    print("Analyzing cluster distribution...")
    results = analyze_cluster_distribution(acc_to_cluster, acc_to_genus, genus_id_to_name, 
                                         args.g1_id, args.g2_id)
    
    # Sort by cluster_id
    results.sort(key=lambda x: x['cluster_id'])
    
    # Write detailed results
    output_file = Path(args.output_dir) / f"{args.g1_name.lower()}_{args.g2_name.lower()}_cluster_distribution.tsv"
    g1_name_lower = args.g1_name.lower()
    g2_name_lower = args.g2_name.lower()
    
    fieldnames = [
        'cluster_id', 'total_genomes', 
        f'{g1_name_lower}_count', f'{g2_name_lower}_count', 'other_count',
        f'{g1_name_lower}_ratio', f'{g2_name_lower}_ratio', 'other_ratio',
        f'{g1_name_lower}_species_nuniq', f'{g2_name_lower}_species_nuniq', 'other_species_nuniq', 'total_species_nuniq',
        'is_mixed', 'merge_type'
    ]
    
    with open(output_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        # Write results with proper column names
        for r in results:
            output_row = {
                'cluster_id': r['cluster_id'],
                'total_genomes': r['total_genomes'],
                f'{g1_name_lower}_count': r['g1_count'],
                f'{g2_name_lower}_count': r['g2_count'],
                'other_count': r['other_count'],
                f'{g1_name_lower}_ratio': f"{r['g1_ratio']:.3f}",
                f'{g2_name_lower}_ratio': f"{r['g2_ratio']:.3f}",
                'other_ratio': f"{r['other_ratio']:.3f}",
                f'{g1_name_lower}_species_nuniq': r['g1_species_nuniq'],
                f'{g2_name_lower}_species_nuniq': r['g2_species_nuniq'],
                'other_species_nuniq': r['other_species_nuniq'],
                'total_species_nuniq': r['total_species_nuniq'],
                'is_mixed': str(r['is_mixed']),
                'merge_type': r['merge_type']
            }
            writer.writerow(output_row)
    
    print(f"Detailed results written to: {output_file}")
    
    # Calculate summary statistics
    total_clusters = len(results)
    merged_clusters = [r for r in results if r['is_mixed']]
    g1_only = [r for r in results if r['g1_count'] > 0 and r['g2_count'] == 0]
    g2_only = [r for r in results if r['g2_count'] > 0 and r['g1_count'] == 0]
    
    total_g1 = sum(r['g1_count'] for r in results)
    total_g2 = sum(r['g2_count'] for r in results)
    total_genomes = sum(r['total_genomes'] for r in results)
    
    g1_ratio_overall = total_g1 / total_genomes if total_genomes > 0 else 0
    g2_ratio_overall = total_g2 / total_genomes if total_genomes > 0 else 0
    
    # Write summary
    summary_file = Path(args.output_dir) / f"{args.g1_name.lower()}_{args.g2_name.lower()}_cluster_distribution_summary.tsv"
    with open(summary_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['metric', 'value'])
        writer.writerow([f'Total clusters with {args.g1_name} or {args.g2_name}', total_clusters])
        writer.writerow([f'Clusters with both genera (merged)', len(merged_clusters)])
        writer.writerow([f'Clusters with {args.g1_name} only', len(g1_only)])
        writer.writerow([f'Clusters with {args.g2_name} only', len(g2_only)])
        writer.writerow([f'Total {args.g1_name} genomes', total_g1])
        writer.writerow([f'Total {args.g2_name} genomes', total_g2])
        writer.writerow([f'Total genomes in relevant clusters', total_genomes])
        writer.writerow([f'{args.g1_name} ratio (overall)', f"{g1_ratio_overall:.3f}"])
        writer.writerow([f'{args.g2_name} ratio (overall)', f"{g2_ratio_overall:.3f}"])
    
    print(f"Summary written to: {summary_file}")
    print(f"\nSummary:")
    print(f"  Total clusters: {total_clusters}")
    print(f"  Merged clusters: {len(merged_clusters)}")
    print(f"  {args.g1_name} only: {len(g1_only)}")
    print(f"  {args.g2_name} only: {len(g2_only)}")
    print(f"  Total {args.g1_name} genomes: {total_g1}")
    print(f"  Total {args.g2_name} genomes: {total_g2}")


if __name__ == '__main__':
    main()

