#!/usr/bin/env python3
"""
Generic visualization script for genus pair relationship across clusters.
"""
import argparse
import csv
import matplotlib.pyplot as plt
import numpy as np


def load_distribution_table(table_file):
    """Load the cluster distribution table."""
    clusters = []
    with open(table_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            clusters.append(row)
    return clusters


def main():
    parser = argparse.ArgumentParser(description='Visualize genus pair relationship')
    parser.add_argument('--input', required=True, help='Input distribution table TSV')
    parser.add_argument('--output', required=True, help='Output PNG file')
    parser.add_argument('--g1-name', required=True, help='Genus 1 name')
    parser.add_argument('--g2-name', required=True, help='Genus 2 name')
    args = parser.parse_args()
    
    clusters = load_distribution_table(args.input)
    
    # Get column names (they should be lowercase genus names)
    g1_col = f"{args.g1_name.lower()}_count"
    g2_col = f"{args.g2_name.lower()}_count"
    
    # Separate merged and non-merged clusters
    merged = [c for c in clusters if c.get('is_mixed', '').lower() == 'true']
    g1_only = [c for c in clusters if c.get('is_mixed', '').lower() == 'false' and int(c.get(g1_col, 0)) > 0]
    g2_only = [c for c in clusters if c.get('is_mixed', '').lower() == 'false' and int(c.get(g2_col, 0)) > 0]
    
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # 1. Stacked bar chart for merged clusters
    ax2 = fig.add_subplot(gs[0, 0])
    
    if merged:
        merged_sorted = sorted(merged, key=lambda x: int(x.get('total_genomes', 0)), reverse=True)
        cluster_ids = [c['cluster_id'] for c in merged_sorted]
        g1_counts = [int(c.get(g1_col, 0)) for c in merged_sorted]
        g2_counts = [int(c.get(g2_col, 0)) for c in merged_sorted]
        
        x_pos = np.arange(len(cluster_ids))
        width = 0.6
        
        bars1 = ax2.bar(x_pos, g1_counts, width, label=args.g1_name, color='#3498db', edgecolor='black', linewidth=0.5)
        bars2 = ax2.bar(x_pos, g2_counts, width, bottom=g1_counts, label=args.g2_name, 
                       color='#9b59b6', edgecolor='black', linewidth=0.5)
        
        # Add value labels
        for i, (g1, g2) in enumerate(zip(g1_counts, g2_counts)):
            total = g1 + g2
            if g1 > 0:
                ax2.text(i, g1/2, str(g1), ha='center', va='center', fontsize=8, fontweight='bold', color='white')
            if g2 > 0:
                ax2.text(i, g1 + g2/2, str(g2), ha='center', va='center', fontsize=8, fontweight='bold', color='white')
            ax2.text(i, total + 0.5, f"n={total}", ha='center', va='bottom', fontsize=7)
        
        ax2.set_xlabel('Cluster ID', fontsize=11, fontweight='bold')
        ax2.set_ylabel('Number of Genomes', fontsize=11, fontweight='bold')
        ax2.set_title('Merged Clusters: Composition Breakdown', fontsize=12, fontweight='bold')
        ax2.set_xticks(x_pos)
        ax2.set_xticklabels(cluster_ids, rotation=45, ha='right')
        ax2.legend(loc='upper right', fontsize=9)
        ax2.grid(axis='y', alpha=0.3)
    else:
        ax2.text(0.5, 0.5, 'No merged clusters found', ha='center', va='center', 
                transform=ax2.transAxes, fontsize=12)
        ax2.set_title('Merged Clusters: Composition Breakdown', fontsize=12, fontweight='bold')
    
    # 2. Pie chart: Overall distribution
    ax3 = fig.add_subplot(gs[0, 1])
    
    total_g1 = sum(int(c.get(g1_col, 0)) for c in clusters)
    total_g2 = sum(int(c.get(g2_col, 0)) for c in clusters)
    total_other = sum(int(c.get('other_count', 0)) for c in clusters)
    
    sizes = [total_g1, total_g2, total_other]
    labels = [args.g1_name, args.g2_name, 'Other']
    colors_pie = ['#3498db', '#9b59b6', '#95a5a6']
    explode = (0.05, 0.1, 0)
    
    if sum(sizes) > 0:
        wedges, texts, autotexts = ax3.pie(sizes, labels=labels, colors=colors_pie, autopct='%1.1f%%',
                                           explode=explode, shadow=True, startangle=90,
                                           textprops={'fontsize': 10, 'fontweight': 'bold'})
    else:
        ax3.text(0.5, 0.5, 'No data', ha='center', va='center', 
                transform=ax3.transAxes, fontsize=12)
    
    ax3.set_title('Overall Genome Distribution\nin Relevant Clusters', fontsize=12, fontweight='bold')
    
    # 3. Cluster size distribution
    ax4 = fig.add_subplot(gs[1, 0])
    
    merged_sizes = [int(c.get('total_genomes', 0)) for c in merged]
    g1_only_sizes = [int(c.get('total_genomes', 0)) for c in g1_only]
    g2_only_sizes = [int(c.get('total_genomes', 0)) for c in g2_only]
    
    all_sizes = merged_sizes + g1_only_sizes + g2_only_sizes
    if all_sizes:
        bins = np.arange(0, max(all_sizes) + 5, 5)
        ax4.hist([merged_sizes, g1_only_sizes, g2_only_sizes], bins=bins, 
                label=['Merged', f'{args.g1_name} only', f'{args.g2_name} only'],
                color=['#e74c3c', '#3498db', '#9b59b6'], alpha=0.7, edgecolor='black', linewidth=0.5)
    else:
        ax4.text(0.5, 0.5, 'No data', ha='center', va='center', 
                transform=ax4.transAxes, fontsize=12)
    
    ax4.set_xlabel('Cluster Size (number of genomes)', fontsize=11, fontweight='bold')
    ax4.set_ylabel('Number of Clusters', fontsize=11, fontweight='bold')
    ax4.set_title('Cluster Size Distribution', fontsize=12, fontweight='bold')
    ax4.legend(loc='upper right', fontsize=9)
    ax4.grid(axis='y', alpha=0.3)
    
    # 4. Summary statistics text
    ax5 = fig.add_subplot(gs[1, 1])
    ax5.axis('off')
    
    total_clusters = len(clusters)
    merged_count = len(merged)
    total_genomes_merged = sum(int(c.get('total_genomes', 0)) for c in merged)
    balanced_count = sum(1 for c in merged if c.get('merge_type', '') == 'Balanced merge')
    
    total_g1_merged = sum(int(c.get(g1_col, 0)) for c in merged)
    total_g2_merged = sum(int(c.get(g2_col, 0)) for c in merged)
    
    if merged:
        largest_balanced = max([c for c in merged if c.get('merge_type', '') == 'Balanced merge'], 
                              key=lambda x: int(x.get('total_genomes', 0)), 
                              default={'cluster_id': 'N/A', 'total_genomes': '0'})
        largest_cluster_id = largest_balanced['cluster_id']
        largest_cluster_size = largest_balanced.get('total_genomes', '0')
    else:
        largest_cluster_id = 'N/A'
        largest_cluster_size = 0
    
    stats_text = f"""
    SUMMARY STATISTICS
    
    Total Clusters: {total_clusters}
    ├─ Merged Clusters: {merged_count} ({merged_count/total_clusters*100:.1f}% if {total_clusters} > 0 else 0)
    │  ├─ Balanced Merges: {balanced_count}
    │  └─ Minority Merges: {merged_count - balanced_count}
    ├─ {args.g1_name} Only: {len(g1_only)}
    └─ {args.g2_name} Only: {len(g2_only)}
    
    Total Genomes in Merged Clusters: {total_genomes_merged}
    ├─ {args.g1_name}: {total_g1_merged} ({total_g1_merged/total_genomes_merged*100:.1f}% if {total_genomes_merged} > 0 else 0)
    └─ {args.g2_name}: {total_g2_merged} ({total_g2_merged/total_genomes_merged*100:.1f}% if {total_genomes_merged} > 0 else 0)
    
    Key Finding:
    Largest balanced merge: Cluster {largest_cluster_id}
    ({largest_cluster_size} genomes)
    """
    
    # Fix percentage calculations
    stats_text = stats_text.replace(f'({merged_count/total_clusters*100:.1f}% if {total_clusters} > 0 else 0)', 
                                    f'{merged_count/total_clusters*100:.1f}%' if total_clusters > 0 else '0%')
    stats_text = stats_text.replace(f'({total_g1_merged/total_genomes_merged*100:.1f}% if {total_genomes_merged} > 0 else 0)', 
                                    f'{total_g1_merged/total_genomes_merged*100:.1f}%' if total_genomes_merged > 0 else '0%')
    stats_text = stats_text.replace(f'({total_g2_merged/total_genomes_merged*100:.1f}% if {total_genomes_merged} > 0 else 0)', 
                                    f'{total_g2_merged/total_genomes_merged*100:.1f}%' if total_genomes_merged > 0 else '0%')
    
    ax5.text(0.1, 0.9, stats_text, transform=ax5.transAxes, fontsize=10,
            verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
    
    plt.suptitle(f'{args.g1_name} and {args.g2_name} Relationship Analysis', 
                fontsize=16, fontweight='bold', y=0.995)
    
    plt.savefig(args.output, dpi=200, bbox_inches='tight', facecolor='white')
    print(f"Visualization saved to: {args.output}")


if __name__ == '__main__':
    main()

