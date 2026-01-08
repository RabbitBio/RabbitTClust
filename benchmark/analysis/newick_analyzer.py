#!/usr/bin/env python3
"""
Newick Tree Analyzer - A complete phylogenetic tree analysis utility

Features:
  - Basic statistics
  - Pairwise distance computation
  - Nearest-neighbor search
  - Subtree extraction
  - Format conversion
  - Distance matrix generation
  - ASCII tree visualization
"""

import sys
import argparse
from Bio import Phylo
from io import StringIO
import os

__version__ = "1.0.0"


class NewickAnalyzer:
    """Newick tree analyzer"""

    def __init__(self, newick_file):
        """Initialize the analyzer."""
        self.newick_file = newick_file
        self.tree = None
        self.terminals = None
        self.nonterminals = None

    def load_tree(self):
        """Load the tree from a Newick file."""
        try:
            self.tree = Phylo.read(self.newick_file, "newick")
            self.terminals = self.tree.get_terminals()
            self.nonterminals = self.tree.get_nonterminals()
            return True
        except Exception as e:
            print(f"Error: failed to read file {self.newick_file}")
            print(f"Details: {e}")
            return False

    def basic_stats(self):
        """Print basic statistics."""
        print("\n" + "=" * 70)
        print("Basic Statistics")
        print("=" * 70)

        num_leaves = len(self.terminals)
        num_internal = len(self.nonterminals)

        print(f"Number of leaf nodes (genomes/sequences): {num_leaves}")
        print(f"Number of internal nodes: {num_internal}")
        print(f"Total nodes: {num_leaves + num_internal}")

        # Tree depth and total length
        depths = self.tree.depths()
        if depths:
            max_depth = max(depths.values())
            min_depth = min(depths.values())
            print(f"\nMaximum depth (root to leaf): {max_depth:.6f}")
            print(f"Minimum depth (root to leaf): {min_depth:.6f}")

        total_length = self.tree.total_branch_length()
        print(f"Total branch length: {total_length:.6f}")

        # Topology
        clade_sizes = [len(clade.clades) for clade in self.nonterminals]
        if clade_sizes:
            avg_children = sum(clade_sizes) / len(clade_sizes)
            is_binary = all(c == 2 for c in clade_sizes)
            print(f"\nTree type: {'Binary' if is_binary else 'Multifurcating'}")
            print(f"Average number of children: {avg_children:.2f}")

        # Branch length stats
        branch_lengths = [
            c.branch_length for c in self.tree.find_clades()
            if c.branch_length is not None
        ]
        if branch_lengths:
            branch_lengths_sorted = sorted(branch_lengths)
            mid = len(branch_lengths_sorted) // 2
            median = branch_lengths_sorted[mid]
            print("\nBranch length statistics:")
            print(f"  Min: {min(branch_lengths):.6f}")
            print(f"  Max: {max(branch_lengths):.6f}")
            print(f"  Mean: {sum(branch_lengths)/len(branch_lengths):.6f}")
            print(f"  Median: {median:.6f}")

    def list_leaves(self, n=20):
        """List leaf node names."""
        print("\n" + "=" * 70)
        print(f"Leaf Node List (first {min(n, len(self.terminals))})")
        print("=" * 70)

        for i, terminal in enumerate(self.terminals[:n]):
            print(f"{i+1:4d}. {terminal.name}")

        if len(self.terminals) > n:
            print(f"\n... {len(self.terminals) - n} more leaf nodes")

    def find_closest_pairs(self, n_pairs=5, sample_size=100):
        """Find the most similar genome pairs (smallest distances)."""
        print("\n" + "=" * 70)
        print(f"Top {n_pairs} Most Similar Genome Pairs")
        print("=" * 70)

        distances = []
        sample = min(sample_size, len(self.terminals))

        print(f"Search scope: first {sample} genomes")

        for i in range(sample):
            for j in range(i + 1, sample):
                dist = self.tree.distance(self.terminals[i], self.terminals[j])
                distances.append((dist, self.terminals[i].name, self.terminals[j].name))

        distances.sort()

        for i, (dist, name1, name2) in enumerate(distances[:n_pairs]):
            print(f"\nPair #{i+1}:")
            print(f"  Genome 1: {name1}")
            print(f"  Genome 2: {name2}")
            print(f"  Distance: {dist:.6f}")

    def find_farthest_pairs(self, n_pairs=5, sample_size=100):
        """Find the most dissimilar genome pairs (largest distances)."""
        print("\n" + "=" * 70)
        print(f"Top {n_pairs} Most Dissimilar Genome Pairs")
        print("=" * 70)

        distances = []
        sample = min(sample_size, len(self.terminals))

        print(f"Search scope: first {sample} genomes")

        for i in range(sample):
            for j in range(i + 1, sample):
                dist = self.tree.distance(self.terminals[i], self.terminals[j])
                distances.append((dist, self.terminals[i].name, self.terminals[j].name))

        distances.sort(reverse=True)

        for i, (dist, name1, name2) in enumerate(distances[:n_pairs]):
            print(f"\nPair #{i+1}:")
            print(f"  Genome 1: {name1}")
            print(f"  Genome 2: {name2}")
            print(f"  Distance: {dist:.6f}")

    def find_neighbors(self, query_name, n_neighbors=10):
        """Find nearest neighbors for a given genome (by name substring match)."""
        print("\n" + "=" * 70)
        print("Nearest Neighbor Search")
        print("=" * 70)

        # Find the query genome
        query = None
        for terminal in self.terminals:
            if query_name in terminal.name:
                query = terminal
                break

        if query is None:
            print(f"Error: no genome name contains '{query_name}'")
            return

        print(f"Query genome: {query.name}")

        # Compute distances
        neighbors = []
        for terminal in self.terminals:
            if terminal != query:
                dist = self.tree.distance(query, terminal)
                neighbors.append((dist, terminal.name))

        neighbors.sort()

        print(f"\nTop {min(n_neighbors, len(neighbors))} nearest neighbors:")
        for i, (dist, name) in enumerate(neighbors[:n_neighbors]):
            print(f"{i+1:3d}. {name} (distance: {dist:.6f})")

    def compute_distance_matrix(self, output_file=None, sample_size=None):
        """Compute a pairwise distance matrix and save it as TSV."""
        print("\n" + "=" * 70)
        print("Distance Matrix")
        print("=" * 70)

        # Determine sample size
        if sample_size is None:
            if len(self.terminals) > 500:
                sample_size = 100
                print(f"Warning: tree is large; computing distance matrix for first {sample_size} genomes only")
            else:
                sample_size = len(self.terminals)

        sample_terminals = self.terminals[:sample_size]

        print(f"Computing {len(sample_terminals)}×{len(sample_terminals)} distance matrix...")

        # Default output filename
        if output_file is None:
            base = os.path.splitext(self.newick_file)[0]
            if sample_size < len(self.terminals):
                output_file = f"{base}_distances_{sample_size}.tsv"
            else:
                output_file = f"{base}_distances.tsv"

        # Compute and write
        with open(output_file, "w") as f:
            # Header
            f.write("Genome\t" + "\t".join(t.name for t in sample_terminals) + "\n")

            # Rows
            for t1 in sample_terminals:
                row = [t1.name]
                for t2 in sample_terminals:
                    if t1 == t2:
                        row.append("0.000000")
                    else:
                        dist = self.tree.distance(t1, t2)
                        row.append(f"{dist:.6f}")
                f.write("\t".join(row) + "\n")

        print(f"✓ Distance matrix saved: {output_file}")
        return output_file

    def extract_subtree(self, genomes, output_file=None):
        """Extract a subtree that contains the specified genomes (by substring match)."""
        print("\n" + "=" * 70)
        print("Subtree Extraction")
        print("=" * 70)

        # Find matching genomes
        matched_terminals = []
        for pattern in genomes:
            for terminal in self.terminals:
                if pattern in terminal.name and terminal not in matched_terminals:
                    matched_terminals.append(terminal)

        if not matched_terminals:
            print("Error: no genomes matched the given patterns")
            return

        print(f"Matched {len(matched_terminals)} genomes:")
        for t in matched_terminals[:10]:
            print(f"  - {t.name}")
        if len(matched_terminals) > 10:
            print(f"  ... {len(matched_terminals) - 10} more")

        # Find the common ancestor
        common_ancestor = self.tree.common_ancestor(matched_terminals)
        subtree_size = common_ancestor.count_terminals()

        print(f"\nSubtree contains {subtree_size} leaf nodes")

        # Default output filename
        if output_file is None:
            base = os.path.splitext(self.newick_file)[0]
            output_file = f"{base}_subtree_{len(matched_terminals)}.newick"

        # Save subtree
        Phylo.write(common_ancestor, output_file, "newick")
        print(f"✓ Subtree saved: {output_file}")

        # Show ASCII for small trees
        if subtree_size <= 30:
            print("\nSubtree visualization:")
            Phylo.draw_ascii(common_ancestor)

        return output_file

    def convert_format(self, output_format):
        """Convert the tree to another format."""
        print("\n" + "=" * 70)
        print(f"Convert to {output_format.upper()} format")
        print("=" * 70)

        base = os.path.splitext(self.newick_file)[0]

        if output_format == "nexus":
            output_file = f"{base}.nexus"
        elif output_format == "phyloxml":
            output_file = f"{base}.xml"
        elif output_format == "nexml":
            output_file = f"{base}.nexml"
        else:
            print(f"Error: unsupported format '{output_format}'")
            return

        try:
            Phylo.write(self.tree, output_file, output_format)
            print(f"✓ Saved: {output_file}")
            return output_file
        except Exception as e:
            print(f"Error: conversion failed - {e}")
            return None

    def show_ascii_tree(self):
        """Show an ASCII visualization of the tree."""
        print("\n" + "=" * 70)
        print("ASCII Tree Visualization")
        print("=" * 70)

        if len(self.terminals) > 50:
            print(f"Tree is large ({len(self.terminals)} leaves). Consider extracting a subtree first.")
            print("Use: --extract-subtree")
        else:
            print()
            Phylo.draw_ascii(self.tree)

    def pairwise_distances(self, genome1, genome2):
        """Compute the distance between two genomes (by substring match)."""
        print("\n" + "=" * 70)
        print("Pairwise Distance")
        print("=" * 70)

        # Find genomes
        t1 = None
        t2 = None

        for terminal in self.terminals:
            if genome1 in terminal.name:
                t1 = terminal
            if genome2 in terminal.name:
                t2 = terminal

        if t1 is None:
            print(f"Error: no genome name contains '{genome1}'")
            return
        if t2 is None:
            print(f"Error: no genome name contains '{genome2}'")
            return

        dist = self.tree.distance(t1, t2)

        print(f"Genome 1: {t1.name}")
        print(f"Genome 2: {t2.name}")
        print(f"Distance: {dist:.6f}")

    def cluster_by_threshold(self, threshold, output_file=None):
        """
        Cluster genomes using a distance threshold.

        Note:
          This is a simple greedy clustering:
          - pick one remaining genome as a seed
          - attach all genomes with distance < threshold to that seed
          - repeat until all genomes are assigned
        """
        print("\n" + "=" * 70)
        print(f"Threshold-based Clustering (threshold = {threshold})")
        print("=" * 70)

        clusters = []
        remaining = set(self.terminals)

        while remaining:
            query = remaining.pop()
            cluster = [query]

            to_remove = []
            for t in remaining:
                if self.tree.distance(query, t) < threshold:
                    cluster.append(t)
                    to_remove.append(t)

            for t in to_remove:
                remaining.remove(t)

            clusters.append(cluster)

        print(f"Found {len(clusters)} clusters")

        # Cluster size stats
        cluster_sizes = sorted([len(c) for c in clusters], reverse=True)
        print("\nCluster size distribution:")
        print(f"  Max: {cluster_sizes[0]}")
        print(f"  Min: {cluster_sizes[-1]}")
        print(f"  Mean: {sum(cluster_sizes)/len(cluster_sizes):.1f}")

        print("\nTop 10 largest clusters:")
        for i, size in enumerate(cluster_sizes[:10]):
            print(f"  Cluster {i+1}: {size} genomes")

        # Save results
        if output_file is None:
            base = os.path.splitext(self.newick_file)[0]
            output_file = f"{base}_clusters_t{threshold}.txt"

        with open(output_file, "w") as f:
            for i, cluster in enumerate(clusters):
                f.write(f">Cluster_{i+1} (size={len(cluster)})\n")
                for terminal in cluster:
                    f.write(f"{terminal.name}\n")
                f.write("\n")

        print(f"\n✓ Clustering result saved: {output_file}")
        return output_file


def main():
    parser = argparse.ArgumentParser(
        description="Newick Tree Analyzer - Phylogenetic tree analysis tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic statistics
  %(prog)s input.newick --stats

  # List leaf nodes
  %(prog)s input.newick --list-leaves

  # Find closest genome pairs
  %(prog)s input.newick --closest-pairs 10

  # Find nearest neighbors of a genome
  %(prog)s input.newick --neighbors GCF_000123456

  # Compute distance matrix
  %(prog)s input.newick --distance-matrix

  # Extract a subtree
  %(prog)s input.newick --extract-subtree GCF_000123 GCF_000456

  # Convert format
  %(prog)s input.newick --convert nexus

  # Show ASCII tree
  %(prog)s input.newick --ascii-tree

  # Pairwise distance
  %(prog)s input.newick --pairwise GCF_000123 GCF_000456

  # Threshold-based clustering
  %(prog)s input.newick --cluster 0.05

  # Run a full analysis (excluding distance matrix and clustering)
  %(prog)s input.newick --all
        """,
    )

    # Required
    parser.add_argument("newick_file", help="Input Newick tree file")

    # Analysis options
    parser.add_argument("--stats", action="store_true", help="Show basic statistics")

    parser.add_argument(
        "--list-leaves",
        type=int,
        metavar="N",
        nargs="?",
        const=20,
        help="List the first N leaf nodes (default: 20)",
    )

    parser.add_argument(
        "--closest-pairs",
        type=int,
        metavar="N",
        nargs="?",
        const=5,
        help="Find the N closest genome pairs (default: 5)",
    )

    parser.add_argument(
        "--farthest-pairs",
        type=int,
        metavar="N",
        nargs="?",
        const=5,
        help="Find the N farthest genome pairs (default: 5)",
    )

    parser.add_argument(
        "--neighbors",
        type=str,
        metavar="GENOME",
        help="Find nearest neighbors for a genome (provide full or partial name)",
    )

    parser.add_argument(
        "--n-neighbors",
        type=int,
        default=10,
        help="Number of neighbors to show (default: 10)",
    )

    parser.add_argument(
        "--distance-matrix",
        action="store_true",
        help="Compute and save a distance matrix",
    )

    parser.add_argument(
        "--matrix-sample",
        type=int,
        metavar="N",
        help="Sample size for distance matrix on large trees",
    )

    parser.add_argument(
        "--extract-subtree",
        type=str,
        nargs="+",
        metavar="PATTERN",
        help="Extract a subtree containing specified genomes (multiple patterns allowed)",
    )

    parser.add_argument(
        "--convert",
        type=str,
        choices=["nexus", "phyloxml", "nexml"],
        help="Convert the tree to another format",
    )

    parser.add_argument("--ascii-tree", action="store_true", help="Show ASCII tree")

    parser.add_argument(
        "--pairwise",
        type=str,
        nargs=2,
        metavar=("GENOME1", "GENOME2"),
        help="Compute distance between two genomes",
    )

    parser.add_argument(
        "--cluster",
        type=float,
        metavar="THRESHOLD",
        help="Cluster genomes using a distance threshold",
    )

    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all analyses (excluding distance matrix and clustering)",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output filename (for certain features)",
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    args = parser.parse_args()

    # Check input file
    if not os.path.exists(args.newick_file):
        print(f"Error: file not found: {args.newick_file}")
        sys.exit(1)

    analyzer = NewickAnalyzer(args.newick_file)

    print("\n" + "=" * 70)
    print(f"Newick Tree Analyzer v{__version__}")
    print("=" * 70)
    print(f"Input file: {args.newick_file}")

    if not analyzer.load_tree():
        sys.exit(1)

    print("✓ Tree loaded successfully")

    # No action check
    no_action = not any(
        [
            args.stats,
            args.list_leaves,
            args.closest_pairs,
            args.farthest_pairs,
            args.neighbors,
            args.distance_matrix,
            args.extract_subtree,
            args.convert,
            args.ascii_tree,
            args.pairwise,
            args.cluster,
            args.all,
        ]
    )

    if no_action:
        print("\nTip: Please specify at least one analysis option.")
        print("Use --help to see all options.")
        print("\nQuick start:")
        print(f"  python {sys.argv[0]} {args.newick_file} --stats")
        sys.exit(0)

    # Run analyses
    if args.all or args.stats:
        analyzer.basic_stats()

    if args.all or args.list_leaves:
        n = args.list_leaves if args.list_leaves else 20
        analyzer.list_leaves(n)

    if args.all or args.closest_pairs:
        n = args.closest_pairs if args.closest_pairs else 5
        analyzer.find_closest_pairs(n)

    if args.all or args.farthest_pairs:
        n = args.farthest_pairs if args.farthest_pairs else 5
        analyzer.find_farthest_pairs(n)

    if args.neighbors:
        analyzer.find_neighbors(args.neighbors, args.n_neighbors)

    if args.distance_matrix:
        analyzer.compute_distance_matrix(args.output, args.matrix_sample)

    if args.extract_subtree:
        analyzer.extract_subtree(args.extract_subtree, args.output)

    if args.convert:
        analyzer.convert_format(args.convert)

    if args.all or args.ascii_tree:
        analyzer.show_ascii_tree()

    if args.pairwise:
        analyzer.pairwise_distances(args.pairwise[0], args.pairwise[1])

    if args.cluster:
        analyzer.cluster_by_threshold(args.cluster, args.output)

    print("\n" + "=" * 70)
    print("Done!")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
