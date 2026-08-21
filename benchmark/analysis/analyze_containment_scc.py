#!/usr/bin/env python3
"""Annotate containment-SCC output with RefSeq species ground truth."""

import argparse
import csv
import os
import re
from collections import Counter, defaultdict


ACCESSION_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


def accession_from_path(path):
    match = ACCESSION_RE.search(os.path.basename(path))
    return match.group(1) if match else ""


def load_ground_truth(path):
    result = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            result[row["assembly_accession"]] = (
                row["species_taxid"], row["organism_name"]
            )
    return result


def load_members(path, ground_truth):
    nodes = {}
    components = defaultdict(list)
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            node = int(row["node_id"])
            component = int(row["scc_id"])
            accession = accession_from_path(row["path"])
            taxid, name = ground_truth.get(accession, ("NA", "NA"))
            record = {
                "node": node,
                "scc": component,
                "path": row["path"],
                "accession": accession,
                "length": int(row["length"]),
                "hashes": int(row["sketch_hashes"]),
                "taxid": taxid,
                "name": name,
            }
            nodes[node] = record
            components[component].append(record)
    return nodes, components


def species_split_rows(nodes):
    species = defaultdict(list)
    for record in nodes.values():
        species[record["taxid"]].append(record)
    rows = []
    for taxid, records in species.items():
        # RefSeq organism_name commonly includes a strain suffix, so taxid is
        # the species label and the most frequent text is only a display name.
        name = Counter(record["name"] for record in records).most_common(1)[0][0]
        counts = Counter(record["scc"] for record in records)
        largest_scc, largest = counts.most_common(1)[0]
        total = len(records)
        rows.append({
            "species_taxid": taxid,
            "organism_name": name,
            "genomes": total,
            "scc_count": len(counts),
            "largest_scc": largest_scc,
            "largest_scc_genomes": largest,
            "largest_fraction": largest / total,
            "split_genomes": total - largest,
        })
    rows.sort(key=lambda row: (row["split_genomes"], row["genomes"]),
              reverse=True)
    return rows


def component_impurity_rows(components, taxid_names):
    rows = []
    for component, records in components.items():
        counts = Counter(record["taxid"] for record in records)
        major_taxid, major_count = counts.most_common(1)[0]
        major_name = taxid_names[major_taxid]
        total = len(records)
        rows.append({
            "scc_id": component,
            "size": total,
            "species_count": len(counts),
            "major_species_taxid": major_taxid,
            "major_organism_name": major_name,
            "major_count": major_count,
            "purity": major_count / total,
            "off_species_genomes": total - major_count,
        })
    rows.sort(key=lambda row: (row["off_species_genomes"], row["size"]),
              reverse=True)
    return rows


def annotate_dag(path, nodes):
    rows = []
    if not path:
        return rows
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            source = nodes[int(row["source_node"])]
            target = nodes[int(row["target_node"])]
            annotated = dict(row)
            annotated.update({
                "source_accession": source["accession"],
                "target_accession": target["accession"],
                "source_species_taxid": source["taxid"],
                "target_species_taxid": target["taxid"],
                "source_organism_name": source["name"],
                "target_organism_name": target["name"],
                "same_species": int(source["taxid"] == target["taxid"]),
                "containment_asymmetry": (
                    float(row["containment"]) -
                    float(row["reciprocal_containment"])
                ),
            })
            rows.append(annotated)
    return rows


def load_assignment(path, components):
    if not path:
        return []
    rows = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            component = int(row["scc_id"])
            rows.append({
                "scc_id": component,
                "status": row["status"],
                "genomes": len(components[component]),
            })
    return rows


def directed_source_rows(dag_rows):
    grouped = defaultdict(list)
    for row in dag_rows:
        grouped[row["source_accession"]].append(row)
    result = []
    for accession, rows in grouped.items():
        first = rows[0]
        cross = [row for row in rows if not row["same_species"]]
        result.append({
            "source_accession": accession,
            "source_species_taxid": first["source_species_taxid"],
            "source_organism_name": first["source_organism_name"],
            "source_length": first["source_length"],
            "outgoing_scc_relations": len(rows),
            "cross_species_relations": len(cross),
            "distinct_target_taxids": len({
                row["target_species_taxid"] for row in rows
            }),
            "max_containment_asymmetry": max(
                float(row["containment_asymmetry"]) for row in rows
            ),
            "min_length_ratio": min(float(row["length_ratio"])
                                    for row in rows),
        })
    result.sort(key=lambda row: (
        row["cross_species_relations"], row["outgoing_scc_relations"],
        row["max_containment_asymmetry"]), reverse=True)
    return result


def ambiguous_sink_rows(nodes, components, dag_rows):
    component_count = len(components)
    outgoing = [set() for _ in range(component_count)]
    indegree = [0] * component_count
    direct_rows = defaultdict(list)
    for row in dag_rows:
        source = int(row["source_scc"])
        target = int(row["target_scc"])
        if target not in outgoing[source]:
            outgoing[source].add(target)
            indegree[target] += 1
        direct_rows[source].append(row)

    queue = [component for component in range(component_count)
             if indegree[component] == 0]
    topological = []
    for component in queue:
        topological.append(component)
        for target in outgoing[component]:
            indegree[target] -= 1
            if indegree[target] == 0:
                queue.append(target)
    if len(topological) != component_count:
        raise ValueError("SCC condensation graph unexpectedly contains a cycle")

    reachable_sinks = [set() for _ in range(component_count)]
    for component in reversed(topological):
        if not outgoing[component]:
            reachable_sinks[component].add(component)
        else:
            for target in outgoing[component]:
                reachable_sinks[component].update(reachable_sinks[target])

    representatives = {}
    for component, records in components.items():
        representatives[component] = max(
            records, key=lambda record: (record["length"], -record["node"])
        )

    result = []
    for component, sinks in enumerate(reachable_sinks):
        if len(sinks) < 2:
            continue
        source = representatives[component]
        sink_records = [representatives[sink] for sink in sorted(sinks)]
        sink_taxids = sorted({record["taxid"] for record in sink_records})
        direct = direct_rows[component]
        direct_taxids = sorted({row["target_species_taxid"] for row in direct})
        if set(sink_taxids) == {source["taxid"]}:
            category = "within_taxid_multiple_containers"
        elif len(sink_taxids) == 1:
            category = "foreign_taxid_container_set"
        elif source["taxid"] in sink_taxids:
            category = "mixed_taxid_ambiguity"
        else:
            category = "cross_taxid_ambiguity"
        result.append({
            "scc_id": component,
            "scc_size": len(components[component]),
            "source_accession": source["accession"],
            "source_taxid": source["taxid"],
            "source_name": source["name"],
            "source_length": source["length"],
            "direct_outdegree": len(outgoing[component]),
            "direct_target_taxids": len(direct_taxids),
            "reachable_sinks": len(sinks),
            "reachable_sink_taxids": len(sink_taxids),
            "category": category,
            "sink_taxids": ";".join(sink_taxids),
            "sink_accessions": ";".join(
                record["accession"] for record in sink_records
            ),
            "sink_names": ";".join(record["name"] for record in sink_records),
            "max_direct_containment": max(
                (float(row["containment"]) for row in direct), default=0.0
            ),
            "max_direct_asymmetry": max(
                (float(row["containment_asymmetry"]) for row in direct),
                default=0.0
            ),
            "min_direct_length_ratio": min(
                (float(row["length_ratio"]) for row in direct), default=1.0
            ),
        })
    result.sort(key=lambda row: (
        row["scc_size"] == 1,
        row["reachable_sink_taxids"],
        row["direct_target_taxids"],
        row["reachable_sinks"],
        row["max_direct_asymmetry"],
    ), reverse=True)
    return result


def write_tsv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]),
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--members", required=True)
    parser.add_argument("--ground-truth", required=True)
    parser.add_argument("--dag")
    parser.add_argument("--assignment")
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()

    ground_truth = load_ground_truth(args.ground_truth)
    nodes, components = load_members(args.members, ground_truth)
    names_by_taxid = defaultdict(Counter)
    for record in nodes.values():
        names_by_taxid[record["taxid"]][record["name"]] += 1
    taxid_names = {
        taxid: names.most_common(1)[0][0]
        for taxid, names in names_by_taxid.items()
    }
    splits = species_split_rows(nodes)
    impurities = component_impurity_rows(components, taxid_names)
    dag_rows = annotate_dag(args.dag, nodes)
    assignment_rows = load_assignment(args.assignment, components)
    source_rows = directed_source_rows(dag_rows)
    ambiguity_rows = ambiguous_sink_rows(nodes, components, dag_rows)

    write_tsv(args.output_prefix + ".species_splits.tsv", splits)
    write_tsv(args.output_prefix + ".scc_impurity.tsv", impurities)
    write_tsv(args.output_prefix + ".dag_annotated.tsv", dag_rows)
    write_tsv(args.output_prefix + ".directed_sources.tsv", source_rows)
    write_tsv(args.output_prefix + ".ambiguous_sinks.tsv", ambiguity_rows)

    mapped = sum(record["taxid"] != "NA" for record in nodes.values())
    split_species = sum(row["scc_count"] > 1 for row in splits)
    impure_sccs = sum(row["species_count"] > 1 for row in impurities)
    same_species_dag = sum(row["same_species"] for row in dag_rows)
    short_to_long_dag = sum(
        float(row["length_ratio"]) <= 0.8 for row in dag_rows
    )
    strongly_asymmetric_dag = sum(
        float(row["containment_asymmetry"]) >= 0.1 for row in dag_rows
    )
    unique_sink_components = sum(
        row["status"] == "unique" for row in assignment_rows
    )
    ambiguous_components = sum(
        row["status"] == "ambiguous" for row in assignment_rows
    )
    unique_sink_genomes = sum(
        row["genomes"] for row in assignment_rows
        if row["status"] == "unique"
    )
    ambiguous_genomes = sum(
        row["genomes"] for row in assignment_rows
        if row["status"] == "ambiguous"
    )
    outgoing_components = {
        int(row["source_scc"]) for row in dag_rows
    }
    unique_non_sink = [
        row for row in assignment_rows
        if row["status"] == "unique" and row["scc_id"] in outgoing_components
    ]
    ambiguous_singletons = sum(
        row["scc_size"] == 1 for row in ambiguity_rows
    )
    cross_taxid_ambiguities = sum(
        row["reachable_sink_taxids"] > 1 for row in ambiguity_rows
    )
    direct_multi_sink_ambiguities = sum(
        row["direct_outdegree"] > 1 for row in ambiguity_rows
    )
    summary = [
        {"metric": "nodes", "value": len(nodes)},
        {"metric": "mapped_nodes", "value": mapped},
        {"metric": "species", "value": len(splits)},
        {"metric": "species_split_across_sccs", "value": split_species},
        {"metric": "sccs", "value": len(components)},
        {"metric": "impure_sccs", "value": impure_sccs},
        {"metric": "dag_edges", "value": len(dag_rows)},
        {"metric": "same_species_dag_edges", "value": same_species_dag},
        {"metric": "cross_species_dag_edges",
         "value": len(dag_rows) - same_species_dag},
        {"metric": "short_to_long_dag_edges_length_ratio_le_0.8",
         "value": short_to_long_dag},
        {"metric": "strongly_asymmetric_dag_edges_delta_ge_0.1",
         "value": strongly_asymmetric_dag},
        {"metric": "unique_sink_components", "value": unique_sink_components},
        {"metric": "unique_sink_genomes", "value": unique_sink_genomes},
        {"metric": "ambiguous_components", "value": ambiguous_components},
        {"metric": "ambiguous_genomes", "value": ambiguous_genomes},
        {"metric": "non_sink_components", "value": len(outgoing_components)},
        {"metric": "unique_non_sink_components", "value": len(unique_non_sink)},
        {"metric": "unique_non_sink_genomes", "value": sum(
            row["genomes"] for row in unique_non_sink
        )},
        {"metric": "ambiguous_singleton_components",
         "value": ambiguous_singletons},
        {"metric": "ambiguous_cross_taxid_components",
         "value": cross_taxid_ambiguities},
        {"metric": "ambiguous_direct_multi_sink_components",
         "value": direct_multi_sink_ambiguities},
    ]
    write_tsv(args.output_prefix + ".taxonomy_summary.tsv", summary)

    print("Top species splits (genomes outside largest SCC):")
    for row in splits[:args.top]:
        print("{split_genomes}\t{genomes}\t{scc_count}\t{species_taxid}\t"
              "{organism_name}".format(**row))
    print("\nTop impure SCCs (off-species genomes):")
    for row in impurities[:args.top]:
        print("{off_species_genomes}\t{size}\t{species_count}\t{scc_id}\t"
              "{major_species_taxid}\t{major_organism_name}".format(**row))
    print("\nTop directed relations (containment asymmetry):")
    for row in sorted(dag_rows,
                      key=lambda item: float(item["containment_asymmetry"]),
                      reverse=True)[:args.top]:
        print("{containment_asymmetry:.4f}\t{length_ratio}\t"
              "{source_accession}\t{target_accession}\t"
              "{source_organism_name}\t{target_organism_name}".format(
                  **row))
    print("\nTop multi-sink ambiguities:")
    for row in ambiguity_rows[:args.top]:
        print("{reachable_sink_taxids}\t{reachable_sinks}\t{scc_size}\t"
              "{source_accession}\t{source_taxid}\t{category}\t"
              "{sink_taxids}\t{source_name}".format(**row))


if __name__ == "__main__":
    main()
