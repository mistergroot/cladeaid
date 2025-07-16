from collections import defaultdict
import gzip
import tax_parsing
from tax_parsing import parse_nodes_dmp
from tax_parsing import parse_names_dmp
from tax_parsing import smart_open

def get_first_level_children(taxid_list, nodes_path, names_path):
    taxid_set = set(taxid_list)
    children_dict = defaultdict(list)

    parent_map, rank_map = parse_nodes_dmp(nodes_path)

    # Track depth for sorting
    taxid_depth = {}

    for candidate in taxid_list:
        current = candidate
        while True:
            parent = parent_map.get(current)
            if parent is None or parent == current:
                break
            if parent in taxid_set:
                children_dict[parent].append(candidate)
                break
            current = parent

    # Remove terminal taxa (i.e., keys with no children)
    children_dict = {k: v for k, v in children_dict.items() if v}

    return children_dict

def propagate_counts(
    taxid_list, nodes_path, names_path,
    observed_read_counts,
    pseudocount=1.0,
    max_iter=100,
    tol=1e-6,
    mash_penalty=False):
    """
    Propagating counts from higher taxonomic levels down to child taxa

    Parameters:
    - parent_to_children: dict of {parent_taxid: [child_taxids]} from get_first_level_children()
    - observed_read_counts: dict of {taxid: observed counts} at leaves or any level
    - pseudocount: value to avoid zero weights
    - max_iter: max number of EM iterations
    - tol: convergence tolerance
    """
    #sorted_dict, named_dict = get_first_level_children(tmp["TaxID"].tolist(), "../data/nodes.dmp.gz", "../data/names.dmp.gz")
    name_map = parse_names_dmp(names_path)
    
    parent_to_children = get_first_level_children(taxid_list, 
                                                  nodes_path, 
                                                  names_path)
    
    # Step 1: Build complete taxid set
    taxid_set = set(observed_read_counts)
    for parent, children in parent_to_children.items():
        taxid_set.add(parent)
        taxid_set.update(children)

    # Step 2: Build reverse mapping from child → parent
    child_to_parent = {}
    for parent, children in parent_to_children.items():
        for child in children:
            child_to_parent[child] = parent

    # Step 3: Initialize current counts
    current_counts = defaultdict(float)
    for taxid in taxid_set:
        current_counts[taxid] = float(observed_read_counts.get(taxid, 0.0))

    # Step 4: Identify ambiguous parents that need redistribution
    ambiguous_nodes = [
        taxid for taxid in parent_to_children
        if taxid in observed_read_counts
    ]

    # Step 5: EM iterations
    for _ in range(max_iter):
        new_counts = defaultdict(float)

        # Copy over fixed assignments for non-ambiguous taxa
        for taxid in taxid_set:
            if taxid not in ambiguous_nodes:
                new_counts[taxid] += current_counts[taxid]

        max_delta = 0.0

        for parent in ambiguous_nodes:
            children = parent_to_children.get(parent, [])
            weights = [current_counts[c] + pseudocount for c in children]
            total_weight = sum(weights)
            if total_weight == 0:
                continue
            parent_count = current_counts[parent]
            for child, weight in zip(children, weights):
                redistributed = parent_count * (weight / total_weight)
                max_delta = max(max_delta, abs(redistributed - new_counts[child]))
                new_counts[child] += redistributed

        current_counts = new_counts
        if max_delta < tol:
            break

    propagated_counts = {k: v for k, v in dict(current_counts).items() if v != 0}
    
    named_dict = {}
    for parent, counts in propagated_counts.items():
        parent_name = name_map.get(parent, f"TAXID_{parent}")
        named_dict[parent_name] = counts

    return propagated_counts, parent_to_children, named_dict
