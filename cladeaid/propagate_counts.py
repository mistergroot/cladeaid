from collections import defaultdict
import gzip
import numpy as np
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

def get_terminal_tip_counts(parent_to_children):
    from collections import defaultdict

    tip_counts = {}

    def count_tips(taxid):
        if taxid not in parent_to_children:
            return 1  # It's a terminal tip
        if taxid in tip_counts:
            return tip_counts[taxid]
        total = sum(count_tips(child) for child in parent_to_children[taxid])
        tip_counts[taxid] = total
        return total

    for parent in parent_to_children:
        count_tips(parent)

    # Sort taxids from most to fewest terminal descendants
    sorted_taxa = sorted(tip_counts.items(), key=lambda x: -x[1])
    return [taxid for taxid, _ in sorted_taxa]

def sum_descendant_counts(taxid, parent_to_children, current_counts):
    """Sum this taxid's count and all counts from its terminal descendants."""
    total = current_counts.get(taxid, 0.0)
    if taxid in parent_to_children:
        for child in parent_to_children[taxid]:
            total += sum_descendant_counts(child, parent_to_children, 
                                           current_counts)
    return total

def get_all_leaves(taxid, parent_to_children):
    """Recursively collect all terminal descendant taxids under the 
    given taxid."""
    if taxid not in parent_to_children:
        return [taxid]  # It's a terminal tip
    else:
        leaves = []
        for child in parent_to_children[taxid]:
            leaves.extend(get_all_leaves(child, parent_to_children))
        return leaves

def min_pairwise_mash(all_tips, distance_matrix, epsilon=1e-8):
    """Compute minimum Mash distance between every pair of children in a 
    dict of all_tips."""
    penalties = {}
    
    for focal, focal_tips in all_tips.items():
        min_dist = float("inf")

        for other, other_tips in all_tips.items():
            if other == focal:
                continue

            for ft in focal_tips:
                for ot in other_tips:
                    if (ft in distance_matrix.index 
                        and ot in distance_matrix.columns):
                        dist = distance_matrix.at[ft, ot]
                        if dist < min_dist:
                            min_dist = dist

        # Fall back to maximum distance in dist_matrix if no distances 
        # found for the specific lineages (this shouldn't happen, but just 
        # in case)
        if min_dist == float("inf"):
            min_dist = np.max(distance_matrix)

        # Invert and penalize
        penalties[focal] = (1.0 / (min_dist + epsilon))

    return penalties

def propagate_counts(
    taxid_list, 
    nodes_path, 
    names_path,
    observed_read_counts,
    pseudocount=1.0,
    max_iter=100,
    tol=1e-6,
    mash_penalty=False,
    distance_matrix=None,
    distance_power=1.0
):
    """
    Propagate ambiguous genus-level counts to species using both observed 
    read counts and optional Mash distance penalties.

    Parameters:
    - taxid_list: list of all taxids involved
    - nodes_path, names_path: taxonomy file paths
    - observed_read_counts: dict {taxid: read count}
    - pseudocount: value added to each child to avoid division by zero
    - mash_penalty: bool, use Mash distance to penalize divergent children. 
    Two dicts are returned, 
    one with naive EM propagation, and one with EM propagation and mash 
    penalty. If mash_penalty=False, two dicts are still returned, but they 
    will be identical.
    - distance_matrix: pd.DataFrame, Mash distance matrix (index/cols = 
    species names or taxids)
    - distance_epsilon: small constant to avoid division by zero
    - distance_power: exponent to control effect of Mash distance 
    (higher = more penalizing)
    
    Returns:
    - naive_propagated_counts: counts using EM based only on child ratios
    - mash_propagated_counts: counts penalized using Mash distance (identical 
    to naive if mash_penalty=False)
    """
    from collections import defaultdict
    import numpy as np

    name_map = parse_names_dmp(names_path)
    parent_to_children = get_first_level_children(taxid_list, nodes_path, 
                                                  names_path)
    ordered_taxa = get_terminal_tip_counts(parent_to_children)
    ambiguous_nodes = [t for t in ordered_taxa if t in observed_read_counts]

    # Step 1: Build full taxid set
    taxid_set = set(observed_read_counts)
    for parent, children in parent_to_children.items():
        taxid_set.add(parent)
        taxid_set.update(children)

    # Step 2: Build reverse mapping
    child_to_parent = {}
    for parent, children in parent_to_children.items():
        for child in children:
            child_to_parent[child] = parent

    # Step 3: Initialize counts
    initial_counts = defaultdict(float)
    for taxid in taxid_set:
        initial_counts[taxid] = float(observed_read_counts.get(taxid, 0.0))

    def run_em(counts, use_mash=False):
        current_counts = counts.copy()
        for _ in range(max_iter):
            new_counts = defaultdict(float)
            max_delta = 0.0

            # Fixed counts
            for taxid in taxid_set:
                if taxid not in parent_to_children:
                    new_counts[taxid] += current_counts[taxid]

            # Redistribute ambiguous counts
            for parent in ambiguous_nodes:
                if parent not in current_counts or current_counts[parent] == 0:
                    continue
                children = parent_to_children.get(parent, [])
                parent_count = current_counts[parent]

                # Calculcating weights for naive propagation (counts of all 
                # children)
                weights = {}
                for child in children:
                    count_part = sum_descendant_counts(child, 
                                                       parent_to_children, 
                                                       observed_read_counts)
                    weights[child] = count_part + pseudocount

                # Finding minimum mash distance between each child (in the 
                # case of genera or higher, we need to recursively find all 
                # leaves, and then compare all leaves in one lineage against 
                # all the leaves in another
                if use_mash and distance_matrix is not None:
                    all_tips = {}
                    for child in children:
                        all_tips[child] = get_all_leaves(child, 
                                                         parent_to_children)
                    penalties = min_pairwise_mash(all_tips, distance_matrix, 
                                                  epsilon=1e-8)
                    for child in weights.keys():
                        weights[child] = weights[child] * (penalties[child] 
                                                           ** distance_power)

                total_weight = sum(weights.values())
                if total_weight == 0:
                    continue

                for child in weights.keys():
                    redistributed = parent_count * (weights[child] / 
                                                    total_weight)
                    max_delta = max(max_delta, 
                                    abs(redistributed - new_counts[child]))
                    new_counts[child] += redistributed

            current_counts = new_counts
            if max_delta < tol:
                break
        return current_counts

    # Run naive EM
    naive_counts = run_em(initial_counts, use_mash=False)

    # Run Mash-penalized EM (or return naive twice if mash_penalty=False)
    if mash_penalty:
        mash_counts = run_em(initial_counts, use_mash=True)
    else:
        mash_counts = naive_counts.copy()

    # Output only nonzero results
    naive_final = {k: v for k, v in dict(naive_counts).items() if v != 0}
    mash_final = {k: v for k, v in dict(mash_counts).items() if v != 0}

    return naive_final, mash_final