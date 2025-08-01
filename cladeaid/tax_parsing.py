import gzip

def smart_open(filepath):
   return gzip.open(filepath, 
                    'rt') if filepath.endswith('.gz') else open(filepath, 'r')

# --- Step 1: Parse taxonomy files ---
def parse_nodes_dmp(path):
    parent = {}
    rank = {}
    with smart_open(path) as f:
        for line in f:
            fields = [x.strip() for x in line.split("|")]
            taxid = int(fields[0])
            parent_id = int(fields[1])
            tax_rank = fields[2]
            parent[taxid] = parent_id
            rank[taxid] = tax_rank
    return parent, rank

def parse_names_dmp(path):
    name_map = {}
    with smart_open(path) as f:
        for line in f:
            fields = [x.strip() for x in line.split("|")]
            taxid = int(fields[0])
            name_txt = fields[1]
            class_txt = fields[3]
            if class_txt == "scientific name":
                name_map[taxid] = name_txt
    return name_map

def load_taxonomy(nodes_fp, names_fp):
    parent_map, taxid_to_rank = parse_nodes_dmp(nodes_fp)
    taxid_to_name = parse_names_dmp(names_fp)
    return taxid_to_name, taxid_to_rank, parent_map

def get_ancestor_at_level(taxid, target_level, parent_map, 
                          taxid_to_rank, taxid_to_name=None):
    current = taxid
    while current in parent_map:
        rank = taxid_to_rank.get(current, "")
        if rank == target_level:
            if taxid_to_name:
                return taxid_to_name.get(current, f"taxid:{current}")
            else:
                return current
        if parent_map[current] == current:
            # reached root
            break
        current = parent_map[current]
    return None