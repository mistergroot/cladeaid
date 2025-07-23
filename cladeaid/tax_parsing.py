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