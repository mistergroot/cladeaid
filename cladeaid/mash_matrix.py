import os
import subprocess
import pandas as pd
from io import StringIO
from Bio import SeqIO
import tax_parsing
import os

def run_mash(ref_list, outprefix = ".", threads = 4):
    if threads < 8:
        jobs = "1"
        jobthreads = str(int(threads))
    else:
        jobs = str(np.floor(threads / 8).astype(int))
        jobthreads = "8"

    ps = subprocess.run(["cat", ref_list], check=True, capture_output=True)
    processNames = subprocess.run(['parallel', '-j', jobs, 'mash', 'sketch', '-p', jobthreads, 
                                   '{}', '-s', '10000'],
                                  input=ps.stdout, capture_output=True)

    ps = subprocess.run(["cat", ref_list], check=True, capture_output=True)
    processNames = subprocess.run(['parallel', '-j', jobs, 'mash', 'dist', '-p', jobthreads, 
                                   '{}', ('$( cat ' + ref_list + ')'), '>>', (outprefix + '.dists')],
                                  input=ps.stdout, capture_output=True)

def make_dist_matrix(ref_list, outprefix, threads, acc2tid):
    if os.path.exists((outprefix + '.dists')):
        print("Mash distance matrix already exists. Skipping computation.")
    else:
        run_mash(ref_list, outprefix, threads)
    df = pd.read_csv((outprefix + '.dists'), sep="\t", header=None,
                     names=["species1", "species2", "distance", "shared", "pvalue"])

    # Get unique species names
    species = sorted(set(df["species1"]).union(set(df["species2"])))

    taxid_map_dist = {}
    for sp in species:
        fasta_file = sp
        with tax_parsing.smart_open(fasta_file) as handle:
            first_record = next(SeqIO.parse(handle, "fasta"))
            ps = subprocess.run(["zgrep", first_record.id, acc2tid], 
                        check=True, capture_output=True, text = True)
            processNames = subprocess.run(['cut', '-f3'],
                                          input=ps.stdout, capture_output=True, text = True)
            taxid_map_dist[sp] = int(processNames.stdout.replace("\n", ""))

    # Initialize distance matrix
    dist_matrix = pd.DataFrame(0.0, index=species, columns=species)

    # Fill the matrix with Mash distances
    for _, row in df.iterrows():
        s1 = row["species1"]
        s2 = row["species2"]
        dist_matrix.loc[s1, s2] = row["distance"]
        dist_matrix.loc[s2, s1] = row["distance"]

    return dist_matrix, taxid_map_dist
