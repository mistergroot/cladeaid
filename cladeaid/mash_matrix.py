import os
import subprocess
import pandas as pd
import numpy as np
from io import StringIO
from Bio import SeqIO
from . import tax_parsing
import os

def run_mash(ref_list, outfile, 
             multifasta=False, threads=4):
    print("Running Mash sketch and distance calculations...")
    if threads < 8:
        jobs = "1"
        jobthreads = str(int(threads))
    else:
        jobs = str(np.floor(threads / 8).astype(int))
        jobthreads = "8"

    ps = subprocess.run(["cat", ref_list], check=True, stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE)
    if multifasta==True:
        processNames = subprocess.run(['parallel', '-j', jobs, 'mash', 
                                       'sketch', '-i', '-p', jobthreads, '{}', 
                                       '-s', '10000'], input=ps.stdout, 
                                       stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE)
    else:
        processNames = subprocess.run(['parallel', '-j', jobs, 'mash', 'sketch', 
                                       '-p', jobthreads, '{}', '-s', '10000'],
                                       input=ps.stdout, stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE)

    ps = subprocess.run(["cat", ref_list], check=True, stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE)
    if multifasta==True:
        processNames = subprocess.run(['parallel', '-j', jobs, 'mash', 'dist', 
                                   '-i', '-p', jobthreads, '{}', 
                                   ('$( cat ' + ref_list + ')'), '>>', 
                                   outfile], 
                                   input=ps.stdout, stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
    else:
        processNames = subprocess.run(['parallel', '-j', jobs, 'mash', 'dist', 
                                   '-p', jobthreads, '{}', 
                                   ('$( cat ' + ref_list + ')'), '>>', 
                                   outfile], 
                                   input=ps.stdout, stdout=subprocess.PIPE, 
                                   stderr=subprocess.PIPE)
    print("Mash sketch and distance calculations completed.")
    
def count_ref_size(fasta_file):
    total_bases = 0
    try:
        with tax_parsing.smart_open(fasta_file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                total_bases += len(record.seq)
    except FileNotFoundError:
            print(f"Error: File not found at {fasta_file}")
    except Exception as e:
            print(f"An error occurred while processing {fasta_file}: {e}")
    
    return total_bases

def make_dist_matrix(ref_list, outfile, acc2tid, 
                     multifasta=False, threads=4):
    if os.path.exists(outfile):
        print("Mash distance matrix already exists. Skipping computation.")
    else:
        run_mash(ref_list, outfile, 
                 multifasta=multifasta, threads=threads)
    df = pd.read_csv(outfile, sep="\t", header=None,
                     names=["species1", "species2", 
                            "distance", "shared", "pvalue"])

    # Get unique species names
    species = sorted(set(df["species1"]).union(set(df["species2"])))

    taxid_map_dist = {}
    if multifasta:
        with open(outfile + ".tmp", 'w') as file:
            for item in species:
                file.write(f"{item}\n")

        ps = subprocess.run(["zgrep", "-f", (outfile + ".tmp"), acc2tid], 
                                check=True, stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE, universal_newlines=True)
        processNames = subprocess.run(['cut', '-f2,3'],
                                    input=ps.stdout, 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE, 
                                    universal_newlines=True)
        taxid_map_dist = {
            key: int(value)
            for line in processNames.stdout.strip().split('\n') if line
            for key, value in [line.split('\t')]
        }

    else:
        sp_to_header = {}
        for sp in species:
            with tax_parsing.smart_open(sp) as handle:
                first_record = next(SeqIO.parse(handle, "fasta"))
                sp_to_header[sp] = first_record.id

        with open(outfile + ".tmp", 'w') as file:
            for header in sp_to_header.values():
                file.write(f"{header}\n")

        ps = subprocess.run(["zgrep", "-f", (outfile + ".tmp"), acc2tid],
                            check=True, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)

        processNames = subprocess.run(['cut', '-f2,3'],
                                    input=ps.stdout,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    universal_newlines=True)

        # Map header → taxid
        header_to_taxid = {
            key: int(value)
            for line in processNames.stdout.strip().split('\n') if line
            for key, value in [line.split('\t')]
        }

        # Now map sp → taxid using the sp_to_header mapping
        taxid_map_dist = {
            sp: header_to_taxid[header]
            for sp, header in sp_to_header.items()
            if header in header_to_taxid
        }

    # Initialize distance matrix
    dist_matrix = pd.DataFrame(0.0, index=species, columns=species)

    # Fill the matrix with Mash distances
    for _, row in df.iterrows():
        s1 = row["species1"]
        s2 = row["species2"]
        dist_matrix.loc[s1, s2] = row["distance"]
        dist_matrix.loc[s2, s1] = row["distance"]

    return dist_matrix, taxid_map_dist