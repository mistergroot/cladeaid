from html import parser
import pysam
from collections import defaultdict
import pandas as pd
import numpy as np
import time
from tqdm import tqdm
import gzip
import argparse
import csv
from . import tax_parsing
import os
from . import bam_processor
from . import propagate_counts
from . import mash_matrix
from Bio import SeqIO

def add_arguments(parser):
    parser.add_argument("--bam", required=True, 
                        help="Input BAM file")
    parser.add_argument("--acc2taxid", required=True, 
                        help="Accession2taxid" \
    " mapping file (.tsv or .gz)")
    parser.add_argument("--nodes", required=True, 
                        help="nodes.dmp taxonomy file (can be .gz)")
    parser.add_argument("--names", required=True, 
                        help="names.dmp taxonomy file (can be .gz)")
    parser.add_argument("--output", required=True, 
                        help="Output prefix")
    parser.add_argument("--verbose", action="store_true", 
                        help="Enable verbose output")
    parser.add_argument("--min_identity", type=float, default=0.0, 
                        help="Minimum identity threshold (default: 0.0)")
    parser.add_argument("--max_reads", type=int, 
                        help="Maximum number of reads to process")
    parser.add_argument("--write_discordant_bam", action="store_true", 
                        help="Write BAM of discordant reads")
    parser.add_argument("--discordant_bam_path", default="discordant.bam", 
                        help="Path to write discordant BAM")
    parser.add_argument("--write_specificity_bam", action="store_true", 
                        help="Write BAM of specificity-overridden reads")
    parser.add_argument("--specificity_bam_path", default="specificity.bam", 
                        help="Path to write specificity BAM")
    parser.add_argument("--extract_read_attributes", action="store_true", 
                        help="Extract and output read attributes")

def run(args):
    results = bam_processor.process_bam_streaming(
        bamfile=args.bam,
        acc2taxid_file=args.acc2taxid,
        nodes_file=args.nodes,
        names_file=args.names,
        max_reads=args.max_reads,
        min_identity=args.min_identity,
        verbose=args.verbose,
        write_discordant_bam=args.write_discordant_bam,
        discordant_bam_path=args.discordant_bam_path,
        write_specificity_bam=args.write_specificity_bam,
        specificity_bam_path=args.specificity_bam_path,
        extract_read_attributes=args.extract_read_attributes
    )

    with open(args.output + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        if args.extract_read_attributes:
            writer.writerow([
                "ReadName", "TaxID", "TaxName", "Rank", "TotalBases", 
                "BestIdentity", "ReadPairConcordance", "R1_Length", 
                "R1_Matches", "R1_Insertions", "R1_Deletions", "R1_Softclips",
                "R1_Hardclips", "R1_AvgQual", "R2_Length", "R2_Matches", 
                "R2_Insertions", "R2_Deletions", "R2_Softclips", 
                "R2_Hardclips", "R2_AvgQual"
            ])
        else:
            writer.writerow([
                "ReadName", "TaxID", "TaxName", "Rank", "TotalBases", 
                "BestIdentity", "ReadPairConcordance"
            ])
        writer.writerows(results)

def main():
    parser = argparse.ArgumentParser(description="Assign taxonomic labels " \
    "using LCA + specificity.")
    add_arguments(parser)
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()