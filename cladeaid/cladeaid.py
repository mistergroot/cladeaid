from html import parser
import pysam
from collections import defaultdict
import pandas as pd
import time
from tqdm import tqdm
import gzip
import argparse
import csv
import tax_parsing
#from tax_parsing import parse_nodes_dmp
#from tax_parsing import parse_names_dmp
#from tax_parsing import smart_open
import bam_processor
import propagate_counts
import mash_matrix

def main():
    parser = argparse.ArgumentParser(description="Assign taxonomic labels to" \
    " BAM reads using LCA + specificity-based logic.")

    parser.add_argument("--bam", required=True, 
                        help="Input BAM file")
    parser.add_argument("--acc2taxid", required=True, 
                        help="Accession2taxid mapping file (.tsv or .gz)")
    parser.add_argument("--nodes", required=True, 
                        help="nodes.dmp taxonomy file (can be .gz)")
    parser.add_argument("--names", required=True, 
                        help="names.dmp taxonomy file (can be .gz)")
    parser.add_argument("--output", required=True, 
                        help="Output prefix - Assigned reads will be output as"
                        " <output>.csv, and, if --estimate_abundance is used" \
                        ", abundances will be written to <output>.abundances")
    parser.add_argument("--verbose", action="store_true", 
                        help="Enable verbose output")
    parser.add_argument("--min_identity", type=float, default=0.0, 
                        help="Optional - Minimum identity threshold " \
                        "(default: 0.0)")
    parser.add_argument("--max_reads", type=int, 
                        help="Optional - maximum number of reads to process")
    parser.add_argument("--write_discordant_bam", action="store_true", 
                        help="Optional - Write BAM of discordant reads")
    parser.add_argument("--discordant_bam_path", default="discordant.bam", 
                        help="Optional - Path to write discordant BAM")
    parser.add_argument("--write_specificity_bam", action="store_true", 
                        help="Optional - Write BAM of specificity-overridden" \
                        " reads")
    parser.add_argument("--specificity_bam_path", default="specificity.bam", 
                        help="Optional - Path to write specificity BAM")
    parser.add_argument("--extract_read_attributes", action="store_true", 
                        help="Optional - Extract and output read length, " \
    "CIGAR, and quality metrics")
    parser.add_argument("--estimate_abundance", action="store_true", 
                        help="Optional - Estimate abundance by " \
                        "proportionally propagating bases from ambiguous " \
                        "taxonomic levels to more specific taxonomic levels." \
                        " This will output a file <output>.abundances with " \
                        "taxon names and assigned bases. By default, this " \
                        "done naively. If there are 3 species in a genus, " \
                        "the reads will be allocated proportional to their" \
                        " species-level base abundances, regardless of their" \
                        " phylogenetic distance. For " \
                        "phylogenetically-informed scaling, use "
                        "--mash_reallocation")
    parser.add_argument("--genome_size_scaling", action="store_true", 
                        help="Optional - Adjust base abundances based on " \
                        "reference genome lengths")
    parser.add_argument("--mash_reallocation", action="store_true", 
                        help="Optional - Use mash to calculate distances and" \
                        " reallocate bases to more specific taxa, taking " \
                        "into account that more distant taxa are less " \
                        "likely to be misassigned")
    parser.add_argument("--reference_genome_list", default="references.list", 
                        help="Optional - Path to a list of reference " \
                        "genomes to use for mash reallocation and adjusting " \
                        "base abundances based on reference genome lengths. " \
                        "Required if --genome_size_scaling and/or "
                        "--mash_reallocation is used")
    parser.add_argument("--threads", type=int, default=4, 
                        help="Optional - Number of threads to use for mash " \
    "distance matrix estimation (default: 4)")

    args = parser.parse_args()

    results = bam_processor.process_bam_streaming(
        bamfile=args.bam,
        acc2taxid_file=args.acc2taxid,
        nodes_file=args.nodes,
        names_file=args.names,
        max_reads=args.max_reads,
        min_identity=args.min_identity,
        verbose=True,
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

    if args.estimate_abundance:
        assignments = pd.read_csv(args.output + ".csv")
        assignments['Bases'] = (assignments.groupby(['TaxID', 'TaxName'])
                                ['TotalBases'].transform('sum'))
        assignments = (assignments.drop_duplicates(subset=['TaxID', 'TaxName'])
                       .reset_index(drop=True))
        observed_read_counts=dict(zip(assignments["TaxID"], 
                                      assignments["Bases"]))
        taxid_list = assignments["TaxID"].tolist()
        if args.mash_reallocation:
            distmatrix, taxrenamemap_dists = mash_matrix.make_dist_matrix(
                args.reference_genome_list, args.output, 
                args.threads, args.acc2taxid)
            naive_abundances, penalized_abundances = (propagate_counts.
                                                      propagate_counts(
                taxid_list=taxid_list, 
                nodes_path=args.nodes, 
                names_path=args.names,
                observed_read_counts=observed_read_counts,
                mash_penalty=True,
                distance_matrix=distmatrix.rename(index=taxrenamemap_dists, 
                                                  columns=taxrenamemap_dists)
                ))
        else:
            naive_abundances, penalized_abundances = (propagate_counts.
                                                      propagate_counts(
                taxid_list=taxid_list, 
                nodes_path=args.nodes, 
                names_path=args.names,
                observed_read_counts=observed_read_counts
                ))
        naive_abundances = pd.DataFrame(naive_abundances.items())
        naive_abundances[2] = naive_abundances[1] / naive_abundances[1].sum()
        naive_abundances.columns = ["TaxID", "Naive_Bases", "Naive_Proportion"]
        penalized_abundances = pd.DataFrame(penalized_abundances.items())
        penalized_abundances[2] = (penalized_abundances[1] / 
                                   penalized_abundances[1].sum())
        penalized_abundances.columns = ["TaxID", "Penalized_Bases", 
                                        "Penalized_Proportion"]
        assignments["Proportion"] = (assignments["Bases"] / 
                                     assignments["Bases"].sum())
        propagated_counts = (assignments[["TaxID", "TaxName", "Rank", "Bases", 
                                          "Proportion"]]
                             .merge(naive_abundances, how = "left", 
                                    on = "TaxID")
                             .merge(penalized_abundances, how = "left", 
                                    on = "TaxID"))
        propagated_counts.to_csv(args.output + ".abundances", index=False)

if __name__ == "__main__":
    main()