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
from . import bam_splitter
from . import mash_matrix
from . import tax_parsing
from Bio import SeqIO

def add_arguments(parser):
    parser.description = (
        "Estimate abundance using naive or mash distance-aware propagation. "
        "Optionally, normalize based on genome size, and output bam file "
        "for each 'confident' taxon."
    )
    parser.add_argument("--csv", required=True, 
                        help="Input CSV from cladeaid.py")
    parser.add_argument("--acc2taxid", required=True, 
                        help="Accession2taxid mapping file (can be .gz)")
    parser.add_argument("--nodes", required=True, 
                        help="nodes.dmp taxonomy file (can be .gz)")
    parser.add_argument("--names", required=True, 
                        help="names.dmp taxonomy file (can be .gz)")
    parser.add_argument("--output", required=True, 
                        help="Output prefix - Assigned reads will be output " \
                        "as <output>.csv, and abundances will be written " \
                        "to <output>.abundances")
    parser.add_argument("--genome_size_scaling", action="store_true", 
                        help="Adjust base abundances based on reference " \
                        "genome lengths")
    parser.add_argument("--bam",
                        help="Input bam file. Only required if "
                        "--output_all_mappings or "
                        "--output_confident_mappings are specified.")
    parser.add_argument("--output_all_mappings", action="store_true", 
                        help="Outputs a bam for each taxon (including high " \
                        "level ones at genus or higher).")
    parser.add_argument("--output_tax_level_mappings", default="none",
                        help="Outputs a bam with reads collapsed to a " \
                        "specified taxonomic level. For example, if 'genus' "
                        "is specified, all reads assigned to that genus, " \
                        "as well as all species-level reads within that " \
                        "genus will be collapsed to that genus and output to" \
                        " <output>_<genus>.bam. Accepted levels are " \
                        "'kingdom', 'phylum', 'class', 'order', 'family', " \
                        "'genus', 'species'.")
    parser.add_argument("--output_leaf_mappings", action="store_true", 
                        help="Outputs a bam for each leaf (the most specific" \
                        " taxon in each lineage) that is considered high " \
                        "confidence (i.e. >0.05 percent abundance. If "
                        "--mash_reallocation was used, confident also " \
                        "refers to species that have no other close " \
                        "relatives with <0.05 mash distance)")
    parser.add_argument("--mash_reallocation", action="store_true", 
                        help="Use mash to calculate distances and reallocate" \
                        " bases to more specific taxa")
    parser.add_argument("--reference_genome_list", default="references.list", 
                        help="Path to list of reference genomes for mash " \
                        "reallocation and genome size adjustment")
    parser.add_argument("--pairwise_dists", default="all.dists",
                        help="Path to pairwise distance file between genomes" \
                        " for mash reallocation")
    parser.add_argument("--multifasta", default=False, action="store_true", 
                        help="Flag if BAMs were mapped against a " \
                        "multi-organism fasta reference")
    parser.add_argument("--normalize", action="store_true",
                        help="Normalize bases assigned to taxa based on " \
                        "genome sizes")
    parser.add_argument("--threads", type=int, default=4, 
                        help="Number of threads for mash distance matrix " \
                        "estimation (default: 4)")

def run(args):
    assignments = pd.read_csv(args.csv)
    assignments['Assigned_Bases'] = (
        assignments.groupby(['TaxID', 'TaxName'])
        ['TotalBases'].transform('sum')
        )
    assignments = (
        assignments.drop_duplicates(subset=['TaxID', 'TaxName'])
        .reset_index(drop=True)
        )
    
    observed_read_counts = dict(zip(assignments["TaxID"], 
                                    assignments["Assigned_Bases"]))
    taxid_list = assignments["TaxID"].tolist()

    if args.mash_reallocation:
        print("Beginning distance matrix calculation and read propagation...")
        if os.path.exists(args.pairwise_dists):
            distmatrix, taxrename = mash_matrix.make_dist_matrix(
                args.reference_genome_list, args.pairwise_dists, 
                args.acc2taxid, multifasta=args.multifasta, 
                threads=args.threads)
        else:
            distmatrix, taxrename = mash_matrix.make_dist_matrix(
                args.reference_genome_list, (args.output + ".dists"), 
                args.acc2taxid, multifasta=args.multifasta, 
                threads=args.threads)
            
        distmatrix = distmatrix.rename(index=taxrename, columns=taxrename)

        grouped_rows = distmatrix.groupby(distmatrix.index).mean()
        grouped_full = grouped_rows.T.groupby(distmatrix.columns).mean().T

        dedup_distmatrix = (grouped_full + grouped_full.T) / 2
        np.fill_diagonal(dedup_distmatrix.values, 0)
        
        naive_abundances, penalized_abundances = (
            propagate_counts.propagate_counts(
            taxid_list=taxid_list, 
            nodes_path=args.nodes, 
            names_path=args.names,
            observed_read_counts=observed_read_counts,
            mash_penalty=True,
            distance_matrix=dedup_distmatrix
            )
        )
        print("Naive abundances and penalized abundances calculated.")

    else:
        print("Naively propagating counts down the tree")
        naive_abundances, penalized_abundances = (
            propagate_counts.propagate_counts(
            taxid_list=taxid_list, 
            nodes_path=args.nodes, 
            names_path=args.names,
            observed_read_counts=observed_read_counts
            )
        )
        print("Naive abundances calculated.")
    
    naive_abundances = pd.DataFrame(naive_abundances.items(), 
                                    columns=["TaxID","Naive_Bases"])
    penalized_abundances = pd.DataFrame(penalized_abundances.items(), 
                                        columns=["TaxID", "Penalized_Bases"])
    
    propagated_counts = (
        assignments[["TaxID", "TaxName", "Rank", "Assigned_Bases"]]
        .merge(naive_abundances, how="left", on="TaxID")
        .merge(penalized_abundances, how="left", on="TaxID"))
    
    propagated_counts["Assigned_Proportion"] = (
        propagated_counts["Assigned_Bases"] / 
        propagated_counts["Assigned_Bases"].sum()
    )
    propagated_counts["Naive_Proportion"] = (
        propagated_counts["Naive_Bases"] / 
        propagated_counts["Naive_Bases"].sum()
    )
    propagated_counts["Penalized_Proportion"] = (
        propagated_counts["Penalized_Bases"] / 
        propagated_counts["Penalized_Bases"].sum()
    )

    if args.normalize: 
        print("Normalizing abundances based on genome sizes.")
        reference_sizes = {}
        refs = open(args.reference_genome_list).read().splitlines()
        for ref in refs:
            reference_sizes[ref] = mash_matrix.count_ref_size(ref)
        reference_sizes = pd.DataFrame(reference_sizes.items())
        reference_sizes[0] = reference_sizes[0].replace(taxrename)
        reference_sizes = (reference_sizes.groupby(0)[1].mean().reset_index())
        reference_sizes = dict(zip(reference_sizes[0], reference_sizes[1]))

        for prop in ["Naive", "Penalized"]:
            propagated_counts[prop + "_Proportion"] = np.nan
            for i in range(len(propagated_counts)):
                if not np.isnan(propagated_counts[prop + "_Bases"][i]):
                    propagated_counts.at[i, prop + "_Proportion"] = (
                        propagated_counts[prop + "_Bases"][i] / 
                        reference_sizes.get(propagated_counts["TaxID"][i], 
                                            np.nan)
                    )
            propagated_counts[prop + "_Proportion"] = (
                propagated_counts[prop + "_Proportion"] / 
                propagated_counts[prop + "_Proportion"].sum()
            )
        print("Normalized abundances calculated based on genome sizes.")
    
    propagated_counts.to_csv(args.output + ".abundances", index=False)
    
    if args.output_all_mappings:
        taxa = propagated_counts.TaxName.unique().tolist()
        bam_splitter.split_bam_by_taxid(args.bam, 
                                        args.output, 
                                        args.csv,
                                        taxa)
        print("BAM files split by taxon.")

    if args.output_leaf_mappings:
        taxa = (propagated_counts.TaxName
                [propagated_counts.Naive_Proportion.notna()].unique().tolist())
        if args.mash_reallocation:
            taxid = (propagated_counts.TaxID
                [propagated_counts.Naive_Proportion.notna()]
                .unique().tolist())
            subset_matrix = dedup_distmatrix.loc[taxid, taxid]
            mask = (subset_matrix >= 0.05) | (subset_matrix == 0)
            valid_species = mask.all(axis=1)
            conf_matrix = subset_matrix.loc[valid_species.index[valid_species], 
                                            valid_species.index[valid_species]]
            taxa = (propagated_counts.TaxName
                [propagated_counts.Naive_Proportion.notna()]
                [propagated_counts.TaxID.isin(conf_matrix.columns.tolist())]
                .unique().tolist())
        bam_splitter.split_bam_by_taxid(args.bam, 
                                        args.output, 
                                        args.csv,
                                        taxa)
        print("BAM files for all confident leaves")
        
    if args.output_tax_level_mappings != "none":
        taxid_to_name, taxid_to_rank, parent_map = (
            tax_parsing.load_taxonomy(args.nodes, args.names)
        )
        taxa = propagated_counts.TaxID.unique().tolist()
        user_specified_level_map = {}
        for tax in taxa:
            user_level = tax_parsing.get_ancestor_at_level(tax, 
                                               args.output_tax_level_mappings, 
                                               parent_map, taxid_to_rank, 
                                               taxid_to_name=taxid_to_name)
            if user_level is not None:
                user_specified_level_map[tax] = user_level
        bam_splitter.split_bam_by_taxid(args.bam, 
                                        args.output, 
                                        args.csv,
                                        user_specified_level_map)
        print(f"BAM files output at the {args.output_tax_level_mappings} level.")

def main():
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()