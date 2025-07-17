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
from tax_parsing import parse_nodes_dmp
from tax_parsing import parse_names_dmp
from tax_parsing import smart_open
import propagate_counts
import mash_matrix

def parse_accession2taxid(acc2taxid_file, bamfile):
    print("🔎 Scanning BAM for reference names...")
    samfile = pysam.AlignmentFile(bamfile, "rb")
    bam_refs = set(samfile.references)
    samfile.close()
    print(f"✅ Found {len(bam_refs):,} reference names in BAM")

    ref_to_taxid = {}
    with smart_open(acc2taxid_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            acc, taxid = parts[1], parts[2]
            if acc in bam_refs:
                ref_to_taxid[acc] = int(taxid)

    refid_to_taxid = {}
    samfile = pysam.AlignmentFile(bamfile, "rb")
    for i, name in enumerate(samfile.references):
        if name in ref_to_taxid:
            refid_to_taxid[i] = ref_to_taxid[name]
    samfile.close()

    print(f"✅ Mapped {len(refid_to_taxid)} reference IDs to taxids")
    return refid_to_taxid

# --- Step 2: Taxonomy logic ---
def lineage(taxid, parent):
    path = []
    while taxid != 1 and taxid in parent:
        path.append(taxid)
        taxid = parent[taxid]
    path.append(1)
    return path

def find_lca(taxid1, taxid2, parent):
    lineage1 = set(lineage(taxid1, parent))
    for anc in lineage(taxid2, parent):
        if anc in lineage1:
            return anc
    return 1

# --- Step 3: Streaming BAM processor ---
def process_bam_streaming(bamfile, acc2taxid_file, nodes_file, names_file,
                          max_reads=None, min_identity=0.0, verbose=True,
                          write_discordant_bam=False, discordant_bam_path="discordant.bam",
                          write_specificity_bam=False, specificity_bam_path="specificity.bam",
                          extract_read_attributes=False):
    start_time = time.time()

    print("🔄 Loading taxonomy data...")
    parent, rank = parse_nodes_dmp(nodes_file)
    names = parse_names_dmp(names_file)
    refid_to_taxid = parse_accession2taxid(acc2taxid_file, bamfile)

    print(f"📖 Opening BAM file: {bamfile}")
    samfile = pysam.AlignmentFile(bamfile, "rb")

    discordant_bam = None
    specificity_bam = None

    if write_discordant_bam:
        discordant_bam = pysam.AlignmentFile(discordant_bam_path, "wb", header=samfile.header)
    if write_specificity_bam:
        specificity_bam = pysam.AlignmentFile(specificity_bam_path, "wb", header=samfile.header)

    output = []
    prev_qname = None
    pair_reads = []

    last_print = time.time()

    processed, paired_count, unpaired_count = 0, 0, 0
    same_taxid_count, discordant_taxid_count, less_specific_taxid_count = 0, 0, 0

    for read in samfile.fetch(until_eof=True):
        if prev_qname is None:
            prev_qname = read.query_name

        if read.query_name != prev_qname:
            if len(pair_reads) > 0:
                has_r1 = any(r.is_read1 for r in pair_reads)
                has_r2 = any(r.is_read2 for r in pair_reads)
                if has_r1 and has_r2:
                    paired_count += 1
                else:
                    unpaired_count += 1

                result = assign_lca_from_pair(pair_reads, refid_to_taxid, parent, rank, names, min_identity, extract_read_attributes)
                if result:
                    output.append(result)
                    if result[6] == "same":
                        same_taxid_count += 1
                    elif result[6] == "discordant":
                        discordant_taxid_count += 1
                        if write_discordant_bam:
                            for r in pair_reads:
                                discordant_bam.write(r)
                    elif result[6] in ("R1 less specific", "R2 less specific"):
                        less_specific_taxid_count +=1
                        if write_specificity_bam:
                            for r in pair_reads:
                                specificity_bam.write(r)
            pair_reads = []

        pair_reads.append(read)
        prev_qname = read.query_name

        processed += 1
        if max_reads and len(output) >= max_reads:
            break
        if time.time() - last_print > 2:
            print(f"⏳ Processed {processed:,} mappings → {len(output):,} reads assigned")
            last_print = time.time()

    if len(pair_reads) > 0:
        has_r1 = any(r.is_read1 for r in pair_reads)
        has_r2 = any(r.is_read2 for r in pair_reads)
        if has_r1 and has_r2:
            paired_count += 1
        else:
            unpaired_count += 1

        result = assign_lca_from_pair(pair_reads, refid_to_taxid, parent, rank, names, min_identity, extract_read_attributes)
        if result:
            output.append(result)
            if result[6] == "same":
                same_taxid_count += 1
            elif result[6] == "discordant":
                discordant_taxid_count += 1
                if write_discordant_bam:
                    for r in pair_reads:
                        discordant_bam.write(r)
            elif result[6] in ("R1 less specific", "R2 less specific"):
                less_specific_taxid_count +=1
                if write_specificity_bam:
                    for r in pair_reads:
                        specificity_bam.write(r)

    samfile.close()
    if discordant_bam:
        discordant_bam.close()
    if specificity_bam:
        specificity_bam.close()

    print(f"✅ Done: {len(output):,} read pairs assigned in {time.time() - start_time:.1f}s")
    print(f"🔢 Read groups with R1 + R2:   {paired_count:,}")
    print(f"🔢 Read groups missing one end: {unpaired_count:,}")
    print(f"🧬 Paired reads with same taxid:       {same_taxid_count:,}")
    print(f"🧬 Paired reads where one read is more specific:       {less_specific_taxid_count:,}")
    print(f"🧬 Paired reads with discordant taxid: {discordant_taxid_count:,}")
    return output

# --- Step 4: LCA read assignment ---
# (rest of code above unchanged)

# (rest of code above unchanged)

# (rest of code above unchanged)

def assign_lca_from_pair(reads, refid_to_taxid, parent, rank, names, min_identity=0.0, extract_read_attributes=False):
    from collections import defaultdict

    best_by_read = defaultdict(list)
    best_scores = {}
    read_metrics = {}

    for read in reads:
        ref_id = read.reference_id
        if ref_id not in refid_to_taxid:
            continue

        taxid = refid_to_taxid[ref_id]
        nm = read.get_tag("NM") if read.has_tag("NM") else 0
        query_len = read.query_length or 0
        softclip = sum(length for op, length in (read.cigartuples or []) if op == 4)
        aligned_len = query_len - softclip
        score = (aligned_len - nm) / query_len if query_len > 0 else 0

        if score < min_identity:
            continue

        key = "R1" if read.is_read1 else "R2" if read.is_read2 else "U"

        if extract_read_attributes:
            entry = (taxid, score, query_len, read)
        else:
            entry = (taxid, score, query_len)

        if key not in best_scores or score > best_scores[key]:
            best_by_read[key] = [entry]
            best_scores[key] = score
        elif score == best_scores[key]:
            best_by_read[key].append(entry)

    if not best_by_read:
        return None

    collapsed_taxids = {}
    total_bases = 0
    best_identity = max(best_scores.values()) if best_scores else 0

    for key, hits in best_by_read.items():
        taxid_list = [h[0] for h in hits]
        bases = sum(h[2] for h in hits)
        total_bases += bases

        if len(taxid_list) == 1:
            collapsed_taxids[key] = taxid_list[0]
        else:
            lca_taxid = taxid_list[0]
            for t in taxid_list[1:]:
                lca_taxid = find_lca(lca_taxid, t, parent)
            collapsed_taxids[key] = lca_taxid

        if extract_read_attributes:
            read = hits[0][3]
            nm = read.get_tag("NM") if read.has_tag("NM") else 0
            insertions = sum(length for op, length in (read.cigartuples or []) if op == 1)
            deletions = sum(length for op, length in (read.cigartuples or []) if op == 2)
            mismatches = max(0, nm - insertions - deletions)

            metrics = {
                "length": read.query_length or 0,
                "mismatches": mismatches,
                "insertions": insertions,
                "deletions": deletions,
                "softclips": sum(length for op, length in (read.cigartuples or []) if op == 4),
                "hardclips": sum(length for op, length in (read.cigartuples or []) if op == 5),
                "avg_qual": round(sum(read.query_qualities or []) / len(read.query_qualities or [1]), 2)
            }

            # For unpaired reads, store in R1 metrics
            read_metrics["R1" if key in ("R1", "U") else "R2"] = metrics

    taxid_r1 = collapsed_taxids.get("R1")
    taxid_r2 = collapsed_taxids.get("R2")
    taxid_u = collapsed_taxids.get("U")
    concordance = "discordant"

    if taxid_r1 and taxid_r2:
        if taxid_r1 == taxid_r2:
            final_taxid = taxid_r1
            concordance = "same"
        elif taxid_r1 in lineage(taxid_r2, parent):
            final_taxid = taxid_r1
            concordance = "R2 less specific"
        elif taxid_r2 in lineage(taxid_r1, parent):
            final_taxid = taxid_r2
            concordance = "R1 less specific"
        else:
            final_taxid = find_lca(taxid_r1, taxid_r2, parent)
            concordance = "discordant"
    else:
        final_taxid = taxid_r1 or taxid_r2 or taxid_u or 1
        concordance = "unpaired"

    tax_name = names.get(final_taxid, "Unknown")
    tax_rank = rank.get(final_taxid, "NA")

    row = [reads[0].query_name, final_taxid, tax_name, tax_rank, total_bases, round(best_identity, 4), concordance]

    if extract_read_attributes:
        def unpack(k):
            m = read_metrics.get(k, {})
            return [
                m.get("length", 0), m.get("mismatches", 0), m.get("insertions", 0), m.get("deletions", 0),
                m.get("softclips", 0), m.get("hardclips", 0), m.get("avg_qual", 0.0)
            ]
        row.extend(unpack("R1") + unpack("R2"))

    return row

# Updated assignment logic is assumed to already include the optional read attribute extraction

# --- Step 5: CLI ---
def main():
    parser = argparse.ArgumentParser(description="Assign taxonomic labels to BAM reads using LCA + specificity-based logic.")

    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--acc2taxid", required=True, help="Accession2taxid mapping file (.tsv or .gz)")
    parser.add_argument("--nodes", required=True, help="nodes.dmp taxonomy file (can be .gz)")
    parser.add_argument("--names", required=True, help="names.dmp taxonomy file (can be .gz)")
    parser.add_argument("--output", required=True, help="Output prefix - Assigned reads will be output as <output>.csv, and, "
    "if --estimate_abundance is used, abundances will be written to <output>.abundances")
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")
    parser.add_argument("--min_identity", type=float, default=0.0, help="Optional - Minimum identity threshold (default: 0.0)")
    parser.add_argument("--max_reads", type=int, help="Optional - maximum number of reads to process")
    parser.add_argument("--write_discordant_bam", action="store_true", help="Optional - Write BAM of discordant reads")
    parser.add_argument("--discordant_bam_path", default="discordant.bam", help="Optional - Path to write discordant BAM")
    parser.add_argument("--write_specificity_bam", action="store_true", help="Optional - Write BAM of specificity-overridden reads")
    parser.add_argument("--specificity_bam_path", default="specificity.bam", help="Optional - Path to write specificity BAM")
    parser.add_argument("--extract_read_attributes", action="store_true", help="Optional - Extract and output read length, " \
    "CIGAR, and quality metrics")
    parser.add_argument("--estimate_abundance", action="store_true", help="Optional - Estimate abundance by proportionally " \
    "propagating bases from ambiguous taxonomic levels to more specific taxonomic levels. This will output a file "
    "<output>.abundances with taxon names and assigned bases. By default, this done naively. If there are 3 species in a genus, " \
    "the reads will be allocated proportional to their species-level base abundances, regardless of their phylogenetic distance. " \
    "For phylogenetically-informed scaling, use --mash_reallocation")
    parser.add_argument("--genome_size_scaling", action="store_true", help="Optional - Adjust base abundances based on " \
    "reference genome lengths")
    parser.add_argument("--mash_reallocation", action="store_true", help="Optional - Use mash to calculate distances and " \
    "reallocate bases to more specific taxa, taking into account that more distant taxa are less likely to be misassigned")
    parser.add_argument("--reference_genome_list", default="references.list", help="Optional - Path to a list of reference genomes " \
    "to use for mash reallocation and adjusting base abundances based on reference genome lengths. Required if "
    "--genome_size_scaling and/or --mash_reallocation is used")
    parser.add_argument("--threads", type=int, default=4, help="Optional - Number of threads to use for mash " \
    "distance matrix estimation (default: 4)")

    args = parser.parse_args()

    results = process_bam_streaming(
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
                "ReadName", "TaxID", "TaxName", "Rank", "TotalBases", "BestIdentity", "ReadPairConcordance",
                "R1_Length", "R1_Matches", "R1_Insertions", "R1_Deletions", "R1_Softclips", "R1_Hardclips", "R1_AvgQual",
                "R2_Length", "R2_Matches", "R2_Insertions", "R2_Deletions", "R2_Softclips", "R2_Hardclips", "R2_AvgQual"
            ])
        else:
            writer.writerow([
                "ReadName", "TaxID", "TaxName", "Rank", "TotalBases", "BestIdentity", "ReadPairConcordance"
            ])
        writer.writerows(results)

    if args.estimate_abundance:
        assignments = pd.read_csv(args.output + ".csv")
        assignments['Bases'] = assignments.groupby(['TaxID', 'TaxName'])['TotalBases'].transform('sum')
        assignments = assignments.drop_duplicates(subset=['TaxID', 'TaxName']).reset_index(drop=True)
        observed_read_counts=dict(zip(assignments["TaxID"], assignments["Bases"]))
        taxid_list = assignments["TaxID"].tolist()
        if args.mash_reallocation:
            distmatrix, taxrenamemap_dists = mash_matrix.make_dist_matrix(args.reference_genome_list,
                                                                      args.output, 
                                                                      args.threads,
                                                                      args.acc2taxid)
        abundances, parenttochildren, named_dict = propagate_counts.propagate_counts(
            taxid_list=taxid_list, 
            nodes_path=args.nodes, 
            names_path=args.names,
            observed_read_counts=observed_read_counts
            )
        propagated_counts = pd.DataFrame(abundances.items())
        propagated_counts[2] = propagated_counts[1] / propagated_counts[1].sum()
        propagated_counts.columns = ["Taxon", "Assigned_Bases", "Proportion"]
        propagated_counts.to_csv(args.output + ".abundances", index=False)

if __name__ == "__main__":
    main()