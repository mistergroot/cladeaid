import pysam
import pandas as pd
from collections import defaultdict
from collections.abc import Iterable
import subprocess
from . import tax_parsing

def split_bam_by_taxid(input_bam_path, outprefix, read_taxid_df, taxa):
    # Open input BAM
    df = pd.read_csv(read_taxid_df)
    if type(taxa) is list:
        df = df[df["TaxName"].isin(taxa)].reset_index(drop=True)
    elif type(taxa) is str:
        df = df[df["TaxName"] == taxa].reset_index(drop=True)
    elif type(taxa) is dict:
        taxnames = list(taxa.keys())
        df = df[df["TaxID"].isin(taxnames)].reset_index(drop=True)
    df = df.set_index("ReadName")

    infile = pysam.AlignmentFile(input_bam_path, "rb")

    # Store BAM writers per taxid
    bam_writers = {}

    # Loop through reads
    if type(taxa) is dict:
        df["UserTaxLevel"] = df["TaxID"].replace(taxa)
        for read in infile.fetch(until_eof=True):
            read_name = read.query_name

            if read_name not in df.index:
                continue  # skip reads not in the taxID mapping

            taxid = str(df.loc[read_name, "UserTaxLevel"]).replace(" ", "_")

            if taxid not in bam_writers:
                # Open a new BAM file for this taxID
                outpath = f"{outprefix}_{taxid}.bam"
                bam_writers[taxid] = pysam.AlignmentFile(outpath, "wb",
                                                        template=infile)
            
            if read.reference_name not in df.loc[read_name, "Contig"]:
                continue

            # Write to corresponding BAM file
            bam_writers[taxid].write(read)

    else:
        for read in infile.fetch(until_eof=True):
            read_name = read.query_name

            if read_name not in df.index:
                continue  # skip reads not in the taxID mapping

            taxid = str(df.loc[read_name, 'TaxName']).replace(" ", "_")

            if taxid not in bam_writers:
                # Open a new BAM file for this taxID
                outpath = f"{outprefix}_{taxid}.bam"
                bam_writers[taxid] = pysam.AlignmentFile(outpath, "wb",
                                                        template=infile)

            if read.reference_name not in df.loc[read_name, "Contig"]:
                continue
            # Write to corresponding BAM file
            bam_writers[taxid].write(read)

    # Close all writers
    for writer in bam_writers.values():
        writer.close()

    infile.close()

    if type(taxa) is list:
        outputs_to_sort = list(set(taxa))
    elif type(taxa) is str:
        outputs_to_sort = [taxa]
    elif type(taxa) is dict:
        outputs_to_sort = list(set(list(taxa.values())))
    
    for tax in outputs_to_sort:
        input_bam_file = f"{outprefix}_{tax.replace(" ", "_")}.bam"
        output_sorted_bam_file = f"{outprefix}_{tax.replace(" ", "_")}.sorted.bam"

        # Sort the BAM file by genomic position
        pysam.sort("-o", output_sorted_bam_file, input_bam_file)

        # Optionally, create an index for the sorted BAM file
        pysam.index(output_sorted_bam_file)

        subprocess.run(["rm", input_bam_file])