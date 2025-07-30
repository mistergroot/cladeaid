import pysam
import pandas as pd
from collections import defaultdict
from collections.abc import Iterable

def split_bam_by_taxid(input_bam_path, outprefix, read_taxid_df, taxa):
    # Open input BAM
    df = pd.read_csv(read_taxid_df)
    if type(taxa) is list:
        df = df[df["TaxName"].isin(taxa)].reset_index(drop=True)
    elif type(taxa) is str:
        df = df[df["TaxName"] == taxa].reset_index(drop=True)
    elif type(taxa) is dict:
        taxnames = (
            [x for k, v in taxa.items() for x in 
             ([k] + (list(v) if isinstance(v, Iterable) and not 
                     isinstance(v, (str, bytes)) else [v]))]
        )
        df = df[df["TaxName"].isin(taxnames)].reset_index(drop=True)
    df = df.set_index("ReadName")
    infile = pysam.AlignmentFile(input_bam_path, "rb")

    # Store BAM writers per taxid
    bam_writers = {}

    # Loop through reads
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

        # Write to corresponding BAM file
        bam_writers[taxid].write(read)

    # Close all writers
    for writer in bam_writers.values():
        writer.close()

    infile.close()