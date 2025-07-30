import pysam
import pandas as pd
from collections import defaultdict

def split_bam_by_taxid(input_bam_path, outprefix, read_taxid_df, taxa):
    # Open input BAM
    df = pd.read_csv(read_taxid_df)
    df = df[df["TaxName"].isin(taxa)].reset_index(drop=True)
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