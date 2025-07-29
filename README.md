[![logo](https://raw.githubusercontent.com/mistergroot/cladeaid/refs/heads/main/cladeaid_logo.png)](https://github.com/mistergroot/cladeaid)

Installation

```
git clone https://github.com/mistergroot/cladeaid
cd cladeaid
conda env create -f environment.yml
pip install -e .
```

Usage

```
$cladeaid --help
usage: cladeaid [-h] {lca,refinery} ...

LCA assignment and downstream refinement tools

positional arguments:
  {lca,refinery}
    lca           Run LCA assignment
    refinery      Run postprocessing/refinement

options:
  -h, --help      show this help message and exit
```

```
$cladeaid lca --help
usage: cladeaid lca [-h] --bam BAM --acc2taxid ACC2TAXID --nodes NODES --names NAMES --output OUTPUT [--verbose]
                    [--min_identity MIN_IDENTITY] [--max_reads MAX_READS] [--write_discordant_bam]
                    [--discordant_bam_path DISCORDANT_BAM_PATH] [--write_specificity_bam]
                    [--specificity_bam_path SPECIFICITY_BAM_PATH] [--extract_read_attributes]

options:
  -h, --help            show this help message and exit
  --bam BAM             Input BAM file
  --acc2taxid ACC2TAXID
                        Accession2taxid mapping file (.tsv or .gz)
  --nodes NODES         nodes.dmp taxonomy file (can be .gz)
  --names NAMES         names.dmp taxonomy file (can be .gz)
  --output OUTPUT       Output prefix
  --verbose             Enable verbose output
  --min_identity MIN_IDENTITY
                        Minimum identity threshold (default: 0.0)
  --max_reads MAX_READS
                        Maximum number of reads to process
  --write_discordant_bam
                        Write BAM of discordant reads
  --discordant_bam_path DISCORDANT_BAM_PATH
                        Path to write discordant BAM
  --write_specificity_bam
                        Write BAM of specificity-overridden reads
  --specificity_bam_path SPECIFICITY_BAM_PATH
                        Path to write specificity BAM
  --extract_read_attributes
                        Extract and output read attributes
```

```
$cladeaid refinery --help
usage: cladeaid refinery [-h] --csv CSV --acc2taxid ACC2TAXID --nodes NODES --names NAMES --output OUTPUT
                         [--genome_size_scaling] [--output_mappings] [--mash_reallocation]
                         [--reference_genome_list REFERENCE_GENOME_LIST] [--pairwise_dists PAIRWISE_DISTS]
                         [--multifasta] [--normalize] [--threads THREADS]

Estimate abundance using naive or mash distance-aware propagation. Optionally, normalize based on genome size, and
output bam file for each 'confident' taxon.

options:
  -h, --help            show this help message and exit
  --csv CSV             Input CSV from cladeaid.py
  --acc2taxid ACC2TAXID
                        Accession2taxid mapping file (can be .gz)
  --nodes NODES         nodes.dmp taxonomy file (can be .gz)
  --names NAMES         names.dmp taxonomy file (can be .gz)
  --output OUTPUT       Output prefix - Assigned reads will be output as <output>.csv, and abundances will be written
                        to <output>.abundances
  --genome_size_scaling
                        Adjust base abundances based on reference genome lengths
  --output_mappings     Outputs a bam for each species with confident mappings
  --mash_reallocation   Use mash to calculate distances and reallocate bases to more specific taxa
  --reference_genome_list REFERENCE_GENOME_LIST
                        Path to list of reference genomes for mash reallocation and genome size adjustment
  --pairwise_dists PAIRWISE_DISTS
                        Path to pairwise distance file between genomes for mash reallocation
  --multifasta          Flag if BAMs were mapped against a multi-organism fasta reference
  --normalize           Normalize bases assigned to taxa based on genome sizes
  --threads THREADS     Number of threads for mash distance matrix estimation (default: 4)
```

To run the test files included in `/examples/`, navigate to the examples directory and run the following:
```
cladeaid lca --bam MixA_bowtie_global.sorted.bam --acc2taxid acc2taxid.gz --names names.dmp.gz --nodes nodes.dmp.gz --output test
cladeaid refinery --csv test.csv --nodes nodes.dmp.gz --names names.dmp.gz --acc2taxid acc2taxid.gz --output test --mash_reallocation --reference_genome_list references.list
```
This should output the files: `test.csv`, `test.abundances`, and `test.dists`.
