# newCloneCalling
Repository detailing updated allowlisting and clonal analysis workflow, name not final


# Overview
## Allowlisting plasmid library to obtain a list of allowed CellTag barcodes
CellTag constructs are made available as lentiviral plasmid libraries. While the theoretical diversity of CellTags is very high (65,536 for the original tags and ~68 billion for the 18N libraries) the real barcode diversity in the plasmid library is limited due to bottlenecks during the process of synthesis, cloning etc. To identify CellTag barcodes present in the plasmid library, we perform an allowlisting step. The detailed methodology for generating sequencing libraries for allowlisting is described in our paper [link]. The following text outlines the computational workflow for processing the sequencing data:
- Obtain Read 1 (R1) fastq files for the 2 replicates and parse celltag reads from each by running the following script: `allowlisting_scripts/parse_fq_allowlisting.sh <sample name> <grep pattern> <path to R1 fastq file>`
- Error correct identified CellTag barcodes using starcode (set distance threshold to 4): `allowlisting_scripts/starcode_collapse.sh <distance threshold> <path to output of parse_fq_allowlisting.sh>`
- Use `allowlisting_scripts/create_allowlist.ipynb` to identify list of barcodes present in both replicates to obtain the final allowlist

The allowlist for the multi-v1 library used in our paper [link] has been provided in this repo at `misc_files/18N-multi-v1-allowlist.csv`. The allowlists for 8N-v1,v2 and v3 libraries have been made available on addgene (Most abundant barcodes): https://www.addgene.org/pooled-library/morris-lab-celltag/

## Parsing single-cell bam files to obtain CellTag reads
The first step of clone calling is to obtain reads containing the CellTag sequence from the single-cell alignment file generated using CellRanger/CellRanger-ATAC. There are 2 modes for this analysis. In the single sample mode, users can directly run the provided Rscript with the required command line arguments. For processing multiple samples at once, we have provided a simple shell script that spawns multiple jobs in parallel, one for each sample to be analyzed. The shell script is compatible with a slurm based job submission system, please modify it as needed, to run with a different compute cluster manager.

### Single sample mode
 - Run `cloneCalling_scripts/sample` with the following command line arguments:
   - Sample name (avoid spaces, this will be prepended to the final output file)
   - Path to cellranger bam file
   - Path to list of valid cell barcodes as identified by cellranger (usually in `outs/filtered_feature_bc_matrix/barcodes.tsv.gz`)
   - single-cell assay: only `rna` and `atac` are currently supported
   - CellTag version to parse: one of `8N-v1`,`8N-v2`,`8N-v3` or `multi-v1` 
   - The clone calling pipeline allows for processing lineage barcodes outside of these 4 standard libraries. To process data from a custom lineage barcode, users can skip the CellTag version argument and instead specify the following three arguments:
     - Lineage barcode grep pattern
     - Expected length of lineage barcode
     - additional secondary grep


## Processing CellTag reads to identify clones
