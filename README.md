# newCloneCalling
Repository detailing updated allowlisting and clonal analysis workflow, name not final


# Overview
## Allowlisting plasmid library to obtain a list of allowed CellTag barcodes
CellTag constructs are made available as lentiviral plasmid libraries. While the theoretical diversity of CellTags is very high (65,536 for the original tags and ~68 billion for the 18N libraries) the real barcode diversity in the plasmid library is limited due to bottlenecks during the process of synthesis, cloning etc. To identify CellTag barcodes present in the plasmid library, we perform an allowlisting step. The detailed methodology for generating sequencing libraries for allowlisting is described in our paper [link]. The following text outlines the computational workflow for processing the sequencing data:
- Obtain Read 1 (R1) fastq files for the 2 replicates and parse celltag reads from each by running the following script: `allowlisting_scripts/parse_fq_allowlisting.sh <sample name> <grep pattern> <path to R1 fastq file>`
- Error correct identified CellTag barcodes using starcode (set distance threshold to 4): `allowlisting_scripts/starcode_collapse.sh <distance threshold> <path to output of parse_fq_allowlisting.sh>`
- Use `allowlisting_scripts/create_allowlist.ipynb` to identify list of barcodes present in both replicates to obtain the final allowlist

The allowlist for the multi-v1 library used in our paper [link] has been provided in this repo at `misc_files/18N-multi-v1-allowlist.csv`. The allowlists for 8N-v1,v2 and v3 libraries have been made available on addgene (Most abundant barcodes): https://www.addgene.org/pooled-library/morris-lab-celltag/

## Single-cell read alignment/CellRanger

## Parsing single-cell bam files to obtain CellTag reads
The first step of clone calling is to obtain reads containing the CellTag sequence from the single-cell alignment file. Currently, we only support outputs from CellRanger/CellRanger-ATAC but would be happy to support alternate single-cell pipelines per user request.

To perform bam parsing, you first need to create a CSV config file.

| sample_id  | bam_file | cell_barcode | assay | celltag_version |
| ------------- | ------------- | ------- | ---- | --- |
| sample1_RNA  | data/sample_1/outs/possorted_genome_bam.bam | data/sample_1/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | RNA | multi-v1 |
| sample1_ATAC  | data/sample_1_atac/outs/possorted_bam.bam | data/sample_1_atac/outs/filtered_peak_bc_matrix/barcodes.tsv | ATAC | multi-v1 |
| sample2_RNA  | data/sample_2/outs/possorted_genome_bam.bam | data/sample_2/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | RNA | multi-v1 |

Each row corresponds to one sample. Column descriptions:
 - `sample_id` : sample name (avoid any spaces)
 - `bam_file` : path to single-cell bam file
 - `cell_barcode` : path to identified filtered cell barcodes
 - `celltag_version`:  one of `8N-v1`,`8N-v2`,`8N-v3` or `multi-v1`. Please check section xx if you would like to run this script with a custom lineage barcode.  


Next, run `cloneCalling_scripts/bam_parsing.R <path to config file>`. This should perform bam parsing for each sample in the config file and store outputs in the `celltag_reads/` folder. Note: The script uses parallelization, it would greatly benefit from a higher processor/core count.


## Processing CellTag reads to identify clones
Once CellTag reads have been parsed for each sample, we perform additional filtering and clone calling. THe jupyter notebook `cloneCalling_scripts/celltag_analysis_single_assay.ipynb` can be run to perform clone calling. Each of the steps have been outlines below. As with the bam parsing scripts, the clone calling notebook has been designed to process multiple files at once. A user might need to process multiple samples together if they expect clonally related cells across samples. This could happen in the following cases:
- A single population of cells has been split across multiple single cell library preps/ ports of the 10x chip.
- Multiple samples have been collected from the same population of cells across time points and the user is interested in identifying clones both within and across time points.
DO NOT process CellTag data obtained from two independently CellTagged cell populations together.

The major steps of clone calling include:
 - Filter, error-correct and allowlist CellTag Reads.
 - Create the Cell x CellTag matrix, binarize, remove cells with too many and too few CellTags.
 - Compute the cell-cell jaccard similarity matrix and identify clones.
 - (Optional) Identify sparse clones and sub-cluster them into smaller high confidence clones.
 - (Optional) Save additional clone metrics.
 - Save the final clone table.
 
Assessing outputs from the clone calling notebook:
 - Sequencing saturation: An estimate of how deeply CellTags have been sequenced. A high sequencing saturation (> 60%) is preferred.
 - Clone table: This is a table that lists cells in clones and their respective clone IDs.
 - QC plots: Plots for various CellTag QC metrics
