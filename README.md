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
The first step of clone calling is to obtain reads containing the CellTag sequence from the single-cell alignment file. Currently, we only support CellRanger/CellRanger-ATAC but would be happy to support alternate single-cell pipelines per user request. The default workflow runs this analysis in batch mode, suitable for processing one or multiple bam files at once. In this a shell script spawns multiple jobs running the bam parsing R script in parallel, one for each sample. Alternatively, users can directly run the Rscipt, if they do not wish to use the batch mode. We have provided 2 versions of the shell script, one for running in bash and another for running on a slurm based cluster. 

To perform bam parsing, you first need to create a CSV config file.

| sample_id  | bam_file | cell_barcode | celltag_version |
| ------------- | ------------- | ------- | ---- |
| Content Cell  | Content Cell  | -- | -- |
| Content Cell  | Content Cell  | -- | -- |

Each row corresponds to one sample. Column descriptions:
 - `sample_id` : sample name (avoid any spaces)
 - `bam_file` : path to single-cell bam file
 - `cell_barcode` : path to identified filtered cell barcodes
 - `celltag_version`:  one of `8N-v1`,`8N-v2`,`8N-v3` or `multi-v1`. Please check section xx if you would like to run this script with a custom lineage barcode.  


Next, run `cloneCalling_scripts/sample <path to config file>`. This should perform bam parsing for each sample in the config file and store outputs in the `celltag_reads/` folder. Log files generated for each sample should be stored in the `logs/` folder.


## Processing CellTag reads to identify clones
