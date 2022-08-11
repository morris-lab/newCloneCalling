# newCloneCalling
Repository detailing updated allowlisting and clonal analysis workflow, name not final


# Overview
## Allowlisting plasmid library to obtain a list of allowed CellTag barcodes
CellTag constructs are made available as lentiviral plasmid libraries. While the theoretical diversity of CellTags is very high (65,536 for the original tags and ~68 billion for the 18N libraries) the real barcode diversity in the plasmid library is limited due to bottlenecks during the process of synthesis, cloning etc. To identify CellTag barcodes present in the plasmid library, we perform an allowlisting step. The detailed methodology for generating sequencing libraries for allowlisting is described in our paper [link]. The following text outlines the computational workflow for processing the sequencing data:
- Obtain Read 1 (R1) fastq files for the 2 replicates and parse celltag reads from each using: `allowlisting_scripts/parse_fq_allowlisting.sh`
- Error correct identified CellTag barcodes using starcode (distance threshold is set to 4): `allowlisting_scripts/starcode_collapse.sh`
- Use `allowlisting_scripts/notebook` to identify list of barcodes in both duplicates to obtain the final allowlist
The allowlist for the multi-v1 library used in our paper [link] has been provided in this repo.

## Parsing single-cell bam files to obtain CellTag reads

## Processing CellTag reads to identify clones
