#!/usr/bin/env bash

#run this script for each R1 fastq file obtained from amplicon sequencing individually.
#grep pattern to use for multi-v1 library: [ACTG]{3}GT[ACTG]{3}CT[ACTG]{3}AG[ACTG]{3}TG[ACTG]{3}CA[ACTG]{3}

if (($# != 3)); then
  echo invalid number of arguments, provide sample name followed by grep pattern and path to R1 fastq file
  exit
fi

grep -Po $2 $3 > $1_celltag_reads.txt
