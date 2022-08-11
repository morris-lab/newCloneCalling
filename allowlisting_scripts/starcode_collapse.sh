#!/usr/bin/env bash

if (($# != 2)); then
  echo invalid number of arguments, please provide distance threshold (recommended value is 4) followed by file path to the output of celltag parsing script
  exit
fi

/scratch/smlab/starcode/starcode -d $1 -t 4 -s -i $2 -o ${2:0:-4}_dist_$1.txt

