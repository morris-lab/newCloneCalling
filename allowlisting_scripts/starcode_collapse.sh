#!/usr/bin/env bash

if (($# != 2)); then
  echo invalid number of arguments, please
  exit
fi

/scratch/smlab/starcode/starcode -d $1 -t 4 -s -i $2 -o dist$1_$2

