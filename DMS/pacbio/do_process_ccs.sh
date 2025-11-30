#!/bin/bash

ncpu=1
targ=SARS-CoV-2-BF7
runs=pbruns.csv

IN=../ccsfastq
OUT=../outputs

python process_ccs.py --seqsfile=BF7_reference.gb \
                      --configfile=parse_specs.yaml \
                      --threads=$ncpu \
                      --pbruns=$runs \
                      --output=$OUT \
                      --target=$targ