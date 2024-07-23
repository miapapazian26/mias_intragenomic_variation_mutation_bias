#!/bin/bash

# Description: ROC Analysis
# Usage: ./runRoc_parallel.sh <SPECIES> <INPUT_DIR> <OUTPUT_DIR>

SPECIES="candida_parapsilosis"
INFOLDER="fasta/revisit_cds_data"
SUFFIX=".max.cds"

INPUT="${INFOLDER}/${SPECIES}${SUFFIX}"
OUTPUT_BASE="Results/ConstMut"

#results saved to a unique subfolder using the species name
OUTPUT="${OUTPUT_BASE}/${SPECIES}"
mkdir -p "$OUTPUT"

nohup Rscript --vanilla Scripts/rocAnalysis.R -i "$INPUT" -o "$OUTPUT" -d 0 -s 1000 -a 20 -t 5 -n 8 --est_csp --est_phi --est_hyp --max_num_runs 2 --codon_table 1 &> "${OUTPUT}/roc_analysis.log" &
wait
