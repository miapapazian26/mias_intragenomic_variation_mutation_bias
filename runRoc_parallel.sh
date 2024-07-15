#!/bin/bash
INPUT="fasta/revisit_cds_data/lachancea_kluyveri.max.cds"
OUTPUT="Results/ConstMut"
nohup Rscript --vanilla Scripts/rocAnalysis.R -i "$INPUT" -o "$OUTPUT" -d 0 -s 20000 -a 20 -t 5 -n 8 --est_csp --est_phi --est_hyp --max_num_runs 2 --codon_table 1 &> Results/ConstMut/$1.Rout &
wait
