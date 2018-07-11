#!/bin/bash
path_1="$1"
path_2="$2"
path_3="$3"
python newfasta.py "$path_1" "$path_3"
blastn  -query new.fasta -subject "$path_2" -outfmt 10 -out blast.csv
python call_hd.py new.fasta "$path_2" blast.csv
python gaps.py result/mismatches/output.txt result/fasta/