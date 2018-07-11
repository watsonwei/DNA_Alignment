# DNA_Alignment
## description
- the results will be stored in result folder
- fasta files in the fasta folder
- alignments files in the alignment folder
- count of length, mismatches in the mismatches folder
- new fasta files will be stored in the newfasta folder
- new alignment files will be stored in txt folder
- new length, new mismatches, new alignment length, old percentage, new percentage will be stored in newoutput.txt


## Step 1:
Install blast
Install python 2.X

## Step 2:
cd to the current folder

## Step 3:
Run the command line:
`$ sh run.sh arg_1 arg_2 arg_3`

arg_1: directory of test sequence fasta file
arg_2: directory of reference sequence fasta file
arg_3: cut length(could be 0 when no cutting)

the columns in output.txt refers to id, length, mismatches(including gaps)
