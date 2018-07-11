from Bio import SeqIO
import sys
f_input=open(sys.argv[1],"r")
f_output=open("new.fasta","wb")
cut_length=int(sys.argv[2])
seqs=SeqIO.parse(f_input,"fasta")
for item in seqs:
    f_output.write("> " + item.id + "\n")
    if cut_length==0:
    	f_output.write(str(item.seq)+"\n")
    else:
    	f_output.write(str(item.seq[cut_length:-cut_length])+"\n")
f_input.close()
f_output.close()



