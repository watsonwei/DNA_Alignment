import csv
from Bio import SeqIO
import os
import shutil
import sys
print sys.argv[1]
f=open(sys.argv[1],"r")
if not os.path.exists("new_result"):
    os.makedirs("new_result")
else:
    shutil.rmtree("new_result")
    os.makedirs("new_result")
if not os.path.exists("new_result/mismatches"):
    os.makedirs("new_result/mismatches")
else:
    shutil.rmtree("new_result/mismatches")
    os.makedirs("new_result/mismatches")
if not os.path.exists("new_result/fasta"):
    os.makedirs("new_result/fasta")
else:
    shutil.rmtree("new_result/fasta")
    os.makedirs("new_result/fasta")
if not os.path.exists("new_result/alignment"):
    os.makedirs("new_result/alignment")
else:
    shutil.rmtree("new_result/alignment")
    os.makedirs("new_result/alignment")
f_output=open("new_result/mismatches/new_output.txt","wb")
count=0
lambda_list= list(SeqIO.parse(open("lambda.fasta"), "fasta"))
s3=str(lambda_list[0].seq)
for line in f:
    line=line.strip().split()
    print line
    count += 1
    print count
    print line
    file = sys.argv[2]+ line[0] + ".fasta"
    records = list(SeqIO.parse(open(file), "fasta"))
    s1 = str(records[0].seq)
    s2 = str(records[1].seq)
    left_cut_index = 0
    right_cut_index = len(s1)
    l_gaps = 0
    r_gaps = 0
    for i in range(10, len(s1), 10):
        for j in range(i - 10, i):
            if s1[j] == '-':
                l_gaps += 1
            elif s2[j] == '-':
                l_gaps += 1
        if l_gaps > 5:
            left_cut_index = i
            l_gaps = 0
        else:
            l_gaps = 0
            for j in range(max(i - 20, 0), i + 10):
                if s1[j] == '-':
                    l_gaps += 1
                elif s2[j] == '-':
                    l_gaps += 1
            if l_gaps > (i + 10 - max(i - 20, 0)) / 2.0:
                left_cut_index = i
                l_gaps = 0
            else:
                while s1[left_cut_index]=="-" or s2[left_cut_index]=="-":
                    left_cut_index+=1
                break
    for m in range(len(s1) - 10, 0, -10):
        for n in range(m, m + 10):
            if s1[n] == '-':
                r_gaps += 1
            elif s2[n] == '-':
                r_gaps += 1
        if r_gaps > 5:
            right_cut_index = m
            r_gaps = 0
        else:
            r_gaps = 0
            for n in range(m - 10, min(m + 20, len(s1))):
                if s1[n] == '-':
                    r_gaps += 1
                elif s2[n] == '-':
                    r_gaps += 1
            if r_gaps > (min(m + 20, len(s1)) - (m - 10)) / 2.0:
                right_cut_index = m
                r_gaps = 0
            else:
                while s1[right_cut_index-1]=="-" or s2[right_cut_index-1]=="-":
                    right_cut_index-=1
                break
    print "left_cut_index", left_cut_index
    print "right_cut_index", right_cut_index
    s1 = s1[left_cut_index:right_cut_index]
    s2 = s2[left_cut_index:right_cut_index]
    f_fasta = open("new_result/fasta/" + line[0] + ".fasta", "wb")
    f_fasta.write(">" + line[0] + "\n")
    f_fasta.write(s1 + "\n")
    f_fasta.write(">NC_001416.1" + "\n")
    f_fasta.write(s2 + "\n")
    mismatches = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            continue
        elif s1[i] == "-":
            mismatches += 1
        elif s2[i] == "-":
            mismatches += 1
        else:
            mismatches += 1
    new_p = float(len(s1) - mismatches) / float(len(s1))
    f_output.write(line[0] + " " + str(len(s1.replace("-", ""))) + " " + str(mismatches) + " " + str(len(s1)) + " " + line[-1] + " " + str(new_p) + "\n")
    f_txt = open("new_result/alignment/" + line[0] + ".txt", "wb")
    for i in range(0,len(s1),50):
        f_txt.write(s2[i:i+50].replace("-"," ")+"\n")
        for j in range(i,min(i+50,len(s1))):
            if s2[j]==s1[j]:
                f_txt.write("|")
            elif s2[j]=="-":
                f_txt.write("-")
            elif s1[j]=="-":
                f_txt.write("-")
            else:
                f_txt.write(" ")
        f_txt.write("\n")
        f_txt.write(s1[i:i+50].replace("-"," ")+"\n")