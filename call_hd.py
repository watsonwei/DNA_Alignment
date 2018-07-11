#!/usr/bin/env python2
from hd import hamming_distance
from Bio import SeqIO
import csv
import os
import shutil
import sys
def readFasta(file):
    records=list(SeqIO.parse(open(file),"fasta"))
    return records
test_file=sys.argv[1]
reference_file=sys.argv[2]
test_list=readFasta(test_file)
reference_list=readFasta(reference_file)

if not os.path.exists("result"):
    os.makedirs("result")
else:
    shutil.rmtree("result")
    os.makedirs("result")
if not os.path.exists("result/alignment"):
    os.makedirs("result/alignment")
else:
    shutil.rmtree("result/alignment")
    os.makedirs("result/alignment")
if not os.path.exists("result/mismatches"):
    os.makedirs("result/mismatches")
else:
    shutil.rmtree("result/mismatches")
    os.makedirs("result/mismatches")
f_output=open("result/mismatches/output.txt","wb")


def complimentary_seq(reference_seq):
    c_reference_seq = ""
    for c in reference_seq:
        if c == "A":
            c_reference_seq += "T"
        if c == "T":
            c_reference_seq += "A"
        if c == "G":
            c_reference_seq += "C"
        if c == "C":
            c_reference_seq += "G"
    return c_reference_seq

# seq is the query sequence
# seq_id is the id of query sequence
# reference_seq is the reference sequence
# qstart is the query start index in blast file
# qend is the query end index in blast file
# sstart is the subject start index in blast file
# send is the subject end index in blast file


def call_hd(seq,seq_id,reference_seq,qstart,qend,sstart,send):
    length = len(seq)
    if sstart < send:
        left_index = max(sstart -1- (qstart-1), 0)
        right_index = send + (length - qend)
        blast_sequence = str(reference_seq[left_index:right_index])
    if sstart > send:
        right_index = sstart + (qstart-1)
        left_index = max(send - 1-(length - qend), 0)
        c_reference_seq=complimentary_seq(reference_seq)
        blast_sequence = c_reference_seq[left_index:right_index]
        blast_sequence = blast_sequence[::-1]
    distance = hamming_distance(blast_sequence, seq, seq_id)
    # f_output.write(seq_id + "  " + str(length) + "  " + str(distance) + "\n")
    # print seq_id + "  " + str(length) + "  " + str(distance)
    print distance
    return seq_id, length, distance


#create result folder to store the sequence alignment txt

# read from csv
list_csv=[]
f_csv=open(sys.argv[3],"r")
reader = csv.reader(f_csv)
for row in reader:
    list_csv.append(row)
f_csv.close()

dic={} #dic after eliminating the duplicates in csv based on the match percentage
for item in list_csv:
    if item[0] not in dic:
        dic[item[0]]=item
    else:
        if float(item[-1])>float(dic[item[0]][-1]):
            dic[item[0]]=item

dic_output={} #store the output statistics
for item in test_list:
    if item.id in dic:
        qstart=int(dic[item.id][6])
        qend=int(dic[item.id][7])
        sstart=int(dic[item.id][-4])
        send=int(dic[item.id][-3])
        print item.id
        seq_id,length,distance=call_hd(str(item.seq), item.id, str(reference_list[0].seq), qstart, qend, sstart,send)
        dic_output[seq_id]=[length,distance]

if not os.path.exists("result/fasta"):
    os.makedirs("result/fasta")
else:
    shutil.rmtree("result/fasta")
    os.makedirs("result/fasta")
 #generate fasta file   
for item in test_list:
    filename="result/alignment/"+item.id+".txt"
    try:
        f = open(filename, "r")
        lines = f.readlines()
        r_string = ""
        m_string = ""
        for i in range(len(lines)):
            line = lines[i].replace(" ", "-").strip()
            if i % 3 == 0:
                m_string = m_string + line
            if i % 3 == 2:
                r_string = r_string + line
        filename_output = "result/fasta/" + item.id + ".fasta"
        f1 = open(filename_output, "wb")
        f1.write("> " + item.id + "\n")
        f1.write(r_string)
        f1.write("\n")
        f1.write("> NC_001416.1" + "\n")
        f1.write(m_string)
        match_percentage=(len(m_string)-dic_output[item.id][1])/float(len(m_string))
        print item.id + "  " + str(dic_output[item.id][0]) + "  " + str(dic_output[item.id][1]) + "  "+ str(len(m_string))+"  "+str(match_percentage)
        f_output.write(item.id + "  " + str(dic_output[item.id][0]) + "  " + str(dic_output[item.id][1]) + "  "+ str(len(m_string))+"  "+str(match_percentage)+"\n")
        
    except Exception as e:
        continue






