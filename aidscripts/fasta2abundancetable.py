#!/usr/bin/env python
# FILE: fasta2adundancetable.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 May 2021
import sys, argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser=argparse.ArgumentParser(prog='fasta2adundancetable.py',  
							   usage="%(prog)s [options] -i folderoffastafiles -o seqtable",
							   description='''Script that create an abundance table file from samples as fasta files.''',
							   epilog="""Written by Duong Vu d.vu@wi.knaw.nl""",
   )

parser.add_argument('-i','--input', required=True, help='the asv output of data.')
parser.add_argument('-o','--out', help='The output prefix for the results.') 
parser.add_argument('-minlen','--minsequencelength', type=int, default=0, help='The minimun sequence length for the selection.') 
parser.add_argument('-minoccurrence','--minoccurrence', type=int, default=0, help='The minimun number of occurrences.') 

args=parser.parse_args()
inputpath= args.input
output= args.out
minlen=args.minsequencelength
minoccurrence=args.minoccurrence

####################MAIN####################
seqtab={}
sequences=[]
for filename in os.listdir(inputpath):
	fullfilename=inputpath + "/" + filename
	samplename=filename.split(".")[0]
	print("Processing " + samplename + "...")
	seqtab.setdefault(samplename,{})
	seqrecords=list(SeqIO.parse(fullfilename, "fasta"))
	for seqrecord in seqrecords:
		seq=str(seqrecord.seq)
		if not (seq in sequences):
			sequences.append(seq)
		if  len(seq) <minlen: #only select sequences with sequence length >=minlen
			continue	
		if not (seq in seqtab[samplename].keys()):
			seqtab[samplename].setdefault(seq,0)
		seqtab[samplename][seq]	= seqtab[samplename][seq] + 1
selectedsequences={}
i=0
for seq in sequences:
	i=i+1
	count=0
	for samplename in seqtab.keys():
		try:
			occurence=seqtab[samplename][seq]
			count=count+occurence
		except KeyError:
			seqtab[samplename].setdefault(seq,0)
			pass
# 		if not (seq in seqtab[samplename].keys()):		
# 			seqtab[samplename].setdefault(seq,0)
#save abundance table
	if count>=minoccurrence:
		selectedsequences.setdefault("asv_" + str(i),seq)
outputfile=open(output + ".seq_as_columnnames.table","w")
outputfile2=open(output + ".table","w")
fastafile=open(output + ".fasta", "w")
header=""
header2=""
for seqid in selectedsequences.keys():
	seq=selectedsequences[seqid]
	header=header + "\t" + seq
	header2=header2 + "\t"  + seqid
	fastafile.write(">" + seqid +"\n")
	fastafile.write(seq + "\n")
outputfile.write(header + "\n")
outputfile2.write(header2 + "\n")
i=0
for samplename in seqtab.keys():
	line=samplename
	for seqid in selectedsequences.keys():
		seq=selectedsequences[seqid]
		line=line + "\t" + str(seqtab[samplename][seq])
	outputfile.write(line + "\n")	
	outputfile2.write(line + "\n")	
outputfile.close()
outputfile2.close()
fastafile.close()
print("Number of samples " + str(len(seqtab.keys())))
print("Number of amplicon sequence variants: " + str(len(selectedsequences)))
print("The abundance tables are saved in file " + output + ".table and " + output + ".seq_as_columnnames.table.")
print("The fasta file is saved in file " + output + ".fasta.")
			