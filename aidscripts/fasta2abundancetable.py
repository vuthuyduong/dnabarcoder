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
parser.add_argument('-o','--out', help='The folder for the results.') 

args=parser.parse_args()
inputpath= args.input
output= args.out

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
		if not (seq in seqtab[samplename].keys()):
			seqtab[samplename].setdefault(seq,0)
		seqtab[samplename][seq]	= seqtab[samplename][seq] + 1
for seq in sequences:
	for samplename in seqtab.keys():
		if not (seq in seqtab[samplename].keys()):		
			seqtab[samplename].setdefault(seq,0)
#save abundance table
outputfile=open(output,"w")
header=""
for seq in sequences:
	header=header + "\t" + seq
outputfile.write(header + "\n")
for samplename in seqtab.keys():
	line=""
	for seq in sequences:
		line=line + "\t" + str(seqtab[samplename][seq])
	outputfile.write(line + "\n")	
outputfile.close()
print("Number of samples " + str(len(seqtab.keys())))
print("Number of amplicon sequence variants: " + str(len(sequences)))
print("The abundance table is saved in file " + output + ".")
			