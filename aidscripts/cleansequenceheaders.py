#!/usr/bin/env python
# FILE: cleansequenceheaders.py, requiring Python 2.7 if classification file contains non-ascii characters
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
import sys, argparse
import numpy as np
import os
from Bio import SeqIO
import json
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='cleansequenceheaders.py',  
							   usage="%(prog)s [options] -i fastafile -o output",
							   description='''Script that keeps unique sequence id from sequence headers. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file.')
parser.add_argument('-o','--out', help='The output fasta file with sequence headers being cleaned. Only the first part of the sequence headers separated by "|" will be kept.')

args=parser.parse_args()
fastafilename= args.input
output=args.out

#fastafilename=sys.argv[1]
#text=sys.argv[2]
#output=sys.argv[3]

outputfile=open(output,"w")
fastafile=open(fastafilename)

seqids=[]
counts=[]
for line in fastafile:
	if line.startswith(">"):
		header=line.rstrip()
		seqid=header
		if "|" in header:
			seqid=header.split("|")[0]	
		if seqid in seqids:
			counts[seqids.index(seqid)]=counts[seqids.index(seqid)]+1
			seqid=seqid + "_" + str(counts[seqids.index(seqid)])
		else:	
			seqids.append(seqid)
			counts.append(1)
		outputfile.write(seqid + "\n")	
	else:		
		outputfile.write(line)
fastafile.close()
outputfile.close()		
i=0
for seqid in seqids:
	if counts[i] >1:
		print("Sequence " + seqid.replace(">","") + " has " + str(counts[i]) + " sequences.")
	i=i+1

