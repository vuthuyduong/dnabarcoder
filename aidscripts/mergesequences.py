#!/usr/bin/env python
# FILE: mergesequences.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020
import sys
import numpy as np
import os, argparse
from Bio import SeqIO
import json
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='mergesequences.py', 
							   usage="%(prog)s [options] -i fastafiles -o output",
							   description='''Script that selects the sequences for of the given taxa. The taxon names in the taxa are separated by ","''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta files to be merged')
parser.add_argument('-o','--out', required=True, help='The fasta output file containing the sequences with unique sequence ids') #optional

args=parser.parse_args()
output=args.out

fastafilenames=args.input.split(",")
seqids=[]
mergedrecords=[]
count=0
for fastafilename in fastafilenames:
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	for seqid in seqrecords.keys():
		if not (seqid in seqids):
			if count >0 and seqid.startswith("K"):
				print(seqid)
			seqids.append(seqid)
			mergedrecords.append(seqrecords[seqid])
	count=count+1		
#save to file:
SeqIO.write(mergedrecords,output,"fasta")
print("All sequences are saved in " + output + ".")

