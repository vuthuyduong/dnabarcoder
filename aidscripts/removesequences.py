#!/usr/bin/env python
# FILE: removesequences.py
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

parser=argparse.ArgumentParser(prog='removesequences.py', 
							   usage="%(prog)s [options] -i fastafile -c classificationfile -t taxa -o output",
							   description='''Script that selects the sequences for of the given taxa. The taxon names in the taxa are separated by ","''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', help='The fasta output file containing the sequences of the given taxa.') #optional
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-t','--taxa', default="", help='the taxa for the selection, separated by ",". If taxa="", all the sequences with the ids given in the classification file will be removed from the input file.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
taxa=args.taxa
output=args.out

#fastafilename=sys.argv[1]
#taxa=sys.argv[2]
#classificationfilename = sys.argv[3]
#output=sys.argv[4]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def SelectSeqIds(classificationfilename,taxa):
	taxalist=[]
	if ";" in taxa:
		taxalist=taxa.split(";")
	elif taxa!="":
		taxalist.append(taxa)
	classificationfile= list(open(classificationfilename, "r"))
	seqids=[]
	for line in classificationfile:
		elements=line.rstrip().split("\t")
		seqid = elements[0].replace(">","").rstrip()
		if taxa!="":
			for taxonname in taxalist:
				if taxonname in elements:
					seqids.append(seqid)		
		else:
			seqids.append(seqid)			
	return seqids

#####main###
seqids=SelectSeqIds(classificationfilename,taxa)
seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
selectedrecords=[]
for seqrec in seqrecords:
	seqid=seqrec.id
	if "|" in seqid:
		seqid=seqid.split("|")[0]
	if not seqid in seqids:
		selectedrecords.append(seqrec)
#save to file:
SeqIO.write(selectedrecords,output,"fasta")

