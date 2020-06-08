#!/usr/bin/env python
# FILE: removesequences.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import numpy as np
import os
from Bio import SeqIO
import json
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils


fastafilename=sys.argv[1]
taxa=sys.argv[2]
classificationfilename = sys.argv[3]
output=sys.argv[4]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def SelectSeqIds(classificationfilename,taxa):
	taxalist=[]
	if ";" in taxa:
		taxalist=taxa.split(";")
	else:
		taxalist.append(taxa)
	classificationfile= list(open(classificationfilename, "r"))
	seqids=[]
	for line in classificationfile:
		elements=line.rstrip().split("\t")
		seqid = elements[0].replace(">","").rstrip()
		for taxonname in taxalist:
			if taxonname in elements:
				seqids.append(seqid)		
	return seqids

#####main###
seqids=SelectSeqIds(classificationfilename,taxa)
seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
selectedrecords=[]
for seqrec in seqrecords:
	if not seqrec.id in seqids:
		selectedrecords.append(seqrec)
#save to file:
SeqIO.write(selectedrecords,output,"fasta")

