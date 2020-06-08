#!/usr/bin/env python
# FILE: addclassificationtosequenceheaders.py, requiring Python 2.7 if classification file contains non-ascii characters
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
classificationfilename=sys.argv[2]
output=sys.argv[3]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadClassification(classificationfilename):
	classificationfile= open(classificationfilename)
	seqids=[]
	classifications=[]
	numberoffeatures=0
	for line in classificationfile:
		elements=line.rstrip().split("\t")
		seqid = elements[0].replace(">","").rstrip()
		seqids.append(seqid)
		classification=(line[line.index("\t")+1:]).replace("\t","|")
		classifications.append(classification)
		if numberoffeatures==0:
			numberoffeatures=classification.count("|")
		else:
			numberoffeatures=min(numberoffeatures,classification.count("|"))
	classificationfile.close()
	return seqids,classifications,numberoffeatures

#####main###

#seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
seqids,classifications,numberoffeatures=LoadClassification(classificationfilename)
outputfile=open(output,"w")
fastafile=open(fastafilename)

for line in fastafile:
	if line.startswith(">"):
		seqid=line.rstrip().replace(">","")
		if "|" in seqid:
			seqid=seqid[0:seqid.index("|")]
		classification=""
		if seqid in seqids:
			classification=unicode(classifications[seqids.index(seqid)],errors='ignore')
		else:
			classification="|" * numberoffeatures + "\n"
		header=">" + seqid + "|" + classification
		outputfile.write(header)	
	else:		
		outputfile.write(line)
	
outputfile.close()		

