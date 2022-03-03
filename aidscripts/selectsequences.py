#!/usr/bin/env python
# FILE: selectsequences.py
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

parser=argparse.ArgumentParser(prog='selectsequences.py', 
							   usage="%(prog)s [options] -i fastafile -c classificationfile -t taxa -o output",
							   description='''Script that selects the sequences for of the given taxa. The taxon names in the taxa are separated by ","''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', help='The fasta output file containing the sequences of the given taxa.') #optional
parser.add_argument('-c','--classification', default="", help='the classification file in tab. format.')
parser.add_argument('-t','--taxa', default="", help='the taxa for the selection, separated by ","')
parser.add_argument('-n','--number', type=int, default=0, help='the maximum number of the sequences to be selected')
parser.add_argument('-p','--classificationposition', type=int, default=0, help='the classification positions for the selection.')
parser.add_argument('-l','--length', type=int, default=0, help='the required minimum length.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
taxa=args.taxa
output=args.out
n=args.number
l=args.length
classificationpos=int(args.classificationposition)

#fastafilename=sys.argv[1]
#taxa=sys.argv[2] #separated by ;
#classificationfilename = sys.argv[3]
#output=sys.argv[4]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def SelectSeqIds(classificationfilename,taxa,classificationpos):
	if not os.path.exists(classificationfilename):
		return [],[]
	taxalist=[]
	classnames=[]
	if "," in taxa:
		taxalist=taxa.split(",")
	elif taxa!="" and taxa!="unidentified":
		taxalist.append(taxa)
	classificationfile= list(open(classificationfilename, "r"))
	seqids=[]
	for line in classificationfile:
		elements=line.rstrip().split("\t")
		seqid = elements[0].replace(">","").rstrip()
		classname=""
		if classificationpos < len(elements):
			classname=elements[classificationpos]
		if classname=="" or classname=="unidentified":
			continue
		if taxa!="":
			for taxonname in taxalist:
				if taxonname!="" and taxonname!="unidentified":
					if taxonname in elements:
						seqids.append(seqid)				
						classnames.append(classname)
		else:
			seqids.append(seqid)
			classnames.append(classname)
	return seqids,classnames

#####main###
seqids,classnames=SelectSeqIds(classificationfilename,taxa,classificationpos)
seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
selectedrecords=[]
selectedclassnames={}
for seqrec in seqrecords:
	if seqrec.id in seqids or len(seqids)==0:
		if n==0:
			if len(str(seqrec.seq))>l:
				selectedrecords.append(seqrec)
		else:
			classname=classnames[seqids.index(seqrec.id)]
			if not classname in selectedclassnames.keys():
				selectedclassnames.setdefault(classname,1)
			if 	selectedclassnames[classname] <n:
				if len(str(seqrec.seq))>l:
					selectedclassnames[classname]=selectedclassnames[classname]+1
					selectedrecords.append(seqrec)
#save to file:
SeqIO.write(selectedrecords,output,"fasta")
print("The selected sequences are saved in " + output + ".")

