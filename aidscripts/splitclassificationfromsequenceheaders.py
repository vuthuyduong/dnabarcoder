#!/usr/bin/env python
# FILE: addclassificationtosequenceheaders.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys, argparse
if sys.version_info[0] >= 3:
	unicode = str
import numpy as np
import os
from Bio import SeqIO
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='splitclassificationtosequenceheaders.py',  
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -o output",
							   description='''Script that remove classification from sequence headers. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be clustered.')
parser.add_argument('-o','--out', help='The fasta filename with sequence headers containing classification.') 
parser.add_argument('-c','--classification', required=True, help='the output file containing classification file in tab. format.')
parser.add_argument('-p','--classificationposition', default="", help='the classification positions of the features to kept for sequence headers.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
output=args.out
classificationpos=args.classificationposition

#fastafilename=sys.argv[1]
#classificationfilename=sys.argv[2]
#output=sys.argv[3]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

#####main###
poslist=[]
if "," in classificationpos:
	texts=args.classificationpos.split(",")
	for t in texts:
		poslist.append(int(t))
elif classificationpos!="":
	poslist.append(int(classificationpos))
if len(poslist)==0:
	poslist.append(0)
#seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
classificationfile=open(classificationfilename,"w")
outputfile=open(output,"w")
fastafile=open(fastafilename)
write=True
seqids=[]
for line in fastafile:
	if line.startswith(">"):
		write=True
		header=line.rstrip().replace(">","")
		classification=header
		seqid=header
		if "|" in header:
			classfication=header.replace("|","\t")
			texts=header.split("|")
			header=""
			seqid=texts[0]
			for pos in poslist:
				if header=="":
					seqid=texts[pos].replace(" ","")
					header=header + seqid + "|"
				else:
					header=header + texts[pos] + "|"
			classification=header.replace("|","\t")
			for i in range(0,len(texts)):
				if i not in poslist:
					classification=classification + texts[i] + "\t"
			header=header[:-1]	
			classification=classification[:-1]
		if not (seqid in seqids):
			seqids.append(seqid)
			classificationfile.write(classification + "\n")	
			outputfile.write(">" + header + "\n")	
		else:
			write=False
	else:
		if write==True:		
			outputfile.write(line)
outputfile.close()		
classificationfile.close()
fastafile.close()
