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

parser=argparse.ArgumentParser(prog='addclassificationtosequenceheaders.py',  
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -o output",
							   description='''Script that adds taxonomic classification to sequence headers. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be clustered.')
parser.add_argument('-o','--out', help='The fasta filename with sequence headers containing classification.') 
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-p','--classificationposition', default="", help='the classification positions for the prediction.')


args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
classificationpos=args.classificationposition
output=args.out

#fastafilename=sys.argv[1]
#classificationfilename=sys.argv[2]
#output=sys.argv[3]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadClassification(classificationfilename,classificationpos):
	classificationfile= open(classificationfilename, errors="ignore")
	poslist=[]
	if "," in classificationpos:
		texts=classificationpos.split(",")
		for t in texts:
			poslist.append(int(t))
	elif classificationpos!="":
		poslist.append(int(classificationpos))	
	seqids=[]
	classifications=[]
	numberoffeatures=0
	for line in classificationfile:
		elements=line.rstrip().split("\t")
		seqid = elements[0].replace(">","").rstrip()
		seqids.append(seqid)
		classification=(line[line.index("\t")+1:])
		texts=classification.split("\t")
		classification=""
		if len(poslist) >0:		
			for pos in poslist:
				pos=pos-1
				text=""
				if pos < len(texts):
					text=texts[pos]
				if text=="":
					text="unidentified"	
				text=text.replace(" ","_")	
				classification=classification + text + "|"
		else:
			for text in texts:
				text=text.rstrip()
				if text=="":
					text="unidentified"
				text=text.replace(" ","_")		
				classification=classification + text + "|"
		classification=classification[:-1] 		
		classifications.append(classification)
		if numberoffeatures==0:
			numberoffeatures=classification.count("|")
		else:
			numberoffeatures=min(numberoffeatures,classification.count("|"))
	classificationfile.close()
	return seqids,classifications,numberoffeatures

#####main###

#seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
seqids,classifications,numberoffeatures=LoadClassification(classificationfilename,classificationpos)
outputfile=open(output,"w")
fastafile=open(fastafilename)

for line in fastafile:
	if line.startswith(">"):
		seqid=line.rstrip().replace(">","")
		if "|" in seqid:
			seqid=seqid[0:seqid.index("|")]
		classification=""
		if seqid in seqids:
			classification=unicode(classifications[seqids.index(seqid)])
		else:
			classification="|" * numberoffeatures
		header=">" + seqid + "|" + classification + "\n"
		outputfile.write(header)	
	else:		
		outputfile.write(line)
	
outputfile.close()		
fastafile.close()
