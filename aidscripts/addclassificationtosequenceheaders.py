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
#parser.add_argument('-p','--classificationposition', default="", help='the classification positions for the prediction.')
parser.add_argument('-ranks','--classificationranks', default="", help='the classification ranks for the prediction, separated by ",".')
parser.add_argument('-fmt','--format', default="", help='the format of the header. There are two options for the fmt: empty or unite.')


args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
#classificationpos=args.classificationposition
classificationranks=args.classificationranks
fmt=args.format
output=args.out

#fastafilename=sys.argv[1]
#classificationfilename=sys.argv[2]
#output=sys.argv[3]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def unite(word,rank):
	word=word.rstrip()
	if rank.lower()=="species":
		word="s__"+word
	elif rank.lower()=="genus":
		word="g__"+word	
	elif rank.lower()=="family":
		word="f__"+word	
	elif rank.lower()=="order":
		word="o__"+word	
	elif rank.lower()=="class":
		word="c__"+word	
	elif rank.lower()=="phylum":
		word="p__"+word	
	elif rank.lower()=="kingdom":
		word="k__"+word	
	return word

def LoadClassification(classificationfilename,poslist,ranklist,fmt):
	classificationfile= open(classificationfilename, errors="ignore")
	seqids=[]
	classifications=[]
	numberoffeatures=0
	for line in classificationfile:
		texts=line.split("\t")
		seqid = texts[0].replace(">","").rstrip()
		seqids.append(seqid)
		classification=""
		i=0
		for pos in poslist:
			text=""
			if pos < len(texts):
				text=texts[pos].rstrip()
			if text=="":
				text="unidentified"	
			if fmt=="":	
				classification=classification + text + "|"
			else:
				text=text.replace(" ","_")	
				classification=classification + unite(text,ranklist[i]) + ";"
			i=i+1	
		classification=classification[:-1] 		
		classifications.append(classification)
		if numberoffeatures==0:
			numberoffeatures=classification.count("|")
		else:
			numberoffeatures=min(numberoffeatures,classification.count("|"))
	classificationfile.close()
	return seqids,classifications,numberoffeatures

def GetPositionList(classificationfilename,ranks):
	ranklist=[]	
	if "," in ranks:
		ranklist=ranks.split(",")
	elif ranks !="":
		ranklist.append(ranks)
	positionlist=[]
	classificationfile=open(classificationfilename)
	header=classificationfile.readline()
	header=header.rstrip()
	classificationfile.close()
	texts=header.split("\t")
	isError=False
	for rank in ranklist:
		if rank in texts:
			pos=texts.index(rank)
			positionlist.append(pos)
		else:
			print("The rank " + rank + " is not given in the classification." )
			isError=True
	return positionlist,ranklist,isError

#####main###

#seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
poslist,ranklist,isError=GetPositionList(classificationfilename,args.classificationranks)	
if isError==True:
	print("The given ranks/features do not exist in the header of the classification file.")
	sys.exit()
if len(ranklist)==0:
	print("Please specify the ranks/features to add to sequence headers.")	
	sys.exit()
seqids,classifications,numberoffeatures=LoadClassification(classificationfilename,poslist,ranklist,fmt)
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
