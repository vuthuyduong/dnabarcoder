#!/usr/bin/env python
# FILE: addclassificationtosequenceheaders.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys, argparse
if sys.version_info[0] >= 3:
	unicode = str
# import numpy as np
# import os
from Bio import SeqIO
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='addclassificationtosequenceheaders.py',  
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -o output",
							   description='''Script that adds taxonomic classification to sequence headers. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be clustered.')
parser.add_argument('-o','--out', help='The fasta filename with sequence descriptions containing classification.') 
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-ranks','--classificationranks', default="", help='the classification ranks for the prediction, separated by ",".')
parser.add_argument('-sep','--separator', default="", help='the separator between the ranks in the header. If sep is not given, the sequence description will have the format: k__kingdom;p__phylum;c__class;o__order;f__family;g__genus;s__species.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')


args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
#classificationpos=args.classificationposition
classificationranks=args.classificationranks
sep=args.separator
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

def LoadClassification(classificationfilename,poslist,ranklist,seqidpos,sep):
	classificationfile= open(classificationfilename, errors="ignore")
	classificationdict={}
	numberoffeatures=0
	for line in classificationfile:
		texts=line.split("\t")
		seqid = texts[seqidpos].rstrip()
		classification=""
		i=0
		for pos in poslist:
			text=""
			if pos >-1 and pos < len(texts):
				text=texts[pos].rstrip()
			if text=="":
				text="unidentified"	
			if sep!="":	
				classification=classification + text + sep
			else:
				text=text.replace(" ","_")	
				classification=classification + unite(text,ranklist[i]) + ";"
			i=i+1	
		classification=classification[:-1] 		
		if numberoffeatures==0:
			numberoffeatures=classification.count(sep)
		else:
			numberoffeatures=min(numberoffeatures,classification.count(sep))
		classificationdict.setdefault(seqid,classification)	
	classificationfile.close()
	return classificationdict,numberoffeatures

def GetPositionList(classificationfilename,ranks,idcolumnname):
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
	i=0
	seqidpos=0
	for text in texts:
		if text.lower()==idcolumnname.lower():
			seqidpos=i
		i=i+1	
	for rank in ranklist:
		pos=-1
		if rank in texts:
			pos=texts.index(rank)
		positionlist.append(pos)
	return seqidpos,positionlist,ranklist,isError

#####main###

#seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
seqidpos,poslist,ranklist,isError=GetPositionList(classificationfilename,args.classificationranks,args.idcolumnname)	

classificationdict,numberoffeatures=LoadClassification(classificationfilename,poslist,ranklist,seqidpos,sep)

seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
for seqrecord in seqrecords:
	classification=""
	if seqrecord.id in classificationdict.keys():
		classification=unicode(classificationdict[seqrecord.id])
	elif sep!="":
		classification=sep * numberoffeatures
	seqrecord.description= classification
#save new file	
SeqIO.write(seqrecords, output, "fasta")	
