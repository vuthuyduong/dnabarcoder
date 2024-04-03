#!/usr/bin/env python
# FILE: removesequences.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020
import sys
import numpy as np
import os, argparse
from Bio import SeqIO
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='removesequences.py', 
							   usage="%(prog)s [options] -i fastafile -c classificationfile -t taxa -o output -rank species",
							   description='''Script that selects the sequences for of the given taxa. The taxon names in the taxa are separated by ","''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', help='The fasta output file containing the sequences of the given taxa.') #optional
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-t','--taxa', default="", help='the taxa for the selection, separated by ",". If taxa="", all the sequences with the ids given in the classification file will be removed from the input file.')
parser.add_argument('-rank','--classificationrank', default="", help='If the value of the classification rank is not empty, the sequences will be removed.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
taxa=args.taxa
output=args.out
classificationrank=args.classificationrank

#fastafilename=sys.argv[1]
#taxa=sys.argv[2]
#classificationfilename = sys.argv[3]
#output=sys.argv[4]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def SelectSeqIds(classificationfilename,taxa,classificationpos,seqidpos):
	if not os.path.exists(classificationfilename):
		return [],[]
	taxalist=[]
	classnames=[]
	if "," in taxa:
		taxalist=taxa.split(",")
	elif taxa!="" and taxa!="unidentified":
		taxalist.append(taxa)
	classificationfile=open(classificationfilename,errors='ignore')
	seqids=[]
	for line in classificationfile:
		elements=line.rstrip().split("\t")
		seqid = elements[seqidpos].replace(">","").rstrip()
		classname = ""
		if classificationpos >= 0 and classificationpos < len(elements):
			classname = elements[classificationpos]	
		if taxa!="":
			if classname == "" or classname == "unidentified":
				continue
			for taxonname in taxalist:
				if taxonname!="" and taxonname!="unidentified":
					if taxonname in elements:
						seqids.append(seqid)				
						classnames.append(classname)
		else:
			seqids.append(seqid)
			classnames.append(classname)
	return seqids,classnames

def SaveNewClassification(seqids,classificationfilename,newclassificationfilename):
	newclassificationfile=open(newclassificationfilename,"w")
	classificationfile=open(classificationfilename,errors='ignore')
	header=next(classificationfile)
	newclassificationfile.write(header)
	for line in classificationfile:
		seqid=line.split("\t")[0].rstrip()
		if seqid in seqids:
			newclassificationfile.write(line)
	newclassificationfile.close()
	classificationfile.close()		

def GetPosition(classificationfilename,rank):
	pos=-1
	classificationfile=open(classificationfilename,errors='ignore')
	header=classificationfile.readline()
	header=header.rstrip()
	classificationfile.close()
	texts=header.split("\t")
	isError=False
	if rank in texts:
		pos=texts.index(rank)
	else:
		print("The rank " + rank + " is not given in the classification." )
		isError=True
	seqidpos = -1
	i = 0
	for text in texts:
		if text.lower() == args.idcolumnname.lower():
			seqidpos = i
		i = i + 1
	if seqidpos == -1:
		print("Please specify the sequence id columnname by using -idcolumnname.")
		isError = True
	return seqidpos,pos,isError

#####main###
classificationpos=-1
newclassificationfilename=GetBase(output) + ".classification"
if classificationfilename!="":
	seqidpos,classificationpos,isError=GetPosition(classificationfilename,classificationrank)
seqids,classnames=SelectSeqIds(classificationfilename,taxa,classificationpos,seqidpos)
seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
selectedrecords=[]
newseqids=list(set(seqrecords.keys())-set(seqids))
for seqid in newseqids:
	try:
		seqrecord=seqrecords[seqid]
		selectedrecords.append(seqrecord)
	except KeyError:
		continue
#save to file:
SeqIO.write(selectedrecords,output,"fasta")
print("The selected sequences are saved in " + output + ".")
#SaveNewClassification(newseqids, classificationfilename, newclassificationfilename)
#print("The classification of remained sequences are saved in " + newclassificationfilename + ".")
