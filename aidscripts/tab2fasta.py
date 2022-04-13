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

parser=argparse.ArgumentParser(prog='tab2fasta.py',  
							   usage="%(prog)s [options] -i tabfile -o output",
							   description='''Script that creates a fasta file from a tab. file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the tab. file.')
parser.add_argument('-o','--out', help='The fasta filename.') 
parser.add_argument('-seqidpos','--seqidpos', required=True, help='The position of seqid in the tab. file.') 
parser.add_argument('-seqpos','--seqpos', type=int, required=True, help='The position of sequences in the tab. file.') 


args=parser.parse_args()
tabfilename= args.input
seqidpos=args.seqidpos
seqpos=args.seqpos
output=args.out

#fastafilename=sys.argv[1]
#classificationfilename=sys.argv[2]
#output=sys.argv[3]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

seqidposlist=[]
if "," in seqidpos:
	texts=seqidpos.split(",")
	for text in texts:
		seqidposlist.append(int(text))
else:
	seqidposlist.append(int(seqidpos))		

outputfile=open(output,"w")
tabfile=open(tabfilename, errors='ignore')
next(tabfile)
seqids=[]
for line in tabfile:
	#line=unicode(line,'latin1')
	texts=line.split("\t")
	seqid=""
	for pos in seqidposlist:
		seqid=texts[pos]
		if seqid!="":
			break
	if seqid =="":
		print(line)
		continue
	seq=texts[seqpos]
	if seq=="" or "fail" in seq:
		print(line)
		continue 
	if seqid in seqids:
		print(line)
		continue
	seqids.append(seqid)
	outputfile.write(">" + seqid + "\n")
	outputfile.write(seq.replace("\n","") + "\n")
outputfile.close()		
tabfile.close()
