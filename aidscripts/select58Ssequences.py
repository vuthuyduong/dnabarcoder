#!/usr/bin/env python
# FILE: select58Ssequences.py, requiring Python 2.7 if classification file contains non-ascii characters
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
import sys, argparse
import numpy as np
import os
from Bio import SeqIO
import json
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='select58Ssequences.py',  
							   usage="%(prog)s [options] -i fastafile  -t ITSxpositionfile -o output",
							   description='''Script that select the 5.8S regions from the sequences. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file.')
parser.add_argument('-t','--ITSxpositions', required=True, help='the ITSx positions file obtained from the ITSx software.')
parser.add_argument('-o','--out', help='The output fasta file with 5.8S sequences.')

args=parser.parse_args()
fastafilename= args.input
ITSxpositionfilename=args.ITSxpositions
output=args.out

ITSxpositionfile=open(ITSxpositionfilename)
positiondict={}
for line in ITSxpositionfile:
	texts=line.split("\t")
	seqid=texts[0]
	if not "5.8S: " in line:
		continue
	texts=line.split("5.8S: ")
	text=texts[1]
	text=text.split("\t")[0]
	if not "-" in text:
		continue
	startpos=int(text.split("-")[0])
	endpos=int(text.split("-")[1])
	positiondict.setdefault(seqid,[startpos,endpos])
ITSxpositionfile.close()	
seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
newseqrecords=[]
for seqrecord in seqrecords:
	newseqrecord=seqrecord
	seq=seqrecord.seq
	if not seqrecord.id in positiondict.keys():
		continue
	position=positiondict[seqid]
	newseq=seq[position[0]-1:position[1]]
	newseqrecord.seq=newseq
	newseqrecords.append(newseqrecord)
#write
SeqIO.write(newseqrecords,output,"fasta")	



