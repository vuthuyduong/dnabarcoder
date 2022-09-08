#!/usr/bin/env python
# FILE: splitFastaFile.py
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

parser=argparse.ArgumentParser(prog='splitFastaFile.py',
							   usage="%(prog)s [options] -i fastafile -n",
							   description='''Script that split a fastafile to n sub fastafiles. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be clustered.')
parser.add_argument('-n','--number', type=int, help='The number of sub fasta files.')

args=parser.parse_args()
fastafilename= args.input
n=args.number

seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
seqnumber = int(len(seqrecords)/n)
i=0
k=0
subseqrecords=[]
for seqrecord in seqrecords:
	if i==seqnumber:
		if len(subseqrecords) >0:
			#write to file
			k=k+1
			subfilename=fastafilename[:-(len(fastafilename)-fastafilename.rindex("."))] + "_" + str(k) + ".fasta" 
			SeqIO.write(subseqrecords,subfilename,"fasta")
		subseqrecords=[]
		i=0
	else:
		i=i+1
		subseqrecords.append(seqrecord)
#save the last file		
if len(subseqrecords) >0:
	#write to file
	k=k+1
	subfilename=fastafilename[:-(len(fastafilename)-fastafilename.rindex("."))] + "_" + str(k) + ".fasta" 
	SeqIO.write(subseqrecords,subfilename,"fasta")		
print("The sub fasta files are saved in " + fastafilename[:-(len(fastafilename)-fastafilename.rindex("."))] + "_" + str(1) + ".fasta" + ",..., " + fastafilename[:-(len(fastafilename)-fastafilename.rindex("."))] + "_" + str(k) + ".fasta" )
		
		