#!/usr/bin/env python
# FILE: addtexttosequenceheaders.py, requiring Python 2.7 if classification file contains non-ascii characters
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys, argparse
import numpy as np
import os
from Bio import SeqIO
import json
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='addTextTosequenceheaders.py',  
							   usage="%(prog)s [options] -i fastafile -t text -o output",
							   description='''Script that adds some default text to sequence headers. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be clustered.')
parser.add_argument('-o','--out', help='The fasta filename with sequence headers containing the given text.') 
parser.add_argument('-t','--text', required=True, help='the text to be added to the sequence headers')

args=parser.parse_args()
fastafilename= args.input
text=args.text
output=args.out

#fastafilename=sys.argv[1]
#text=sys.argv[2]
#output=sys.argv[3]

outputfile=open(output,"w")
fastafile=open(fastafilename)

for line in fastafile:
	if line.startswith(">"):
		line=line.rstrip()
		header=line + text + "\n"
		outputfile.write(header)	
	else:		
		outputfile.write(line)
fastafile.close()
outputfile.close()		

