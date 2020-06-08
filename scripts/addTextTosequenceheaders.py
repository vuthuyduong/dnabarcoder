#!/usr/bin/env python
# FILE: addtexttosequenceheaders.py, requiring Python 2.7 if classification file contains non-ascii characters
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import numpy as np
import os
from Bio import SeqIO
import json
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils


fastafilename=sys.argv[1]
text=sys.argv[2]
output=sys.argv[3]

outputfile=open(output,"w")
fastafile=open(fastafilename)

for line in fastafile:
	if line.startswith(">"):
		line=line.rstrip()
		header=line + text + "\n"
		outputfile.write(header)	
	else:		
		outputfile.write(line)
	
outputfile.close()		

