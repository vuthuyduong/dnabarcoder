#!/usr/bin/env python
# FILE: getLabelsFromNexusTrees.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys, argparse
if sys.version_info[0] >= 3:
	unicode = str
import numpy as np
import os
from Bio import SeqIO
#from nexus import NexusReader
import phylotreelib as pt
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='getLabelsFromNexusTrees.py',
							   usage="%(prog)s [options] -i nexusfile -o output",
							   description='''Script that gets labels from a nexus tree. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the nexus format file.')
parser.add_argument('-o','--out', help='The fasta filename.') 


args=parser.parse_args()
inputfilename= args.input
output=args.out
#######MAIN#########
#n = NexusReader.from_file(inputfilename)
with pt.Nexustreefile(inputfilename) as treefile:
	mytree = treefile.readtree()
outputfile = open(output,"w")
outputfile.write("ID\n")
for label in mytree.leaves:
	outputfile.write(label + "\n")
outputfile.close()