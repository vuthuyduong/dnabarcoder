#!/usr/bin/env python
# FILE: overview.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020

import sys
import os, argparse
from Bio import SeqIO

parser=argparse.ArgumentParser(prog='updateTaxonomicClassification.py',  
							   usage="%(prog)s [options] -i classificationfilename -t taxon -u newtaxon -o outputname",
							   description='''The script that filter the classification for the sequences of the given fasta file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='The taxonomic classification file.')
parser.add_argument('-t','--taxon', required=True, help='The current taxon.')
parser.add_argument('-u','--newtaxon', required=True, help='The updated taxon.')
parser.add_argument('-o','--out', help='The output filename. If this is not given then the original file will be updateted.')

args=parser.parse_args()
classificationfilename= args.input
taxon= args.taxon
newtaxon= args.newtaxon
outputfilename=args.out

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]  

if outputfilename==None:
	outputfilename=GetBase(classificationfilename) + ".taxaupdated.out"
classificationfile=open(classificationfilename)
header=next(classificationfile)
outputfile=open(outputfilename,"w")
outputfile.write(header)
for line in classificationfile:
	texts=line.split("\t")
	newline=""
	for text in texts:
		if text==taxon:
			text=newtaxon
		newline=newline + text + "\t"	
	newline=newline[:-1]
	outputfile.write(newline)
classificationfile.close()
outputfile.close()
if args.out==None:
	os.system("cp " + outputfilename + " " + classificationfilename)
	os.system("rm " + outputfilename)
	outputfilename=classificationfilename
print("The taxa is updated in  file " + outputfilename  + ".")
