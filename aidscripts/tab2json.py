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
import json
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='tab2fasta.py',  
							   usage="%(prog)s [options] -i tabfile -o output",
							   description='''Script that creates a fasta file from a tab. file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the tab. file.')
parser.add_argument('-o','--out', help='The fasta filename.') 
#parser.add_argument('-seqidpos','--seqidpos', required=True, help='The position of seqid in the tab. file.') 
#parser.add_argument('-seqpos','--seqpos', type=int, required=True, help='The position of sequences in the tab. file.') 
parser.add_argument('-idcolumnname','--idcolumnname', default="ID", help='the column name of sequence id in the input file.')


args=parser.parse_args()
tabfilename= args.input
#seqidpos=args.seqidpos
#seqpos=args.seqpos
output=args.out
idcolumnname=args.idcolumnname
#fastafilename=sys.argv[1]
#classificationfilename=sys.argv[2]
#output=sys.argv[3]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]
	

outputfile=open(output,"w")
tabfile=open(tabfilename, errors='ignore')
header=next(tabfile)
subkeys=header.split("\t")

p_id=0
i=0
for key in subkeys:
	if idcolumnname.lower()==key.rstrip().lower():
		p_id=i
	i=i+1	
	
mydict={}	
for line in tabfile:
	line=line.replace("\"","")
	#line=unicode(line,'latin1')
	texts=line.split("\t")
	key=texts[p_id].rstrip()
	if key=="":
		continue
	mydict.setdefault(key,{})
	i=0
	for text in texts:
		if i!=p_id:
			mydict[key].setdefault(subkeys[i].rstrip(),text.rstrip())
		i=i+1
#save mydict
with open(output, 'w') as f:
	json.dump(mydict,f,indent=4)
	
tabfile.close()
