#!/usr/bin/env python
# FILE: abundancetable2fasta.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 May 2021
import sys, argparse
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

parser=argparse.ArgumentParser(prog='abundancetable2fasta.py',  
							   usage="%(prog)s [options] -i asvfilename -o outputprefix",
							   description='''Script that convert asv (amplicon sequence variants) to a fasta file and an abundance table file.''',
							   epilog="""Written by Duong Vu d.vu@wi.knaw.nl""",
   )

parser.add_argument('-i','--input', required=True, help='the asv output of data.')
parser.add_argument('-o','--out', help='The folder for the results.') 

args=parser.parse_args()
inputfilename= args.input
outputprefix= args.out

####################MAIN####################
fastafilename=outputprefix + ".fasta"

asvtable =outputprefix + ".table"
tmpfile=open(inputfilename)
firstline=next(tmpfile)
if firstline.startswith("\t"):
	firstline=firstline[1:]
asvs=firstline.rstrip().split("\t")
#write asvs to a fastafile
fastafile=open(fastafilename,"w")
i=0
table={}
for asv in asvs:
	i=i+1
	asvid="asv_" + str(i)
	fastafile.write(">" + asvid + "\n")
	fastafile.write(asv + "\n")
	table.setdefault(asvid,{})
fastafile.close()
sampleids=[]
for line in tmpfile:
	texts=line.rstrip().split("\t")
	sampleid=texts[0]
# 	if "_" in sampleid:
# 		sampleid=sampleid.split("_")[0]
	sampleids.append(sampleid)
	for i in range(1,len(texts)):
		abundance=texts[i]
		asvid="asv_" + str(i)
		if int(abundance) > 0:
			table[asvid].setdefault(sampleid,abundance)
tmpfile.close()
#save abundance table	
tablefile=open(asvtable,"w")
header=""
for sampleid in sampleids:
	header=header + "\t" + sampleid
header=header+"\n"
tablefile.write(header)
i=0
for asv in asvs:
	i=i+1
	asvid="asv_" + str(i)
	line=asvid
	for sampleid in sampleids:
		if sampleid in table[asvid].keys():
			line=line+ "\t" + table[asvid][sampleid]
		else:
			line=line+ "\t0"
	line=line+"\n"
	tablefile.write(line)
tablefile.close()
print("The asvs are saved in file " + fastafilename + ".")
print("The abundance table is saved in file " + asvtable + ".")
