#!/usr/bin/env python
# FILE: addclassificationfromsequenceheaders.py
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

parser=argparse.ArgumentParser(prog='splitclassificationfromsequenceheaders.py',
							   usage="%(prog)s [options] -i fastafile -prefix prefix",
							   description='''Script that remove classification from sequence headers. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be clustered.')
parser.add_argument('-prefix','--prefix', help='The prefix of the output fasta and classification filenames.')
parser.add_argument('-rank','--classificationranks', default="", help='the classification ranks to compute distribution, separated by ",".')
parser.add_argument('-sep','--separator', default="|", help='the separator that separates the features in sequence headers.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.prefix + ".classification"
output=args.prefix + ".fasta"
output_withclassification=args.prefix + "_classification.fasta"
sep=args.separator

#fastafilename=sys.argv[1]
#classificationfilename=sys.argv[2]
#output=sys.argv[3]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadClassificationFromDescription(seqrecords):
	classificationdict = {}
	for seqid in seqrecords.keys():
		description = seqrecords[seqid].description
		species = "unidentified"
		genus = "unidentified"
		family = "unidentified"
		order = "unidentified"
		bioclass = "unidentified"
		phylum = "unidentified"
		kingdom = "unidentified"
		sh=""
		if " " in description:
			description = description.split(" ")[1]
		texts=[]
		if args.separator!="":
			texts = description.split(args.separator)
		else:
			texts=[description]
		for text in texts:
			text = text.rstrip()
			if text.startswith("SH"):
				sh=text
			taxa = text.split(";")
			for taxon in taxa:
				if taxon[-2:]=="__" or taxon[-15:].lower()==" incertae sedis" or taxon[-15:].lower()=="_incertae_sedis" or taxon[-3]==" sp" or taxon[-3:]=="_sp" or ("unculture" in taxon):
					continue
				if taxon.startswith("k__"):
					kingdom = taxon.replace("k__", "")
				elif taxon.startswith("p__"):
					phylum = taxon.replace("p__", "")
				elif taxon.startswith("c__"):
					bioclass = taxon.replace("c__", "")
				elif taxon.startswith("o__"):
					order = taxon.replace("o__", "")
				elif taxon.startswith("f__"):
					family = taxon.replace("f__", "")
				elif taxon.startswith("g__"):
					genus = taxon.replace("g__", "")
				elif taxon.startswith("s__") and (" " in taxon.replace("s__", "") or "_" in taxon.replace("s__", "")):
					species = taxon.replace("s__", "")
					species = species.replace("_", " ")
		classification =kingdom + "\t" + phylum + "\t" + bioclass + "\t" + order + "\t" + family + "\t" + genus + "\t" + species + "\t" + sh
		taxonomy = "k__" + kingdom + ";p__" + phylum + ";c__" + bioclass + ";o__" + order + ";f__" + family + ";g__" + genus +  ";s__" + species.replace(" ","_")
		newseqid=seqid
		if args.separator in seqid:
			newseqid=seqid.split(args.separator)[0]
		classificationdict.setdefault(newseqid,[classification,str(seqrecords[seqid].seq),taxonomy])
	return classificationdict

#####main###
#load train seq records
seqrecords = SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
classificationdict=LoadClassificationFromDescription(seqrecords)
classificationfile=open(classificationfilename,"w")
outputfile=open(output,"w")
outputfile_withclassification=open(output_withclassification,"w")
classificationfile.write("id\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tSH\n")
for seqid in classificationdict.keys():
	classificationfile.write(seqid + "\t" + classificationdict[seqid][0] + "\n")
	outputfile.write(">" + seqid + "\n")
	outputfile.write(classificationdict[seqid][1] + "\n")
	outputfile_withclassification.write(">" + seqid + " " +  classificationdict[seqid][2] + "\n")
	outputfile_withclassification.write(classificationdict[seqid][1] + "\n")
classificationfile.close()
outputfile.close()
outputfile_withclassification.close()
print("The new fasta file and tab-deliminated classification files are saved in " + output + ", " + classificationfilename + " and " + output_withclassification + ".")
