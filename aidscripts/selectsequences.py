#!/usr/bin/env python
# -*- coding: utf-8 -*-
# FILE: selectsequences.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020
import os, argparse
from Bio import SeqIO
import multiprocessing
import random

nproc = multiprocessing.cpu_count()
# from keras.utils import np_utils

parser = argparse.ArgumentParser(prog='selectsequences.py',
                                 usage="%(prog)s [options] -i fastafile -c classificationfile -rank species -t taxa -o output",
                                 description='''Script that selects the sequences for of the given taxa. The taxon names in the taxa are separated by ","''',
                                 epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
                                 )

parser.add_argument('-i', '--input', required=True, help='the fasta file')
parser.add_argument('-o', '--out', required=True,
                    help='The fasta output file containing the sequences of the given taxa.')  # optional
parser.add_argument('-c', '--classification', default="", help='the classification file in tab. format.')
parser.add_argument('-t', '--taxa', default="", help='the taxa for the selection, separated by ","')
parser.add_argument('-rank', '--classificationrank', default="", help='the classification rank for the selection.')
parser.add_argument('-maxseqnopergroup','--maxseqnopergroup', type=int, default=5000, help='The maximum number of the sequences for each group at the given taxonomic rank.')
parser.add_argument('-unique', '--unique', default="no", help='select only unique sequences if yes otherwise no.')
parser.add_argument('-l', '--length', type=int, default=0, help='the required minimum length.')
parser.add_argument('-idcolumnname', '--idcolumnname', default="ID",
                    help='the column name of sequence id in the classification file.')

args = parser.parse_args()
fastafilename = args.input
classificationfilename = args.classification
taxa = args.taxa
output = args.out
maxseqnopergroup = args.maxseqnopergroup
l = args.length
classificationrank = args.classificationrank


# fastafilename=sys.argv[1]
# taxa=sys.argv[2] #separated by ;
# classificationfilename = sys.argv[3]
# output=sys.argv[4]

def GetBase(filename):
    return filename[:-(len(filename) - filename.rindex("."))]


def GetPosition(classificationfilename, rank):
    pos = -1
    seqidpos = -1
    classificationfile = open(classificationfilename,errors='ignore')
    header = classificationfile.readline()
    header = header.rstrip()
    classificationfile.close()
    texts = header.rstrip().split("\t")
    i = 0
    for text in texts:
        if text.lower() == args.idcolumnname.lower():
            seqidpos = i
        i = i + 1
    if seqidpos == -1:
        print("Please specify the sequence id columnname by using -idcolumnname.")
        isError = True
    isError = False
    if rank in texts:
        pos = texts.index(rank)
    elif rank!="":
        print("The rank " + rank + " is not given in the classification.")
        isError = True
    return seqidpos, pos, isError


def LoadClassification(classificationfilename, taxa, classificationpos, seqidpos):
	if not os.path.exists(classificationfilename):
		return {}
	classnames = {}
	classification = {}
	taxalist = []
	if "," in taxa:
		taxalist = taxa.split(",")
	elif taxa != "" and taxa != "unidentified":
		taxalist.append(taxa)
	classificationfile = open(classificationfilename,encoding='latin1')
	header = next(classificationfile)
	for line in classificationfile:
		elements = line.rstrip().split("\t")
		seqid = elements[seqidpos].rstrip()
		classname = ""
		if (classificationpos >= 0 and classificationpos < len(elements)):
			classname = elements[classificationpos]
		else:
			#classname=line.rstrip()
			classname = seqid
		if taxa != "":
			found=False
			for taxonname in taxalist:
				if taxonname != "" and taxonname != "unidentified":
					for element in elements:
						if taxonname.lower() == element.lower():
							found=True
							break
			if found==True:
				classnames.setdefault(seqid, classname)
				classification.setdefault(seqid, line)
		elif classname !="" and classname !="unidentified":
			classnames.setdefault(seqid, classname)
			classification.setdefault(seqid, line)
		else:
			classification.setdefault(seqid, line)
	return classnames, classification, header


def GetTaxonName(description, rank, taxa):
	found=False
	taxalist=[]
	if "," in taxa:
		taxalist = taxa.split(",")
	elif taxa != "" and taxa != "unidentified":
		taxalist.append(taxa)
	else:
		found=True
	taxonname = ""
	species = ""
	genus = ""
	family = ""
	order = ""
	bioclass = ""
	phylum = ""
	kingdom = ""
	seqid=description
	if " " in description:
		seqid=description.split(" ")[0]
		description = description.split(" ")[1]
	texts = description.split("|")
	for text in texts:
		text = text.rstrip()
		taxa = text.split(";")
		for taxon in taxa:
			if "__" in taxon:
				if taxon.split("__")[1] in taxalist:
					found=True
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
	if rank.lower() == "species":
		taxonname = species
	elif rank.lower() == "genus":
		taxonname = genus
	elif rank.lower() == "family":
		taxonname = family
	elif rank.lower() == "order":
		taxonname = order
	elif rank.lower() == "class":
		taxonname = bioclass
	elif rank.lower() == "phylum":
		taxonname = phylum
	elif rank.lower() == "kingdom":
		taxonname = kingdom
	else:
		taxonname=seqid
		#taxonname=kingdom + "\t" + phylum + "\t" + bioclass + "\t" + order + "\t" + family + "\t" + genus + "\t" + species
	if found==False:
		taxonname=""
	return taxonname

def SelectClassName(seqid, description, rank, taxa, classnames):
	classname = ""
	if classnames == {}:
		classname = GetTaxonName(description, rank,taxa)
	else:
		try:
			classname = classnames[seqid]
		except KeyError:
			pass
	# if taxa != "" and classname != "":
	# 	if "\t" in classname:
	# 		names=list(set(taxalist) & set(classname.split("\t")))
	# 		if len(names) ==0:
	# 			classname=""
	# 	elif not (classname in taxalist):
	# 		classname = ""
	return classname


#####main###
classificationpos = -1
seqidpos = -1
classnames = {}
classification = {}
header = ""
if classificationfilename != "":
    seqidpos, classificationpos, isError = GetPosition(classificationfilename, classificationrank)
    # if isError==False:
    # os.sys.exit()
    classnames, classification, header = LoadClassification(classificationfilename, taxa, classificationpos, seqidpos)
seqrecords = SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
selectedrecords = []
selectedclassnames = {}
newclassificationfilename = ""
if classificationfilename != "":
    if "." in output:
        newclassificationfilename = output[0:output.rindex(".")] + ".classification"
    else:
        newclassificationfilename = output + ".classification"

newclassificationfile=None
if newclassificationfilename != "":
	newclassificationfile=open(newclassificationfilename,"w")	
	newclassificationfile.write(header)
#select unique sequences	
uniquesequences = {}
classificationdict={}
for seqid in seqrecords.keys():
	seqrec = seqrecords[seqid]
	description = seqrec.description
	classname=SelectClassName(seqid,description,classificationrank,taxa,classnames)
	if classname != "" and classname !="unidentified":
		if args.unique == "yes":  # select only unique sequences
			try:
				seqrec = uniquesequences[str(seqrec.seq)]
			except KeyError:
				uniquesequences.setdefault(str(seqrec.seq), seqrec)
				try:
					classificationdict[classname].append(seqrec)
				except KeyError:	
					classificationdict.setdefault(classname,[seqrec])
		else:
			try:
				classificationdict[classname].append(seqrec)
			except KeyError:
				classificationdict.setdefault(classname,[seqrec])
				
#select max seq no per group
selectedrecords=[]
for classname in classificationdict.keys():
	seqrecords=classificationdict[classname]
	indexes=range(len(seqrecords))
	if maxseqnopergroup < len(indexes):
		indexes=random.sample(range(len(seqrecords)), k=maxseqnopergroup)
	for index in indexes:
		seqrec=seqrecords[index]
		selectedrecords.append(seqrec)
		if newclassificationfilename != "":
			newclassificationfile.write(classification[seqrec.id])		

# save to file:
SeqIO.write(selectedrecords, output, "fasta")
#save new classification
if newclassificationfilename != "":
    newclassificationfile.close()
if len(selectedrecords) > 0:
    print("The selected sequences are saved in " + output + ".")
    if newclassificationfilename != "":
        print("The new classification filename is saved in " + newclassificationfilename + ".")
else:
    print("No sequences are selected.")
