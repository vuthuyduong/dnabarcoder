#!/usr/bin/env python
# -*- coding: utf-8 -*-
# FILE: updateClassification.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020

import sys
import os, argparse
from Bio import SeqIO

parser=argparse.ArgumentParser(prog='updateClassification.py',
							   usage="%(prog)s [options] -i classificationfilename -t taxon -u newtaxon -o outputname",
							   description='''The script that updates new classifications for the sequences of the given classification file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='The current taxonomic classification file.')
parser.add_argument('-t','--taxon', default="", help='The current taxon.')
parser.add_argument('-nt','--newtaxon', help='The updated taxon.')
parser.add_argument('-nf','--newclassificationfile', help='The file of updated taxa.')
parser.add_argument('-o','--out', help='The output filename.')

args=parser.parse_args()
classificationfilename= args.input
taxon= args.taxon
newtaxon= args.newtaxon
newclassificationfilename= args.newclassificationfile
outputfilename=args.out

def LoadClassification(newclassificationfilename):
	classification={}
	classificationfile=open(newclassificationfilename)
	seqid="NA"
	strainid="NA"
	species = "NA"
	genus = "NA"
	family = "NA"
	order = "NA"
	bioclass = "NA"
	phylum = "NA"
	kingdom = "NA"
	header = next(classificationfile)
	for line in classificationfile:
		seqid, strainid, species, genus, family, order, bioclass, phylum, kingdom=GetTaxa(header,line)
		if seqid!="NA" and seqid!="":
			classification.setdefault(seqid,"k__" + kingdom + ";p__" + phylum + ";c__" + bioclass + ";o__" + order + ";f__" + family + ";g__" + genus + ";s__" + species)
		if strainid!="NA" and strainid!="":
			classification.setdefault(seqid,"k__" + kingdom + ";p__" + phylum + ";c__" + bioclass + ";o__" + order + ";f__" + family + ";g__" + genus + ";s__" + species)
		if species!="NA" and species!="":
			classification.setdefault(species,"k__" + kingdom + ";p__" + phylum + ";c__" + bioclass + ";o__" + order + ";f__" + family + ";g__" + genus + ";s__" + species)
		if genus!="NA" and genus!="":
			classification.setdefault(genus,"k__" + kingdom + ";p__" + phylum + ";c__" + bioclass + ";o__" + order + ";f__" + family + ";g__" + genus )
		if family!="NA" and family!="":
			classification.setdefault(family,"k__" + kingdom + ";p__" + phylum + ";c__" + bioclass + ";o__" + order + ";f__" + family)
		if order!="NA" and order!="":
			classification.setdefault(order,"k__" + kingdom + ";p__" + phylum + ";c__" + bioclass + ";o__" + order)
		if bioclass!="NA" and bioclass!="":
			classification.setdefault(bioclass,"k__" + kingdom + ";p__" + phylum + ";c__" + bioclass )
		if phylum!="NA" and phylum!="":
			classification.setdefault(phylum,"k__" + kingdom + ";p__" + phylum)
		if kingdom!="NA" and kingdom!="":
			classification.setdefault(kingdom,"k__" + kingdom)
	classificationfile.close()
	return classification

def GetClassification(key,classification):
	strainid="NA"
	species="NA"
	genus="NA"
	family="NA"
	order="NA"
	bioclass="NA"
	phylum="NA"
	kingdom="NA"
	newtaxa=classification[key]
	taxonlist=newtaxa.split(";")
	for taxon in taxonlist:
		if taxon.startswith("s__"):
			species=taxon.replace("s__","")
		elif taxon.startswith("g__"):
			genus=taxon.replace("g__","")
		elif taxon.startswith("f__"):
			family=taxon.replace("f__","")
		elif taxon.startswith("o__"):
			order=taxon.replace("o__","")
		elif taxon.startswith("c__"):
			bioclass=taxon.replace("c__","")
		elif taxon.startswith("p__"):
			phylum=taxon.replace("p__","")
		elif taxon.startswith("k__"):
			kingdom=taxon.replace("k__","")
	return species,genus,family,order,bioclass,phylum,kingdom

def GetTaxa(header,line):
	seqid = "NA"
	strainid = "NA"
	species = "NA"
	genus = "NA"
	family = "NA"
	order = "NA"
	bioclass = "NA"
	phylum = "NA"
	kingdom = "NA"
	texts = header.split("\t")
	p_seqid=-1
	p_strainid=-1
	p_s=-1
	p_g=-1
	p_f=-1
	p_o=-1
	p_c=-1
	p_p=-1
	p_k=-1
	i = 0
	for text in texts:
		text = text.rstrip()
		if text.lower() == "id" or text.lower() == "seqid" or text.lower() == "seq id" or text.lower() == "sequence id" or text.lower() == "sequenceid":
			p_seqid = i
		elif text.lower() == "strain id" or text.lower() == "strainid":
			p_strainid = i
		elif text.lower() == "species":
			p_s = i
		elif text.lower() == "genus":
			p_g = i
		elif text.lower() == "family":
			p_f = i
		elif text.lower() == "order":
			p_o = i
		elif text.lower() == "class":
			p_c = i
		elif text.lower() == "phylum":
			p_p = i
		elif text.lower() == "kingdom":
			p_k = i
		i = i + 1
	texts = line.split("\t")
	if p_seqid > -1:
		seqid = texts[p_seqid].rstrip()
	if p_strainid > -1:
		strainid = texts[p_strainid].rstrip()
	if p_s > -1:
		species = texts[p_s].rstrip()
	if p_g > -1:
		genus = texts[p_g].rstrip()
	if p_f > -1:
		family = texts[p_f].rstrip()
	if p_o > -1:
		order = texts[p_o].rstrip()
	if p_c > -1:
		bioclass = texts[p_c].rstrip()
	if p_p > -1:
		phylum = texts[p_p].rstrip()
	if p_k > -1:
		kingdom = texts[p_k].rstrip()
	return seqid,strainid,species,genus,family,order,bioclass,phylum,kingdom

def UpdateTaxa(header,line,newclassification):
	seqid,strainid,species,genus,family,order,bioclass,phylum,kingdom = GetTaxa(header,line)
	newspecies = "NA"
	newgenus = "NA"
	newfamily = "NA"
	neworder = "NA"
	newclass = "NA"
	newphylum = "NA"
	newkingdom = "NA"
	if seqid!="NA" and seqid!="" and seqid in newclassification.keys():
		newspecies, newgenus, newfamily, neworder, newclass, newphylum, newkingdom = GetClassification(seqid,newclassification)
	elif strainid!="NA" and strainid!="" and strainid in newclassification.keys():
		newspecies, newgenus, newfamily, neworder, newclass, newphylum, newkingdom = GetClassification(strainid,newclassification)
	elif species!="NA" and species!="" and species in newclassification.keys():
		newstrainid,newspecies, newgenus, newfamily, neworder, newclass, newphylum, newkingdom = GetClassification(species,newclassification)
	elif genus!="NA" and genus!="" and genus in newclassification.keys():
		newspecies, newgenus, newfamily, neworder, newclass, newphylum, newkingdom = GetClassification(genus,newclassification)
	elif family!="NA" and family!="" and family in newclassification.keys():
		newspecies, newgenus, newfamily, neworder, newclass, newphylum, newkingdom = GetClassification(family,newclassification)
	elif order!="NA" and order!="" and order in newclassification.keys():
		newspecies, newgenus, newfamily, neworder, newclass, newphylum, newkingdom = GetClassification(order,newclassification)
	elif bioclass!="NA" and bioclass!="" and bioclass in newclassification.keys():
		newspecies, newgenus, newfamily, neworder, newclass, newphylum, newkingdom = GetClassification(bioclass,newclassification)
	elif phylum!="NA" and phylum!="" and phylum in newclassification.keys():
		newspecies, newgenus, newfamily, neworder, newclass, newphylum, newkingdom = GetClassification(phylum,newclassification)
	elif kingdom!="NA" and kingdom!="" and kingdom in newclassification.keys():
		newspecies, newgenus, newfamily, neworder, newclass, newphylum, newkingdom = GetClassification(kingdom,newclassification)
	newline = ""
	texts=line.split("\t")
	for text in texts:
		if text.lower()==species.lower() and newspecies!="NA":
			text=newspecies
		elif text.lower()==genus.lower() and newgenus!="NA":
			text=newgenus
		elif text.lower()==family.lower() and newfamily!="NA":
			text=newfamily
		elif text.lower()==order.lower() and neworder!="NA":
			text=neworder
		elif text.lower()==bioclass.lower() and newclass!="NA":
			text=newclass
		elif text.lower()==phylum.lower() and newphylum!="NA":
			text=newphylum
		elif text.lower()==kingdom.lower() and newkingdom!="NA":
			text=newkingdom
		newline=newline + text + "\t"
	newline=newline[:-1]
	return newline

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]  
####################################MAIN#############################################
if outputfilename==None:
	outputfilename=GetBase(classificationfilename) + ".updated"
newclassification=LoadClassification(newclassificationfilename)
classificationfile=open(classificationfilename, errors='ignore')
header=next(classificationfile)
outputfile=open(outputfilename,"w")
outputfile.write(header)
for line in classificationfile:
	newline=UpdateTaxa(header,line,newclassification)
	texts=newline.split("\t")
	if taxon!="":
		newline = ""
		for text in texts:
			if text==taxon:
				text=newtaxon
			newline=newline + text + "\t"
		newline=newline[:-1]
	outputfile.write(newline)
classificationfile.close()
outputfile.close()

print("The taxa are updated in  file " + outputfilename  + ".")
