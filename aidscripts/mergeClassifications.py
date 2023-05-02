#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 11:56:08 2021

@author: duong
"""
import sys
import re
if sys.version_info[0] >= 3:
	unicode = str
import json
import os, argparse
from Bio import SeqIO
import random
import multiprocessing
parser=argparse.ArgumentParser(prog='mergeClassification.py',
							   usage="%(prog)s [options] -i classificationfiles -o output",
							   description='''Script that compares two classifications.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--classificationfilenames', required=True, help='the classification filenames (the outputs of dnabarcoder.py classify, for example), separated by commas')
parser.add_argument('-o','--out', required=True, help='The output.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')
parser.add_argument('-scorecolumnname','--scorecolumnname',default="score", help='the column name of similarity scores in the classification file.')

args=parser.parse_args()
output=args.out

def has_numbers(inputString):
    return bool(re.search(r'\d', inputString))

def SelectBetterClassification(taxonomy1,taxonomy2):
	taxonomy={}
	rank=""
	if taxonomy1["species"]!="unidentified" or taxonomy2["species"]!="unidentified":
		rank="species"
		if taxonomy1["species"]!="unidentified" and taxonomy2["species"]!="unidentified":
			if taxonomy1["score"] >= taxonomy2["score"]:
				taxonomy =taxonomy1
			else:
				taxonomy = taxonomy2
		elif taxonomy1["species"]!="unidentified":
			taxonomy = taxonomy1
		else:
			taxonomy = taxonomy2
	elif taxonomy1["genus"]!="unidentified" or taxonomy2["genus"]!="unidentified":
		rank = "genus"
		if taxonomy1["genus"]!="unidentified" and taxonomy2["genus"]!="unidentified":
			if taxonomy1["score"] >= taxonomy2["score"]:
				taxonomy =taxonomy1
			else:
				taxonomy = taxonomy2
		elif taxonomy1["genus"]!="unidentified":
			taxonomy = taxonomy1
		else:
			taxonomy = taxonomy2
	elif taxonomy1["family"]!="unidentified" or taxonomy2["family"]!="unidentified":
		rank = "family"
		if taxonomy1["family"]!="unidentified" and taxonomy2["family"]!="unidentified":
			if taxonomy1["score"] >= taxonomy2["score"]:
				taxonomy =taxonomy1
			else:
				taxonomy = taxonomy2
		elif taxonomy1["family"]!="unidentified":
			taxonomy = taxonomy1
		else:
			taxonomy = taxonomy2
	elif taxonomy1["order"]!="unidentified" or taxonomy2["order"]!="unidentified":
		rank = "order"
		if taxonomy1["order"]!="unidentified" and taxonomy2["order"]!="unidentified":
			if taxonomy1["score"] >= taxonomy2["score"]:
				taxonomy =taxonomy1
			else:
				taxonomy = taxonomy2
		elif taxonomy1["order"]!="unidentified":
			taxonomy = taxonomy1
		else:
			taxonomy = taxonomy2
	elif taxonomy1["class"]!="unidentified" or taxonomy2["class"]!="unidentified":
		rank = "class"
		if taxonomy1["class"]!="unidentified" and taxonomy2["class"]!="unidentified":
			if taxonomy1["score"] >= taxonomy2["score"]:
				taxonomy =taxonomy1
			else:
				taxonomy = taxonomy2
		elif taxonomy1["class"]!="unidentified":
			taxonomy = taxonomy1
		else:
			taxonomy = taxonomy2
	elif taxonomy1["phylum"]!="unidentified" or taxonomy2["phylum"]!="unidentified":
		rank = "phylum"
		if taxonomy1["phylum"]!="unidentified" and taxonomy2["phylum"]!="unidentified":
			if taxonomy1["score"] >= taxonomy2["score"]:
				taxonomy =taxonomy1
			else:
				taxonomy = taxonomy2
		elif taxonomy1["phylum"]!="unidentified":
			taxonomy = taxonomy1
		else:
			taxonomy = taxonomy2
	elif taxonomy1["kingdom"]!="unidentified" or taxonomy2["kingdom"]!="unidentified":
		rank = "kingdom"
		if taxonomy1["kingdom"]!="unidentified" and taxonomy2["kingdom"]!="unidentified":
			if taxonomy1["score"] >= taxonomy2["score"]:
				taxonomy =taxonomy1
			else:
				taxonomy = taxonomy2
		elif taxonomy1["kingdom"]!="unidentified":
			taxonomy = taxonomy1
		else:
			taxonomy = taxonomy2
	else:
		taxonomy=taxonomy2
	taxonomy["rank"]=rank
	return taxonomy

def LoadClassification(classificationdict,classificationfilename):
	error=False
	classificationfile=open(classificationfilename)
	header=next(classificationfile)
	p_s=-1
	p_g=-1
	p_f=-1
	p_o=-1
	p_c=-1
	p_p=-1
	p_k=-1
	p_score=-1
	p_id=-1
	p_cutoff=-1
	p_confidence=-1
	p_refid=-1
	p_rank=-1
	i=0
	for text in header.split("\t"):
		text=text.rstrip()
		if text.lower()==args.idcolumnname.lower():
			p_id=i
		elif "reference id" in text.lower() or "referenceid" in text.lower() or "refid" in text.lower():
			p_refid=i
		elif text.lower()=="cutoff":
			p_cutoff=i
		elif text.lower()=="confidence":
			p_confidence=i
		elif text.lower() == "rank":
			p_rank = i
		elif text.lower()=="species":
			p_s=i
		elif text.lower()=="genus":
			p_g=i
		elif text.lower()=="family":
			p_f=i
		elif text.lower()=="order":
			p_o=i
		elif text.lower()=="class":
			p_c=i
		elif text.lower()=="phylum":
			p_p=i
		elif text.lower()=="kingdom":
			p_k=i
		elif text.lower()== args.scorecolumnname:
			p_score=i
		i=i+1
	if p_id==-1:
		print("Please check the id column name in the classification file " + classificationfilename + ".")
		error=True
	for line in classificationfile:
		texts=line.split("\t")
		seqid=texts[p_id].rstrip()
		score = 0
		cutoff=0
		confidence=0
		refid=""
		rank=""
		if p_score > -1 and p_score < len(texts):
			score = float(texts[p_score].rstrip())
		if p_cutoff > -1 and p_cutoff < len(texts):
			cutoff=-1
			try:
				cutoff = float(texts[p_cutoff].rstrip())
			except:
				pass
		if p_confidence > -1 and p_confidence < len(texts):
			try:
				confidence = float(texts[p_confidence].rstrip())
			except:
				pass
		if p_rank > -1 and p_rank < len(texts):
			rank = texts[p_rank].rstrip()
		if p_refid > -1 and p_refid < len(texts):
			refid = texts[p_refid].rstrip()
		species="unidentified"
		if p_s >-1 and p_s < len(texts):
			species = texts[p_s].rstrip()
			if species=="" or has_numbers(species)==True:
				species="unidentified"
		genus = "unidentified"
		if p_g > -1 and p_g < len(texts):
			genus = texts[p_g].rstrip()
			if genus=="" or has_numbers(genus)==True:
				genus = "unidentified"
		family = "unidentified"
		if p_f > -1 and p_f < len(texts):
			family= texts[p_f].rstrip()
			if family=="" or has_numbers(family)==True:
				family = "unidentified"
		order = "unidentified"
		if p_o > -1 and p_o < len(texts):
			order = texts[p_o].rstrip()
			if order=="" or has_numbers(order)==True:
				order = "unidentified"
		bioclass = "unidentified"
		if p_c > -1 and p_c < len(texts):
			bioclass = texts[p_c].rstrip()
			if bioclass=="" or has_numbers(bioclass)==True:
				bioclass = "unidentified"
		phylum = "unidentified"
		if p_p > -1 and p_p < len(texts):
			phylum = texts[p_p].rstrip()
			if phylum=="" or has_numbers(phylum)==True:
				phylum = "unidentified"
		kingdom="unidentified"
		if p_k > -1 and p_k < len(texts):
			kingdom = texts[p_k].rstrip()
			if kingdom=="" or has_numbers(kingdom)==True:
				kingdom = "unidentified"
		taxonomy = {}
		taxonomy.setdefault("score", score)
		taxonomy.setdefault("cutoff", cutoff)
		taxonomy.setdefault("confidence", confidence)
		taxonomy.setdefault("rank", rank)
		taxonomy.setdefault("referenceid", refid)
		taxonomy.setdefault("species", species)
		taxonomy.setdefault("genus", genus)
		taxonomy.setdefault("family", family)
		taxonomy.setdefault("order", order)
		taxonomy.setdefault("class", bioclass)
		taxonomy.setdefault("phylum", phylum)
		taxonomy.setdefault("kingdom", kingdom)
		if seqid in classificationdict.keys():
			existingtaxonomy=classificationdict[seqid]
			taxonomy=SelectBetterClassification(existingtaxonomy,taxonomy)
			classificationdict[seqid]=taxonomy
		else:
			classificationdict.setdefault(seqid, taxonomy)
	classificationfile.close()
	return error

def SaveClassification(classificationdict,output):
	outputfile=open(output,"w")
	outputfile.write("id\tReferenceID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\trank\tscore\tcutoff\tconfidence\n")
	for id in classificationdict.keys():
		taxonomy=classificationdict[id]
		outputfile.write(id + "\t")
		outputfile.write(taxonomy["referenceid"] + "\t")
		outputfile.write(taxonomy["kingdom"] + "\t")
		outputfile.write(taxonomy["phylum"] + "\t")
		outputfile.write(taxonomy["class"] + "\t")
		outputfile.write(taxonomy["order"] + "\t")
		outputfile.write(taxonomy["family"] + "\t")
		outputfile.write(taxonomy["genus"] + "\t")
		outputfile.write(taxonomy["species"] + "\t")
		outputfile.write(taxonomy["rank"] + "\t")
		outputfile.write(str(taxonomy["score"]) + "\t")
		outputfile.write(str(taxonomy["cutoff"]) + "\t")
		outputfile.write(str(taxonomy["confidence"]) + "\n")
	print("The merged classification is saved in file " + output + ".")
###########MAIN########################
classificationnames=[]
if "," in args.classificationfilenames:
	classificationnames = args.classificationfilenames.split(",")
else:
	classificationnames.append(args.classificationfilenames)
classificationdict={}
#merge classifications
for classificationname in classificationnames:
	error=LoadClassification(classificationdict,classificationname)
#save merged classification
SaveClassification(classificationdict,output)

