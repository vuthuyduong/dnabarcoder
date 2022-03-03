#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 17 13:34:05 2021

@author: duong
"""
import sys
import os, argparse
from Bio import SeqIO

parser=argparse.ArgumentParser(prog='cleanClassification.py',  
							   usage="%(prog)s [options] -i classificationfilename -o outputname",
							   description='''The script that clean classification adding missing classification to the sequences. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )


parser.add_argument('-i','--input', required=True, help='The taxonomic classification file.')
parser.add_argument('-o','--out',  help='The output folder.')

args=parser.parse_args()
classificationfilename= args.input
output=args.out

def GetTaxonomicClassification(classificationpos,header,texts):
	classification=""
	p_s=len(texts)
	p_g=len(texts)
	p_f=len(texts)
	p_o=len(texts)
	p_c=len(texts)
	p_p=len(texts)
	p_k=len(texts)
	i=0
	for text in header.split("\t"):
		text=text.rstrip()
		if text.lower()=="species":
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
		i=i+1 
	species=""
	genus=""
	family=""
	order="" 
	bioclass=""
	phylum=""
	kingdom=""
	if p_s< len(texts):
		species=texts[p_s].rstrip()
	if p_g< len(texts):
		genus=texts[p_g].rstrip()
	if p_f< len(texts):
		family=texts[p_f].rstrip()
	if p_o< len(texts):
		order=texts[p_o].rstrip()
	if p_c< len(texts):
		bioclass=texts[p_c].rstrip()
	if p_p< len(texts):
		phylum=texts[p_p].rstrip()
	if p_k< len(texts):
		kingdom=texts[p_k].rstrip()
	if classificationpos==0:
		classificationpos=p_s
	level=7	
	if classificationpos==p_s:	
		level=0
	elif classificationpos==p_g:	
		level=1
	elif classificationpos==p_f:	
		level=2
	elif classificationpos==p_o:	
		level=3
	elif classificationpos==p_c:	
		level=4
	elif classificationpos==p_p:	
		level=5
	elif classificationpos==p_k:	
		level=6
	taxonname=""
	if classificationpos <len(texts):
		taxonname=texts[classificationpos]		
	if level <7 and kingdom!="":
		taxonname=kingdom
		rank="kingdom"
		classification=kingdom
	if level <6 and phylum!="":
		taxonname=phylum
		rank="phylum"
		classification=kingdom +";"+phylum 
	if level <5 and bioclass!="":
		taxonname=bioclass
		rank="class"
		classification=kingdom +";"+phylum +";"+bioclass
	if level <4 and order!="":
		taxonname=order
		rank="order"
		classification=kingdom +";"+phylum +";"+bioclass +";"+order
	if level <3 and family!="":
		taxonname=family
		rank="family"
		classification=kingdom +";"+phylum +";"+bioclass +";"+order+";"+family 
	if level <2 and genus!="":
		taxonname=genus
		rank="genus"
		classification=kingdom +";"+phylum +";"+bioclass +";"+order+";"+family + ";"+ genus
	if level <1 and species!="":
		taxonname=species
		rank="species"
		classification=kingdom +";"+phylum +";"+bioclass +";"+order+";"+family + ";"+ genus+";"+species 	
	return classification,taxonname,rank	

def LoadClassification(classificationfilename,header,pos):
	classificationdict={}
	classificationfile= open(classificationfilename)
	next(classificationfile)
	for line in classificationfile:
		texts=line.split("\t")
		if len(texts) < len(header.split("\t")):
			continue
		#only consider the good format
		if pos < len(texts):
			classname=texts[pos]
			if classname!="":
				classification,classname,rank=GetTaxonomicClassification(pos,header,texts)
				if classname in classificationdict.keys():
					if len(classificationdict[classname].split("\t")) < len(classification.split("\t")):
						print(line)
						print(classificationdict[classname])
						print(classification)
						classificationdict[classname]=classification
				else:
					classificationdict.setdefault(classname,classification)
	classificationfile.close()	
	return classificationdict



##############MAIN############################
classificationfile=open(classificationfilename)
headers=next(classificationfile)
classificationfile.close()
texts=headers.rstrip().split("\t")
p_s=len(texts)
p_g=len(texts)
p_f=len(texts)
p_o=len(texts)
p_c=len(texts)
p_p=len(texts)
p_k=len(texts)
i=0
for text in headers.split("\t"):
	text=text.rstrip()
	if text.lower()=="species":
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
	i=i+1 
	
speciesdict=LoadClassification(classificationfilename,headers,p_s)
genusdict=LoadClassification(classificationfilename,headers,p_g)
familydict=LoadClassification(classificationfilename,headers,p_f)
orderdict=LoadClassification(classificationfilename,headers,p_o)
classdict=LoadClassification(classificationfilename,headers,p_c)
phylumdict=LoadClassification(classificationfilename,headers,p_p)
kingdomdict=LoadClassification(classificationfilename,headers,p_k)

classificationfile=open(classificationfilename)
next(classificationfile)
outputfile=open(output,"w")	
outputfile.write(headers)
for line in classificationfile:
	texts=line.split("\t")
	species=""
	genus=""
	family=""
	order=""
	bioclass=""
	phylum=""
	kingdom=""
	for i in range(0,len(texts)):
		text=texts[len(texts)-i-1].rstrip()
		if text in speciesdict.keys():
			classificationtext=speciesdict[text]
			classification=classificationtext.split(";")
			if species=="":
				species = text
			if genus=="":	
				genus=classification[5]
			if family=="":	
				family=classification[4]
			if order=="":	
				order=classification[3]
			if bioclass=="":	
				bioclass=classification[2]
			if phylum=="":
				phylum=classification[1]
			if kingdom=="":
				kingdom=classification[0]
		elif text in genusdict.keys():
			classificationtext=genusdict[text]
			classification=classificationtext.split(";")
			if genus=="":	
				genus=text
			if family=="":	
				family=classification[4]
			if order=="":	
				order=classification[3]
			if bioclass=="":	
				bioclass=classification[2]
			if phylum=="":
				phylum=classification[1]
			if kingdom=="":
				kingdom=classification[0]
		elif text in familydict.keys():
			classificationtext=familydict[text]
			classification=classificationtext.split(";")
			if family=="":	
				family=text
			if order=="":	
				order=classification[3]
			if bioclass=="":	
				bioclass=classification[2]
			if phylum=="":
				phylum=classification[1]
			if kingdom=="":
				kingdom=classification[0]
		elif text in orderdict.keys():	
			classificationtext=orderdict[text]
			classification=classificationtext.split(";")
			if order=="":	
				order=text
			if bioclass=="":	
				bioclass=classification[2]
			if phylum=="":
				phylum=classification[1]
			if kingdom=="":
				kingdom=classification[0]
		elif text in classdict.keys():	
			classificationtext=classdict[text]
			classification=classificationtext.split(";")
			if bioclass=="":	
				bioclass=text
			if phylum=="":
				phylum=classification[1]
			if kingdom=="":
				kingdom=classification[0]
		elif text in phylumdict.keys():	
			classificationtext=phylumdict[text]
			classification=classificationtext.split(";")
			if phylum=="":
				phylum=text
			if kingdom=="":
				kingdom=classification[0]		
	i=0		
	newclassification=""
	for i in range(0,len(headers.split("\t"))):
		text=""
		if i <len(texts):
			text=texts[i].rstrip()
		if i==p_s:
			text=species
		if i==p_g:
			text=genus
		if i==p_f:
			text=family
		if i==p_o:
			text=order
		if i==p_c:
			text=bioclass
		if i==p_p:
			text=phylum
		if i==0:
			newclassification=text
		else:
			newclassification=newclassification + "\t" + text
		i=i+1	
	outputfile.write(newclassification+"\n")
outputfile.close()