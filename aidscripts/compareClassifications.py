#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 11:56:08 2021

@author: duong
"""
import sys
if sys.version_info[0] >= 3:
	unicode = str
import json
import os, argparse
from Bio import SeqIO
import random
import multiprocessing
parser=argparse.ArgumentParser(prog='compareClassification.py',  
							   usage="%(prog)s [options] -c1 cbsclassification1 -c2 classification2 -o output",
							   description='''Script that compares two classifications.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-c1','--classificationfilename1', required=True, help='the classification filename 1')
parser.add_argument('-c2','--classificationfilename2', required=True, help='the classification filename 2')
parser.add_argument('-o','--out', required=True, help='The output.')

args=parser.parse_args()
classificationfilename1= args.classificationfilename1
classificationfilename2= args.classificationfilename2
output=args.out

def LoadClassificationAtPos(classificationfilename,pos,p_score):
	classificationfile=open(classificationfilename)
	next(classificationfile)
	taxa={}
	for line in classificationfile:
		texts=line.rstrip().split("\t")
		seqid=texts[0]
		taxonname=""
		score=0
		if pos <len(texts):
			taxonname=texts[pos]
			taxonname=taxonname.replace("_"," ")
			if ("uncultured" in taxonname) or ("unidentified" in taxonname) or (" sp " in taxonname) :
				taxonname="unidentified"
		if p_score <len(texts):
			score=float(texts[p_score])
		if taxonname=="" or taxonname=="unidentified":
			continue
		taxa.setdefault(seqid,[taxonname,score])
	classificationfile.close()
	return taxa

def LoadClassification(classificationfilename):
	kingdoms={}
	phyla={}
	classes={}
	orders={}
	families={}
	genera={}
	species={}
	classificationfile=open(classificationfilename)
	header=next(classificationfile)
	classificationfile.close()
	p_s=-1
	p_g=-1
	p_f=-1
	p_o=-1
	p_c=-1
	p_p=-1
	p_k=-1
	p_score=-1
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
		elif text.lower()=="score":
			p_score=i	
		i=i+1 
	species=LoadClassificationAtPos(classificationfilename,p_s,p_score)	
	genera=LoadClassificationAtPos(classificationfilename,p_g,p_score)	
	families=LoadClassificationAtPos(classificationfilename,p_f,p_score)	
	orders=LoadClassificationAtPos(classificationfilename,p_o,p_score)	
	classes=LoadClassificationAtPos(classificationfilename,p_c,p_score)		
	phyla=LoadClassificationAtPos(classificationfilename,p_p,p_score)	
	kingdoms=LoadClassificationAtPos(classificationfilename,p_k,p_score)	
	return kingdoms,phyla,classes,orders,families,genera,species

def Compare(rank,classes1,classes2,classification1,classification2,output):
	outputfile=open(output,"a")
	taxafile=open(output + "." + rank,"w")
	number1=len(classes1.keys())
	number2=len(classes2.keys())
	seqids=list(set(classes1.keys()) & set(classes2.keys()))
	number=len(seqids)
	taxafile.write("Sequence ID\tClassification of " + classification1 + "\tScore given in " + classification1 + "\tClassification of" + classification2 + "\tScore given in " + classification2 +"\tCompared" + "\n")
	same=0
	badseqids1=[]
	badseqids2=[]
	seqids1=[]
	seqids2=[]
	onlyseqids1=[]
	onlyseqids2=[]
	notonlyseqids1=[]
	notonlyseqids2=[]
	count=0
	for seqid in classes1.keys():
		taxon1=classes1[seqid][0]
		score1=classes1[seqid][1]
		taxon2=""
		score2=0
		compared=""	
		if seqid in classes2.keys():
			taxon2=classes2[seqid][0]
			score2=classes2[seqid][1]
			if taxon1==taxon2:
				same=same+1
				compared="same"
				seqids1.append(seqid)	
				seqids2.append(seqid)	
				notonlyseqids1.append(seqid)	
				notonlyseqids2.append(seqid)	
			else:
				compared="diff"
				if score1>=score2:
					onlyseqids1.append(seqid)	
					seqids1.append(seqid)	
				if score2>=score1:
					onlyseqids2.append(seqid)	
					seqids2.append(seqid)		
				if score1<score2:
					badseqids1.append(seqid)
					notonlyseqids1.append(seqid)	
				elif score2<score1:
					count=count+1
					badseqids2.append(seqid)
					notonlyseqids2.append(seqid)	
		else:
			onlyseqids1.append(seqid)	
			seqids1.append(seqid)	
		taxafile.write(seqid + "\t" + taxon1 + "\t" +str(score1) + "\t" + taxon2 + "\t" + str(score2) + "\t" + compared + "\n")
	for seqid in classes2.keys():
		if seqid in classes1.keys():
			continue
		onlyseqids2.append(seqid)	
		seqids2.append(seqid)	
		taxon1=""
		score1=0
		taxon2=classes2[seqid][0]
		score2=classes2[seqid][1]
		compared=""	
		taxafile.write(seqid + "\t" + taxon1 + "\t" + str(score1) + "\t" + taxon2 + "\t" + str(score2) +"\t" + compared + "\n")	
	taxafile.close()	
	difference=number-same
	onlyby1=number1-number
	onlyby2=number2-number
	higherscoreby1=count
	lessscoreby1=difference-count
	higherscoreby2=difference-count
	lessscoreby2=count
	outputfile.write(rank + "\t" + str(number1)+ "\t" + str(number2)+ "\t" + str(number) + "\t" + str(same) + "\t" + str(difference) + "\t" + str(higherscoreby1) + "\t" + str(higherscoreby2) + "\t" + str(lessscoreby1) + "\t" + str(lessscoreby2) + "\t" + str(onlyby1) + "\t" + str(onlyby2) + "\n" )
	outputfile.close()
	print(rank + "\t" + str(number1)+ "\t" + str(number2)+ "\t" + str(number) + "\t" + str(same) + "\t" + str(difference) + "\t" + str(higherscoreby1) + "\t" + str(higherscoreby2) + "\t" + str(lessscoreby1) + "\t" + str(lessscoreby2) + "\t" + str(onlyby1) + "\t" + str(onlyby2) + "\n" )
	print("The comparison at the " + rank + " taxonomic level is saved in file " +  output + "." + rank + "." )
	return seqids1,seqids2,onlyseqids1,onlyseqids2
	
def CreateNewClassificationFile(filename,newfilename,skeys,selectedsseqids,gkeys,selectedgseqids,fkeys,selectedfseqids,okeys,selectedoseqids,ckeys,selectedcseqids,pkeys,selectedpseqids,kkeys,selectedkseqids):
	newfile=open(newfilename,"w")
	file=open(filename)
	header=next(file)
	newfile.write(header)
	for line in file:
		seqid=line.split("\t")[0]
		isok=0
		if (seqid in skeys):
			if (seqid in selectedsseqids):
				isok=1
		elif (seqid in gkeys):
			if seqid in selectedgseqids:
				isok=1
		elif (seqid in fkeys):
			if (seqid in selectedfseqids):
				isok=1
		elif (seqid in okeys):
			if (seqid in selectedoseqids):
				isok=1
		elif (seqid in ckeys):
			if (seqid in selectedcseqids):
				isok=1
		elif (seqid in pkeys):
			if (seqid in selectedpseqids):
				isok=1
		elif (seqid in kkeys):
			if (seqid in selectedkseqids):
				isok=1
		if isok==1:
			newfile.write(line)	
	newfile.close()
	file.close()
###

classification1=os.path.basename(classificationfilename1)
classification2=os.path.basename(classificationfilename2)
kingdoms1,phyla1,classes1,orders1,families1,genera1,species1=LoadClassification(classificationfilename1)
kingdoms2,phyla2,classes2,orders2,families2,genera2,species2=LoadClassification(classificationfilename2)
outputfile=open(output,"w")	
outputfile.write("Rank\tSequence number of " + classification1+ "\tSequence number of " + classification2 + "\tNumber of common sequences\tNumber of the same classification\tNumber of different classifications\tHigher scores by first classification\tHigher scores by second classification\tLess scores by first classification\tLess scores by second classification\tOnly by first classification\tOnly by second classification\n")
outputfile.close()
sseqids1,sseqids2,onlysseqids1,onlysseqids2=Compare("species",species1,species2,classification1,classification2,output)
gseqids1,gseqids2,onlygseqids1,onlygseqids2=Compare("genus",genera1,genera2,classification1,classification2,output)
fseqids1,fseqids2,onlyfseqids1,onlyfseqids2=Compare("family",families1,families2,classification1,classification2,output)
oseqids1,oseqids2,onlyoseqids1,onlyoseqids2=Compare("order",orders1,orders2,classification1,classification2,output)
cseqids1,cseqids2,onlycseqids1,onlycseqids2=Compare("class",classes1,classes2,classification1,classification2,output)
pseqids1,pseqids2,onlypseqids1,onlypseqids2=Compare("phylum",phyla1,phyla2,classification1,classification2,output)
kseqids1,kseqids2,onlykseqids1,onlykseqids2=Compare("kingdom",kingdoms1,kingdoms2,classification1,classification2,output)
outputfile.close()
#write better classification to file
newclassificationfilename1=classificationfilename1 + ".better"
newclassificationfilename2=classificationfilename2 + ".better"
CreateNewClassificationFile(classificationfilename1,newclassificationfilename1,species1.keys(),sseqids1,genera1.keys(),gseqids1,families1.keys(),fseqids1,orders1.keys(),oseqids1,classes1.keys(),cseqids1,phyla1.keys(),pseqids1,kingdoms1.keys(),kseqids1)
CreateNewClassificationFile(classificationfilename2,newclassificationfilename2,species2.keys(),sseqids2,genera2.keys(),gseqids2,families2.keys(),fseqids2,orders2.keys(),oseqids2,classes2.keys(),cseqids2,phyla2.keys(),pseqids2,kingdoms2.keys(),kseqids2)

##write only classification to file
newclassificationfilename1_by_1_only=classificationfilename1 + ".only"
newclassificationfilename2_by_2_only=classificationfilename2 + ".only"
CreateNewClassificationFile(classificationfilename1,newclassificationfilename1_by_1_only,species1.keys(),onlysseqids1,genera1.keys(),onlygseqids1,families1.keys(),onlyfseqids1,orders1.keys(),onlyoseqids1,classes1.keys(),onlycseqids1,phyla1.keys(),onlypseqids1,kingdoms1.keys(),onlykseqids1)
CreateNewClassificationFile(classificationfilename2,newclassificationfilename2_by_2_only,species2.keys(),onlysseqids2,genera2.keys(),onlygseqids2,families2.keys(),onlyfseqids2,orders2.keys(),onlyoseqids2,classes2.keys(),onlycseqids2,phyla2.keys(),onlypseqids2,kingdoms2.keys(),onlykseqids2)

print("The comparison is saved in file " + output + ".")
print("The better classifications are saved in files " + newclassificationfilename1 + " and " + newclassificationfilename2 + ".")
print("The only classifications are saved in files " + newclassificationfilename1_by_1_only + " and " + newclassificationfilename2_by_2_only + ".")

