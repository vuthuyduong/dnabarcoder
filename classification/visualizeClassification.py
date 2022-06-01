#!/usr/bin/env python
# FILE: visualizeClassification.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
#from sklearn.model_selection import StratifiedKFold
#from sklearn.preprocessing import StandardScaler
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
import json
from Bio import SeqIO
import random

parser=argparse.ArgumentParser(prog='visualizeClassification.py',  
							   usage="%(prog)s [options] -i assignmentfile -c classification -o kronareport",
							   description='''Script that visualizes using Krona for the classification results.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the assignment file')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', default ="", help='the classification file of the reference sequences.')

args=parser.parse_args()
predictionfilename=args.input
classificationfilename=args.classification
outputpath=args.out
if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]
	
def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	if "." in basename:
		basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + filename
	return path

def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

def GetTaxonomicClassification(level,header,texts):
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
	species="unidentified"
	genus="unidentified"
	family="unidentified"
	order="unidentified" 
	bioclass="unidentified"
	phylum="unidentified"
	kingdom="unidentified"
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
	taxonname=""
	rank=""
	if level <7 and kingdom!="unidentified":
		taxonname=kingdom
		rank="kingdom"
		classification="k__" + kingdom +";p__unidentified;c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified"	
	if level <6 and phylum!="unidentified":
		taxonname=phylum
		rank="phylum"
		classification="k__" + kingdom +";p__"+phylum +";c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified"
	if level <5 and bioclass!="unidentified":
		taxonname=bioclass
		rank="class"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass+";o__unidentified;f__unidentified;g__unidentified;s__unidentified"
	if level <4 and order!="unidentified":
		taxonname=order
		rank="order"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order + ";f__unidentified;g__unidentified;s__unidentified"
	if level <3 and family!="unidentified":
		taxonname=family
		rank="family"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family +";g__unidentified;s__unidentified"
	if level <2 and genus!="unidentified":
		taxonname=genus
		rank="genus"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family + ";g__"+ genus + "s__unidentified"
	if level <1 and species!="unidentified":
		taxonname=species
		rank="species"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family + ";g__"+ genus+";s__"+species 	
	return classification,taxonname,rank

def LoadClassification(classificationfilename,idcolumnname):
	classificationdict={}
	if not os.path.exists(classificationfilename):
		return classificationdict
	classificationfile=open(classificationfilename)
	header=next(classificationfile)
	seqidpos=-1
	isError=False
	texts=header.rstrip().split("\t")
	i=0
	for text in texts:
		if text.lower()==idcolumnname.lower():
			seqidpos=i
		i=i+1
	if 	seqidpos==-1:
		print("Please specify the sequence id columnname by using -idcolumnname.")
		isError=True	
	for line in classificationfile:
		texts=line.split("\t")
		seqid = texts[seqidpos].rstrip()
		classificationdict.setdefault(seqid,"")
		classification,taxonname,rank=GetTaxonomicClassification(0,header,texts)
		classificationdict[seqid]=classification
	classificationfile.close()	
	return classificationdict,isError

def LoadClassificationFromDescription(fastafilename):
	classificationdict={}
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	for seqid in seqrecords.keys():
		description=seqrecords[seqid].description
		species="s__unidentified"
		genus="g__unidentified"
		family="f__unidentified"
		order="o__unidentified"
		bioclass="c__unidentified"
		phylum="p__unidentified"
		kingdom="k__unidentified"
		if " " in description:
			description=description.split(" ")[1]
		texts=description.split("|")
		for text in texts:
			taxa=text.split(";")	
			for taxon in taxa:
				if taxon.startswith("k__"):
					kingdom=taxon
				elif taxon.startswith("p__"):
					phylum=taxon
				elif taxon.startswith("c__"):
					bioclass=taxon
				elif taxon.startswith("o__"):
					order=taxon
				elif taxon.startswith("f__"):
					family=taxon
				elif taxon.startswith("g__"):
					genus=taxon
				elif taxon.startswith("s__") and (" " in taxon.replace("s__","") or "_" in taxon.replace("s__","")):
					species=taxon	
		classification=kingdom + ";" + phylum + ";" + bioclass + ";" + order + ";" + family + ";" + genus + ";" + species
		classificationdict[seqid]=classification
	return classificationdict

def GetLevel(rank):
	level=-1
	if rank=="species":	
		level=6
	elif rank=="genus":	
		level=5
	elif rank=="family":	
		level=4
	elif rank=="order":	
		level=3
	elif rank=="class":	
		level=2
	elif rank=="phylum":
		level=1
	elif rank=="kingdom":
		level=0	
	return level

def LoadPrediction(predictionfilename,classificationdict):
	classificationdict={}
	predictionfile= open(predictionfilename, "r")
	next(predictionfile)
	for line in predictionfile:
		texts=line.split("\t")
		classification=texts[3]
		classification=classification.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","")
		if ";" in classification:
			texts=classification.split(";")
			classification=texts[0]
			if classification=="":
				classification="unidentified"
			for i in range(1,len(texts)):
				text=texts[i]
				if text=="":
					text="unidentified"
				classification=classification + "\t" + text	
		if classification in classificationdict.keys():
			classificationdict[classification]=classificationdict[classification] + 1
		else:
			classificationdict.setdefault(classification, 1)	
	return classificationdict

def KronaPieCharts(classification,kronareport,kronahtml):
	kronareportfile=open(kronareport,"w")
	for classname in classification.keys():
		kronareportfile.write(str(classification[classname]) + "\t" + classname + "\n")
	kronareportfile.close()	
	#create kronahtml
	command="ImportText.pl " + kronareport + " -o " + kronahtml
	#print(command)
	os.system(command)
	os.system("firefox " + kronahtml) 
	
if __name__ == "__main__":
	classificationdict=LoadClassification(classificationfilename)
	if is_fasta(classificationfilename):
		classificationdict=LoadClassificationFromDescription(classificationfilename)
	else:	
		classificationdict,isError=LoadClassification(classificationfilename,args.idcolumnname)
		if isError==True:
			sys.exit()
	classificationdict= LoadPrediction(predictionfilename,classificationdict)
	#making krona report
	kronareport = GetWorkingBase(predictionfilename) + ".krona.report"
	kronahtml=GetWorkingBase(kronareport) + ".html"
	KronaPieCharts(classificationdict,kronareport,kronahtml)
	print("The krona report and html are saved in files " + kronareport + " and " + kronahtml + ".") 
