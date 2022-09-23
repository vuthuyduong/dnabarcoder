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
#from sklearn.metrics import precision_recall_fscore_support
#from sklearn.metrics import cohen_kappa_score
#from sklearn.metrics import matthews_corrcoef
#from sklearn.metrics import confusion_matrix
#from sklearn.metrics import accuracy_score
import json
from Bio import SeqIO
import random

parser=argparse.ArgumentParser(prog='visualizeClassification.py',  
							   usage="%(prog)s [options] -i assignmentfile -c classification -o kronareport",
							   description='''Script that visualizes using Krona for the classification results.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the assignment/classification file')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')

args=parser.parse_args()
predictionfilename=args.input
outputpath=args.out
if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

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
		classification= kingdom +"\tunidentified\tunidentified\tunidentified\tunidentified\tunidentified\tunidentified"
	if level <6 and phylum!="unidentified":
		taxonname=phylum
		rank="phylum"
		classification= kingdom +"\t"+phylum +"\tunidentified\tunidentified\tunidentified\tunidentified\tunidentified"
	if level <5 and bioclass!="unidentified":
		taxonname=bioclass
		rank="class"
		classification=kingdom +"\t"+phylum +"\t"+bioclass+"\tunidentified\tunidentified\tunidentified\tunidentified"
	if level <4 and order!="unidentified":
		taxonname=order
		rank="order"
		classification=kingdom +"\t"+phylum +"\t"+bioclass +"\t"+ order + "\tunidentified\tunidentified\tunidentified"
	if level <3 and family!="unidentified":
		taxonname=family
		rank="family"
		classification="\t" + kingdom +"\t"+phylum +"\t"+bioclass +"\t"+ order+"\t"+family +"\tunidentified\tunidentified"
	if level <2 and genus!="unidentified":
		taxonname=genus
		rank="genus"
		classification=kingdom +"\t"+phylum +"\t"+bioclass +"\t"+ order+"\t"+family + "\t"+ genus + "\tunidentified"
	if level <1 and species!="unidentified":
		taxonname=species
		rank="species"
		classification=kingdom +"\t"+phylum +"\t"+bioclass +"\t"+ order+"\t"+family + "\t"+ genus+"\t"+species
	return classification

def GetTaxonomicClassificationFromDescription(texts):
	species = "unidentified"
	genus = "unidentified"
	family = "unidentified"
	order = "unidentified"
	bioclass = "unidentified"
	phylum = "unidentified"
	kingdom = "unidentified"
	for text in texts:
		text=text.rstrip()
		taxa=[]
		if ";" in text:
			taxa = text.split(";")
		else:
			taxa.append(text)
		for taxon in taxa:
			if taxon.startswith("k__"):
				kingdom = taxon.replace("k__","")
			elif taxon.startswith("p__"):
				phylum = taxon.replace("p__","")
			elif taxon.startswith("c__"):
				bioclass = taxon.replace("c__","")
			elif taxon.startswith("o__"):
				order = taxon.replace("o__","")
			elif taxon.startswith("f__"):
				family = taxon.replace("f__","")
			elif taxon.startswith("g__"):
				genus = taxon.replace("g__","")
			elif taxon.startswith("s__") and (" " in taxon.replace("s__", "") or "_" in taxon.replace("s__", "")):
				species = taxon.replace("s__","")
	classification = kingdom + "\t" + phylum + "\t" + bioclass + "\t" + order + "\t" + family + "\t" + genus + "\t" + species
	return classification

def LoadPrediction(predictionfilename):
	classificationdict={}
	predictionfile= open(predictionfilename, "r")
	header=next(predictionfile)
	header=header.lower()
	areTaxaSeparatedByTab=False
	if (" species " in header) or (" genus " in header) or (" family " in header) or (" order " in header) or (" class " in header) or (" phylum " in header):
		areTaxaSeparatedByTab = True
	for line in predictionfile:
		texts=line.split("\t")
		classification=""
		if areTaxaSeparatedByTab==True:
			classification = GetTaxonomicClassification(0, header, texts)
		else:
			classification = GetTaxonomicClassificationFromDescription(texts)
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
	classificationdict = LoadPrediction(predictionfilename)
	#making krona report
	kronareport = GetBase(predictionfilename) + ".krona.report"
	kronahtml=GetBase(kronareport) + ".html"
	KronaPieCharts(classificationdict,kronareport,kronahtml)
	print("The krona report and html are saved in files " + kronareport + " and " + kronahtml + ".") 
