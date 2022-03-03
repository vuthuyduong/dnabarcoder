#!/usr/bin/env python
# FILE: addFeaturesToTaxa.py, requiring Python 2.7 if classification file contains non-ascii characters
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys, argparse
import numpy as np
import os
from Bio import SeqIO
import json
import multiprocessing
nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='addFeaturesTotTaxa.py',  
							   usage="%(prog)s [options] -i taxafile -ip taxapositionintheinputfile -f featurefile -fp featurepositionsinthefeaturefile -o output",
							   description='''Script that adds features of the taxa given in the features file to the taxa in the taxafile. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the input file containing taxa in a tab limited format.')
parser.add_argument('-ip','--inputpos', required=True, type =int, help='the position in the input file containing the taxa.')
parser.add_argument('-o','--out', help='the output file.') 
parser.add_argument('-f','--features', required=True, help='the features file containing features associated with some taxa in a tab. limited format.')
parser.add_argument('-fp','--featurepositions', required=True, default="", help='the positions of the feature in the feature file.')
parser.add_argument('-unique','--unique', default="", help='If unique=="" then the most dominant feature is selected, otherwise all features are selected.')

args=parser.parse_args()
inputfilename= args.input
featurefilename=args.features
pos=args.inputpos
featurepositions=args.featurepositions
unique=args.unique
output=args.out

def LoadFeatures(featurefilename,featurepositions):
	featuredict={}
	featurefile=open(featurefilename,errors="ignore")
	header=next(featurefile)
	poslist=[]
	if "," in featurepositions:
		texts=featurepositions.split(",")
		for text in texts:
			poslist.append(int(text))
	else:
		poslist.append(int(featurepositions))
	featurelabellist=[]
	for featurepos in poslist:
		featurelabel=header.split("\t")[featurepos]
		featurelabel=featurelabel.rstrip()
		featurelabellist.append(featurelabel)
	for line in featurefile:
		texts=line.split("\t")
		for text in texts:
			text=text.rstrip()
			if text=="":
				continue
			if not (text in featuredict.keys()):
				featuredict.setdefault(text,{})
			subfeaturedict=featuredict[text]
			for featurepos in poslist:
				if not (featurepos in subfeaturedict.keys()):
					subfeaturedict.setdefault(featurepos,{})
				feature=""	
				if featurepos < len(texts):	
					feature=texts[featurepos].rstrip()
				if feature!="":
					if (not (feature in subfeaturedict[featurepos].keys())):
						subfeaturedict[featurepos].setdefault(feature,0)
					subfeaturedict[featurepos][feature]=subfeaturedict[featurepos][feature]+1
	featurefile.close()				
	return featuredict,featurelabellist,poslist

#####MAIN
featuredict,featurelabellist,poslist=LoadFeatures(featurefilename,featurepositions)	
outputfile=open(output,"w")
inputfile=open(inputfilename,errors="ignore")
header=next(inputfile)
header=header.replace("\n","")
for label in featurelabellist:
	header=header + "\t" + label 
outputfile.write(header + "\n")
for line in inputfile:
	if pos>=len(line.split("\t")):
		continue
	taxon=line.split("\t")[pos].rstrip()
	line=line.replace("\n","")
	if taxon in featuredict.keys():
		subfeaturedict=featuredict[taxon]
		for featurepos in poslist:
			feature=""
			featurelist=subfeaturedict[featurepos]
			if unique=="":
				i=0
				for f in featurelist.keys():
					feature=feature + f
					if i<len(featurelist.keys())-1:
						feature=feature+","
			else:
				m=0
				for f in featurelist.keys():
					if m<featurelist[f]:
						m=featurelist[f]
				for f in featurelist.keys():
					if m==featurelist[f]:
						feature=feature + f + ","
				feature=feature[:-1]		
			line=line + "\t" + feature 
	outputfile.write(line + "\n")
inputfile.close()
outputfile.close()	
print("The output is saved in file " + output + ".")	

