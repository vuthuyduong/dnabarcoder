#!/usr/bin/env python
# FILE: computeVariationLocally.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import numpy as np
import os
from Bio import SeqIO
import json
import random
import multiprocessing
nproc=multiprocessing.cpu_count()

referencename=sys.argv[1]
classificationfilename = sys.argv[2]
classificationposition =int(sys.argv[3])
mincoverage=300
if len(sys.argv) >4:
	mincoverage = float(sys.argv[4])
jsonvariationfilename =""
if len(sys.argv) >5:
	jsonvariationfilename  = sys.argv[5]
	
def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def LoadClassification(seqIDs,seqrecords,classificationfilename,pos):
	classification=[""]*len(seqIDs)
	classes=[]
	classnames=[]
	preclassids=[]
	preclassnames=[]
	postclassids=[]
	postclassnames=[]
	
	if classificationfilename == "":
		return classification
	records= list(open(classificationfilename, "r"))
	for line in records:
		if line.startswith("#"):
			continue 
		elements=line.split("\t")
		seqid = elements[0].replace(">","").rstrip()
		if seqid in seqIDs:
			index=seqIDs.index(seqid)
			classname=elements[pos].rstrip()
			if classname in classnames:
				classid=classnames.index(classname)
				classes[classid].append(seqrecords[index])
			else:
				classnames.append(classname)
				seqs=[]
				seqs.append(seqrecords[index])
				classes.append(seqs)
			preclassname=""	
			if pos>=1:
				preclassname=elements[pos-1].rstrip()
			postclassname=""
			if pos+1<len(elements):
				postclassname=elements[pos+1].rstrip()
			if  preclassname !="" and preclassname in preclassnames:
				classid=preclassnames.index(preclassname)
				preclassids[classid].append(seqid)
			elif preclassname !="":
				preclassnames.append(preclassname)
				seqs=[]
				seqs.append(seqid)
				preclassids.append(seqs)	
			if postclassname !="" and postclassname in postclassnames:
				classid=postclassnames.index(postclassname)
				postclassids[classid].append(seqid)
			elif postclassname!="":
				postclassnames.append(postclassname)
				seqs=[]
				seqs.append(seqid)
				postclassids.append(seqs)		
	return classification,classes,classnames,preclassids,postclassids

def GetSeqIndex(seqname,seqrecords):
	i=0
	for seqrecord in seqrecords:
		if (seqname == seqrecord.id):
			return i
		i = i + 1
	return -1

def ComputeBLASTscoreMatrix(fastafilename,records,mincoverage):
	scorematrix = [[0 for x in range(len(records))] for y in range(len(records))] 
	seqids=[]
	for rec in records:
		seqids.append(rec.id)
	#blast
	makedbcommand = "makeblastdb -in " + fastafilename + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + fastafilename + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	os.system(blastcommand)
	
	#read blast output
	if not os.path.isfile("out.txt"):
		return scorematrix
	blastoutputfile = open("out.txt")
	refid = ""
	score=0
	queryid=""
	for line in blastoutputfile:
		if line.rstrip()=="":
			continue
		words = line.split("\t")
		queryid = words[0].rstrip()
		refid = words[1].rstrip()
		i = seqids.index(queryid)
		j = seqids.index(refid)
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		score=sim
		if coverage < mincoverage:
			score=float(score * coverage)/mincoverage
		if scorematrix[i][j] < score:
			scorematrix[i][j]=score
			scorematrix[j][i]=score
	os.system("rm out.txt")
	return scorematrix

def ComputeVariation(reffilename,mincoverage):
	#load sequeces from the fasta files
	records = list(SeqIO.parse(reffilename, "fasta"))
	scorematrix=ComputeBLASTscoreMatrix(reffilename,records,mincoverage)
	scorelist=[]
	for i in range(0,len(scorematrix)-2):
		for j in range(i+1,len(scorematrix)-1):
			if i!=j:
				scorelist.append(scorematrix[i][j])
	threshold=1
	minthreshold=1		
	if len(scorelist) >0:
		x=np.array(scorelist)
		minthreshold=np.min(x)
		threshold=np.median(x)
	return threshold,minthreshold

def ComputeVariations(variationfilename,classes,classnames,mincoverage):
	#create json dict
	variations={}
	i=0
	for taxonname in classnames:
		sequences=classes[i]
		if len(sequences) >0:
			threshold=0
			minthreshold=0
			if len(sequences) < 100:
				fastafilename=taxonname.replace(" ","_") + ".fasta"
				SeqIO.write(sequences,fastafilename,"fasta")
				threshold,minthreshold=ComputeVariation(fastafilename,mincoverage)
				os.system("rm " + fastafilename)
			else:
				threshold,minthreshold=EvaluateVariation(taxonname,sequences,mincoverage)
			currentvariation={'Threshold': threshold, 'MinThreshold': minthreshold,'NumberOfSequences': len(sequences)}
			variations[taxonname]=currentvariation
		i=i+1	
	#write to file
	with open(variationfilename,"w") as json_file:
		json.dump(variations,json_file,encoding='latin1')	
	return variations

def EvaluateVariation(taxonname,sequences,mincoverage):
	thresholds=[]
	minthresholds=[] 
	for i in range(0,10):
		n=int(len(sequences)/10)
		selectedindexes=random.sample(range(0, len(sequences)), k=n)
		selectedsequences=[]
		for index in selectedindexes:
			selectedsequences.append(sequences[index])
		fastafilename=taxonname.replace(" ","_") + ".fasta"
		SeqIO.write(selectedsequences,fastafilename,"fasta")
		threshold,minthreshold=ComputeVariation(fastafilename,mincoverage)
		os.system("rm " + fastafilename)
		thresholds.append(threshold)
		minthresholds.append(minthreshold)
	threshold = np.median(np.array(thresholds))
	minthreshold =  np.median(np.array(minthresholds))
	return threshold,minthreshold

def EvaluateVariations(variationfilename,classes,classnames,mincoverage):
	variations={}
	thresholds=[]
	minthresholds=[] 
	for i in range(0,len(classnames)):
		thresholds.append([])
		minthresholds.append([])
	for i in range(0,10):
		selectedclasses=SelectRandomSequences(classes)
		tmpvariations=ComputeVariations(selectedclasses,classnames,mincoverage)
		j=0
		for classname in classnames:
			thresholds[j].append(tmpvariations[classname]['Threshold'])
			minthresholds[j].append(tmpvariations[classname]['MinThreshold'])	
			j=j+1
	j=0		
	for classname in classnames:
		threshold = np.median(np.array(thresholds[j]))
		minthreshold =  np.median(np.array(minthresholds[j]))
		seqno=len(classes[j])
		currentvariation={'Threshold': threshold, 'MinThreshold': minthreshold,'NumberOfSequences': seqno}
		variations[classname]=currentvariation
		j=j+1
	#write to file
	with open(variationfilename,"w") as json_file:
		json.dump(variations,json_file,encoding='latin1')
	return variations
		
def IndexSequences(filename):
	indexedfilename = GetBase(filename) + ".indexed.fasta"
	fastafile = open(filename)
	indexedfile = open(indexedfilename, "w")
	i=0
	for line in fastafile:
		if line.startswith('>'):
			indexedfile.write(">" + str(i) + "|" + line.rstrip()[1:] + "\n")
			i=i+1
		else:
			indexedfile.write(line)    
	fastafile.close()
	indexedfile.close()
	return indexedfilename

def SelectRandomSequences(classes):
	selectedclasses=[]
	selectedids=[]
	for tmpclass in classes:
		#select randomly 1/10 of the sequences to compare
		if len(tmpclass) >=20:
			n=int(len(tmpclass)/10)
			selectedids=random.sample(range(0, len(tmpclass)), k=n)
			selectedclass=[]
			for i in selectedids:
					selectedclass.append(tmpclass[i])
			selectedclasses.append(selectedclass)		
		else:
			selectedclasses.append(tmpclass)
	return selectedclasses

def SelectSequencesUsingChildClassification(childclassids,classes):
	selectedclasses=[]
	selectedids=[]
	for childclass in childclassids:
		i=random.randint(0,len(childclass)-1)
		selectedids.append(childclass[i])
	for tmpclass in classes:
		selectedclass=[]
		for seq in tmpclass:
			if seq.id in selectedids:
				selectedclass.append(seq)
		if len(selectedclass)==0:
			selectedclass=tmpclass
		selectedclasses.append(selectedclass)		
	return selectedclasses

def SaveVariationInTabFormat(output,variation):
	outputfile=open(output,"w")
	outputfile.write("Taxonname\tThreshold\tMin threshold\tNumber of sequences\n")	
	for classname in variation.keys():
		threshold=variation[classname]['Threshold']
		minthreshold=variation[classname]['MinThreshold']
		seqno=variation[classname]['NumberOfSequences']
		outputfile.write(classname + "\t" + str(threshold) + "\t" + str(minthreshold) + "\t" + str(seqno) + "\n")
	outputfile.close()
##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]
#load train seq records
referencerecords =  list(SeqIO.parse(referencename, "fasta"))
referenceIDs=[]
for seq in referencerecords:
	referenceIDs.append(seq.id)
if jsonvariationfilename == "":
	jsonvariationfilename = GetBase(referencename) + "." + str(classificationposition) + ".variation"

#Load classes, classification:
referenceclassification,classes,classnames,preclassids,postclassids= LoadClassification(referenceIDs,referencerecords,classificationfilename, classificationposition)
#selectedclasses=[]
#if len(preclassids) > len(postclassids):
#	selectedclasses=SelectSequencesForComputingVariation(preclassids,classes)
#else:
#	selectedclasses=SelectSequencesForComputingVariation(postclassids,classes)
variations=ComputeVariations(jsonvariationfilename,classes,classnames,mincoverage)
SaveVariationInTabFormat(jsonvariationfilename + ".txt",variations)
print("The result is saved in the json file  " + jsonvariationfilename + " and tab file " + jsonvariationfilename + ".txt." )


