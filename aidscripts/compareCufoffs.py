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

parser=argparse.ArgumentParser(prog='compareCutoffs.py',  
							   usage="%(prog)s [options] -i cutofffiles -o output",
							   description='''Script that compares the cutoffs from different cutoff files.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--filenames', required=True, help='the cutoff filenames separated by a comma')
parser.add_argument('-o','--out', required=True, help='The output.')

args=parser.parse_args()
filenames= args.filenames.split(",")
output=args.out

def LoadCutoffs(cutoffdict,cutofffilename):
	cutoffs={}
	if cutoffsfilename!="" and cutoffsfilename!=None:
		with open(cutoffsfilename) as cutoffsfile:
			cutoffs = json.load(cutoffsfile)
		for rank in cutoffs.keys():
			if not (rank in cutoffdict.keys()):
				cutoffdict.setdefault(rank,{})
			combineddatasets=cutoffdict[rank]
			datasets=cutoffs[rank]
			for datasetname in datasets.keys():
				dataset=datasets[datasetname]
				cutoff=dataset["cut-off"]
				confidence=dataset["confidence"]
				if not (datasetname in combineddatasets.keys()):
					combineddatasets.setdefault(datasetname,{})
				tmp=combineddatasets[datasetname]
				if not (cutofffilename in tmp.keys()):
					tmp.setdefault(cutofffilename,{})	
				tmpdataset=tmp[cutofffilename]	
				tmpdataset["cut-off"]=cutoff	
				tmpdataset["confidence"]=confidence

########MAIN###############

header="Rank\tTaxon name"
cutoffdict={}
for cutoffsfilename in filenames:
	LoadCutoffs(cutoffdict,cutoffsfilename)	
	header=header + "\t" + os.path.basename(cutoffsfilename) + " cut-off\t" + os.path.basename(cutoffsfilename) + " confidence"
header=header + "\n"
outputfile=open(output,"w")
outputfile.write(header)
confidencecounts=[0]*len(filenames)
cutoffcounts=[0]*len(filenames)
texts={}
cutofftexts={}
count=0

for rank in cutoffdict.keys():
	combineddatasets=cutoffdict[rank]
	for datasetname in combineddatasets.keys():
		tmp=combineddatasets[datasetname]
		line=rank + "\t" + datasetname
		maxconfidence=0
		maxcutoff=0
		bestcutoff=0
		bestregion=""
		maxregion=""
		check=0
		text=rank + "\t" + datasetname + "\t"
		confidencelist=[]
		cutofflist=[]
		for filename in filenames:
			cutoff=""
			confidence=""
			if filename in tmp.keys():
				check=check+1
				dataset=tmp[filename]
				cutoff=str(dataset["cut-off"])
				confidence=str(dataset["confidence"])
				confidencelist.append(float(confidence))
				cutofflist.append(float(cutoff))
				if float(confidence) > maxconfidence:
					maxconfidence=float(confidence)
					bestregion=filename
					bestcutoff=float(cutoff)
				if float(cutoff) > maxcutoff:
					maxcutoff=float(cutoff)
					maxregion=filename
			line=line+"\t" + cutoff+"\t"+confidence
		line=line+"\n"	
		outputfile.write(line)	
		if check==len(filenames) :
			#absolute confidence
			#text=text + bestregion  + "\t" + str(bestcutoff) + "\t" + str(maxconfidence)
			if confidencelist.count(maxconfidence)== 1:
				if not (bestregion in texts.keys()):
					texts.setdefault(bestregion,[line])
				else:
					texts[bestregion].append(line)
				confidencecounts[filenames.index(bestregion)]=confidencecounts[filenames.index(bestregion)] + 1
			if cutofflist.count(maxcutoff)== 1:
				if not (maxregion in cutofftexts.keys()):
					cutofftexts.setdefault(maxregion,[line])
				else:
					cutofftexts[maxregion].append(line)
				cutoffcounts[filenames.index(maxregion)]=cutoffcounts[filenames.index(maxregion)] + 1 	 
			count=count+1
		
outputfile.close()		
print("The comparison is saved in file " + output + ".")
print("Number of taxa have all regions available: " + str(count))
i=0
for filename in filenames:
	print("Number of taxa that " + filename + " has the highest resolving power: " + str(confidencecounts[i]) + " (" + str(round(confidencecounts[i]*100/count,2)) + "% )")
	if filename in texts.keys():
		suboutputname=output + "." + os.path.basename(filename) + ".higherconfidence"
		suboutputfile=open(suboutputname,"w")
		#suboutputfile.write("Rank\tTaxon name\tcut-off\tconfidence\n")
		suboutputfile.write(header)
		for line in texts[filename]:
			suboutputfile.write(line)
		suboutputfile.close()	
		print("The taxa that " + filename + " has the highest resolving power are saved in file " + suboutputname)
	i=i+1
i=0
for filename in filenames:
	print("Number of taxa that " + filename + " has the highest similarity cut-offs: " + str(cutoffcounts[i]) + " (" + str(round(cutoffcounts[i]*100/count,2)) + "% )")
	if filename in cutofftexts.keys():
		suboutputname=output + "." + os.path.basename(filename) + ".highercutoff"
		suboutputfile=open(suboutputname,"w")
		#suboutputfile.write("Rank\tTaxon name\tcut-off\tconfidence\n")
		suboutputfile.write(header)
		for line in cutofftexts[filename]:
			suboutputfile.write(line)
		suboutputfile.close()	
		print("The taxa that " + filename + " has the highest cutoff are saved in file " + suboutputname)
	i=i+1	
outputfile.close()
	
		