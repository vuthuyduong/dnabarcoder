#!/usr/bin/env python
# FILE: mergeCutoffs.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020


import sys, argparse
import json
#from copy import deepcopy

parser=argparse.ArgumentParser(prog='mergeCutoffs.py',  
							   usage="%(prog)s [options] -i listofdict -o outputname",
							   description='''The script that merges a list of dictionaries to one dictionary. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='A list of dictionaries separated by commas.')
parser.add_argument('-o','--out', help='The merged dictionary.')


args=parser.parse_args()
dictionarylist= args.input
outputfilename=args.out


def SaveCutoffs(mergeddict,outputfilename):
	#save mergeddict		
	with open(outputfilename,"w") as json_file:
		if sys.version_info[0] >= 3:
			json.dump(mergeddict,json_file,indent=2)
		else:
			json.dump(mergeddict,json_file,encoding='latin1',indent=2)			
	#save as tab. format
	textoutput=outputfilename+".txt"
	textfile=open(textoutput,"w")
	textfile.write("Rank\tDataset\tcut-off\tconfidence\tsequence number\tgroup number\tmin alignment length\tfasta filename\tclassification filename\n")
	cutoff=0
	confidence=0
	SeqNo=0
	GroupNo=0
	fastafilename=""
	classificationfilename=""
	#classificationposition=-1
	for rank in mergeddict.keys():
		datasets=mergeddict[rank]
		for datasetname in datasets.keys():
			dataset=datasets[datasetname]
			cutoff=0
			if "cut-off" in dataset.keys():
				cutoff=dataset["cut-off"]
			confidence=0
			if "confidence" in dataset.keys():
				confidence=dataset["confidence"]
			SeqNo=0
			if "sequence number" in dataset.keys():
				SeqNo=dataset["sequence number"]
			GroupNo=0	
			if "group number" in dataset.keys():
				GroupNo=dataset["group number"]
			minalignmentlength=0
			if "min alignment length" in dataset.keys():
				minalignmentlength=dataset["min alignment length"]	
			fastafilename=""	
			if "fasta filename" in dataset.keys():	
				fastafilename=dataset["fasta filename"]
			classificationfilename=""	
			if 	"classification filename" in dataset.keys():
				classificationfilename=dataset["classification filename"]
			#classificationposition=""	
			#if 	"classification position" in dataset.keys():
			#	classificationposition=dataset["classification position"]
			textfile.write(rank+"\t" + datasetname + "\t"+str(cutoff)+"\t"+str(confidence)+"\t"+str(SeqNo)+"\t"+str(GroupNo)+"\t"+ str(minalignmentlength) + "\t" + fastafilename+"\t"+classificationfilename+"\n")
	textfile.close()
	print("The outputs are saved in " + outputfilename + " and " + textoutput + ".")
####MAIN####
dictionaries=[]
if "," in dictionarylist:
	dictionaries=dictionarylist.split(",")
else:
	dictionaries.append(dictionarylist)	
mergeddict={}	
for dictionaryname in dictionaries:
	cutoffdict={}
	#load classes
	with open(dictionaryname,encoding='latin1') as json_file:
		cutoffdict = json.load(json_file)
	if cutoffdict!={}:	
		for rank in cutoffdict.keys():
			if not (rank in mergeddict.keys()):
				mergeddict.setdefault(rank,{})
			mergeddatasets=mergeddict[rank]	
			datasets=cutoffdict[rank]
			for datasetname in datasets.keys():
				dataset=datasets[datasetname]
				confidence=0
				if "confidence" in dataset.keys():
					confidence=float(dataset["confidence"])
				if not (datasetname in mergeddatasets.keys()):
					mergeddatasets.setdefault(datasetname,{})
				mergeddataset=mergeddatasets[datasetname]
				mergedconfidence=0
				if "confidence" in mergeddataset.keys():
					mergedconfidence=float(mergeddataset["confidence"])
				if mergedconfidence < confidence:
					mergeddataset=dataset.copy()
				mergeddatasets[datasetname]=mergeddataset	
SaveCutoffs(mergeddict,outputfilename)


