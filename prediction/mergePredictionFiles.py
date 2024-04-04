#!/usr/bin/env python
# FILE: mergeCutoffs.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020


import sys, argparse
import json
#from copy import deepcopy

parser=argparse.ArgumentParser(prog='mergePredictionFiles.py',  
							   usage="%(prog)s [options] -i listofjson_prediction_files -o outputname",
							   description='''The script that merges a list of prediction files (.predicted) in json format, separated by commas, to one prediction file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='A list of prediction filenames (.predicted) separated by commas.')
parser.add_argument('-o','--out', help='The merged file of the prediction files.')
parser.add_argument('-mingroupno','--mingroupno', type=int, default=10, help='The minimum number of groups needed for prediction.')
parser.add_argument('-minseqno','--minseqno', type=int, default=30, help='The minimum number of sequences needed for prediction.')
parser.add_argument('-maxproportion','--maxproportion', type=float, default=1, help='Only predict when the proportion of the sequences the largest group of the dataset is less than maxproportion. This is to avoid the problem of inaccurate prediction due to imbalanced data.')

args=parser.parse_args()
dictionarylist= args.input
outputfilename=args.out

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]  

def SavePrediction(predictiondict,outputnamewithfmeasures,outputnamewithoutfmeasures):
	#save the whole prediction file
	with open(outputnamewithfmeasures,"w") as json_file:
		if sys.version_info[0] >= 3:
			json.dump(predictiondict,json_file,indent=2)
		else:
			json.dump(predictiondict,json_file,encoding='latin1',indent=2)
	#save only the final results (without fmeasures)
	finalresults=predictiondict.copy()
	for rank in finalresults.keys():
		datasets=finalresults[rank]
		datasetnames=list(datasets.keys()).copy()
		for datasetname in datasetnames:
			dataset=datasets[datasetname]
			if "fmeasures" in dataset.keys(): #remove prediction before saving cutoffs
				del dataset["fmeasures"]
			seqno=0 
			if "sequence number" in dataset.keys():
				seqno=dataset["sequence number"]
			groupno=0
			if "group number" in dataset.keys():
				groupno=dataset["group number"]
			minalignmentlength=0
			if "min alignment length" in dataset.keys():
				minalignmentlength=dataset["min alignment length"]
			maxproportion =0
			if "max proportion" in dataset.keys():
				maxproportion=dataset["max proportion"]	
			#if groupno < minGroupNo or seqno < minSeqNo or (minalignmentlength>0 and minalignmentlength <mincoverage):	#delete the cutoffs that dont have enough sequences and groups for prediction
			if groupno < args.mingroupno or seqno < args.minseqno or maxproportion > args.maxproportion :	#remove the cutoffs that dont have enough sequences and groups for prediction
				del datasets[datasetname]			
	with open(outputnamewithoutfmeasures,"w") as json_file:
		if sys.version_info[0] >= 3:
			json.dump(finalresults,json_file,indent=2)
		else:
			json.dump(finalresults,json_file,encoding='latin1',indent=2)
	#save as tab. format
	textoutput=outputnamewithoutfmeasures+".txt"
	textfile=open(textoutput,"w")
	textfile.write("Rank\tDataset\tcut-off\tconfidence\tsequence number\tgroup number\tmin alignment length\tmax proportion\tfasta filename\tclassification filename\n")
	for rank in finalresults.keys():
		datasets=finalresults[rank]
		for datasetname in datasets.keys():
			dataset=datasets[datasetname]
			cutoff=0 
			if "cut-off" in dataset.keys():
				cutoff=dataset["cut-off"]
			confidence=0 
			if "confidence" in dataset.keys():
				confidence=dataset["confidence"]	
			seqno=0 
			if "sequence number" in dataset.keys():
				seqno=dataset["sequence number"]
			groupno=0
			if "group number" in dataset.keys():
				groupno=dataset["group number"]
			minalignmentlength=0
			if "min alignment length" in dataset.keys():
				minalignmentlength=dataset["min alignment length"]
			fastafilename=""
			maxproportion =0
			if "max proportion" in dataset.keys():
				maxproportion=dataset["max proportion"]
			if "fasta filename" in dataset.keys():
				fastafilename=dataset["fasta filename"]
			classificationfilename=""
			if "classification filename" in dataset.keys():	
				classificationfilename=dataset["classification filename"]
			#classificationposition=dataset["classification position"]
			textfile.write(rank+"\t" + datasetname + "\t"+str(cutoff)+"\t"+str(confidence)+"\t"+str(seqno)+"\t"+str(groupno)+"\t"+ str(minalignmentlength) + "\t" + str(maxproportion) + "\t" + fastafilename+"\t"+classificationfilename+"\n")
	textfile.close()
	
####MAIN####
dictionaries=[]
if "," in dictionarylist:
	dictionaries=dictionarylist.split(",")
else:
	dictionaries.append(dictionarylist)	
mergeddict={}	
for dictionaryname in dictionaries:
	predictiondict={}
	#load classes
	with open(dictionaryname,encoding='latin1') as json_file:
		predictiondict = json.load(json_file)
	if predictiondict!={}:	
		for rank in predictiondict.keys():
			if not (rank in mergeddict.keys()):
				mergeddict.setdefault(rank,{})
			mergeddatasets=mergeddict[rank]	
			datasets=predictiondict[rank]
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
outputwithoutfmeasures=GetBase(outputfilename) + ".cutoffs.json"					
SavePrediction(mergeddict,outputfilename,outputwithoutfmeasures)


