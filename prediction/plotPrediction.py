#!/usr/bin/env python
# -*- coding: utf-8 -*-
# FILE: plotPrediction.py
# AUTHOR: Duong Vu
# CREATE DATE: 08 oct 2020

import os, argparse
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
plt.rc('font',size=6)
import numpy as np
import multiprocessing
import json


nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='plotPredict.py', 
							   usage="%(prog)s [options] -i listofpredictionfiles -t taxonname -st startingthreshold -et endthreshold -rank taxonomic level",
							   description='''Script that plot predictions from different predictionfiles for a taxonomic clade''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the list of prediction files, separated by commas.')
parser.add_argument('-o','--out', required=True, help='The output folder.')
parser.add_argument('-biomarkers','--biomarkers',default="", help='The list of biomarkers associated with the prediction files.')
parser.add_argument('-labelstyle','--labelstyle', default='normal', help='The label style to be displayed: normal, italic, or bold.')
parser.add_argument('-rank','--classificationranks', default="", help='the classification ranks for the prediction, separated by ",".')
parser.add_argument('-t','--taxa', default="All", help='the list of taxa that we want to load the predictions for, separated by commas.')
parser.add_argument('-mingroupno','--mingroupno', type=int, default=10, help='The minimum number of groups needed for prediction.')
parser.add_argument('-minseqno','--minseqno', type=int, default=30, help='The minimum number of sequences needed for prediction.')
parser.add_argument('-maxproportion','--maxproportion', type=float, default=1, help='Only predict when the proportion of the sequences the largest group of the dataset is less than maxproportion. This is to avoid the problem of inaccurate prediction due to imbalanced data.')
parser.add_argument('-mincutoff','--mincutoff', type=float, default=0, help='The minimum cutoff for selection.')
parser.add_argument('-st','--startingthreshold', type=float, default=0, help='starting threshold')
parser.add_argument('-et','--endthreshold', type=float, default=0, help='ending threshold')
parser.add_argument('-display','--display',default="", help='If display=="yes" then the plot figure is displayed.')
parser.add_argument('-label','--label',default="", help='The label to display in the figure.')
#parser.add_argument('-labelstyle','--labelstyle', default='normal', help='The label style to be displayed: normal, italic, or bold.')


args=parser.parse_args()
figoutput=args.out
taxa=args.taxa
classificationranks=args.classificationranks


def LoadPrediction(predictionfilename):
	existingprediction={}
	#load classes
	with open(predictionfilename,encoding='latin1') as json_file:
		existingprediction = json.load(json_file)
	return existingprediction
			
def LoadPredictionForGivenRankAndDataset(prediction_datasetname):
	thresholds=[]
	fmeasures=[]	
	optthreshold=0
	if 'cut-off' in prediction_datasetname.keys():
		optthreshold=prediction_datasetname['cut-off']
	bestFmeasure=0
	if 'confidence' in prediction_datasetname.keys():
		bestFmeasure=prediction_datasetname['confidence']
	seqno=0
	if 'sequence number' in prediction_datasetname.keys():
		seqno=prediction_datasetname['sequence number']
	groupno=0	
	if 'group number' in prediction_datasetname.keys():
		groupno=prediction_datasetname['group number']
	minalignmentlength=0	
	if 'min alignment length' in prediction_datasetname.keys():	
		minalignmentlength=prediction_datasetname['min alignment length']	
	minalignmentlength=0	
	if 'min alignment length' in prediction_datasetname.keys():	
		minalignmentlength=prediction_datasetname['min alignment length']		
	maxproportion=0	
	if 'max proportion' in prediction_datasetname.keys():	
		maxproportion=prediction_datasetname['max proportion']		
	fmeasuredict={}
	if 'fmeasures' in prediction_datasetname.keys():
		fmeasuredict=prediction_datasetname['fmeasures']
	for t in fmeasuredict.keys():
		if float(t)>=args.startingthreshold and ((float(t)<=args.endthreshold and args.endthreshold>0) or args.endthreshold==0):
			thresholds.append(float(t))
			fmeasures.append(fmeasuredict[t])
	#sorting
	keydict = dict(zip(fmeasures,thresholds))
	fmeasures.sort(key=keydict.get)
	thresholds.sort()
	return thresholds,fmeasures,optthreshold,bestFmeasure,seqno,groupno,minalignmentlength,maxproportion

def PlotPrediction(label,thresholdlist,fmeasurelist,optthresholds,bestFmeasures,features,datasetnames, biomarkerlist, figoutput):
	if len(thresholdlist) >5:
		fig, ax = plt.subplots(figsize=(6,3))
	else:	
		fig, ax = plt.subplots(figsize=(4.3,3)) 
	ax.set_xlabel("Cut-off")
	ax.set_ylabel('F-measure')
	colors = plt.cm.rainbow(np.linspace(0, 1, len(thresholdlist)))
	labels=[]
	i=0
	for thresholds in thresholdlist:
		biomarker=biomarkerlist[i]
		fmeasures=fmeasurelist[i]
		optthreshold=optthresholds[i]
		bestFmeasure=bestFmeasures[i]
		ax.plot(np.array(thresholds), np.array(fmeasures),color=colors[i])
		#labels.append(features[i] + " " + biomarker + " cut-off for " + datasetnames[i] + ": "  + str(round(optthresholds[i],4)))
		labels.append(features[i] + " " + biomarker + " cut-off " + str(round(optthreshold,4)) +  " with highest F-measure " +  str(round(bestFmeasure,4)) + " for " + datasetnames[i])
		#ax.text(round(optthreshold,4), 0.97, round(bestFmeasure,4), transform=ax.get_xaxis_transform(), horizontalalignment='center', size='x-small', color=colors[i])
		i=i+1
	#ax.set_title(datasetname + " Predicting similarity cut-offs for sequence identification")	
	ax.set_title(label)	
	plt.legend(labels, loc="lower left", bbox_transform=plt.gcf().transFigure,prop={'style': args.labelstyle})
	plt.tight_layout()
	plt.savefig(figoutput, dpi = 500)
	print("The prediction plot is saved in file " + figoutput + ".")
	if args.display=="yes":
		plt.show()	
	
def PlotResults(prefix,optthresholds,bestFmeasures,features,datasetnames,localfigoutput):
	if len(optthresholds)==0:
		return
	bestFmeasures,optthresholds,datasetnames,features = zip(*sorted(zip(bestFmeasures,optthresholds,datasetnames,features)))	
	labels=[]
	i=0
	for datasetname in datasetnames:
		if len(set(features))==1:
			labels.append(datasetnames[i]) 
		else:
			labels.append(features[i] + " cut-off for " + datasetnames[i] )
		i=i+1
	x = np.arange(len(labels))  # the label locations
	#width = 0.35  # the width of the bars
	fig, ax = plt.subplots(figsize=(6,2.2)) 
	ax2=ax.twinx()
	#ax.bar(x,optthresholds,color='b')
	ax.plot(np.array(labels), np.array(optthresholds), color='b')
	ax2.plot(np.array(labels), np.array(bestFmeasures), color='r')
	ax.set_ylabel('Cut-off',color='b')
	ax2.set_ylabel("F-measure",color="r")
	if len(set(features))==1:
		ax.set_title(prefix + ": " + features[0] + " cut-offs and F-measures predicted for different groups")
	else:
		ax.set_title(prefix + ": cut-offs and F-measures predicted for different groups")
	
	ax.set_xticks(x)
	#ax.set_xticklabels(labels,rotation=90, style=args.labelstyle)
	ax.set_xticklabels(labels,rotation=90, style='italic')
	ax.legend()
	plt.tight_layout()
	plt.savefig(localfigoutput, dpi = 500)
	print("The prediction plot is saved in file " + localfigoutput + ".")
	if args.display=="yes":
		plt.show()	
	
if __name__ == "__main__":
	
	predictiondict={}
	
	predictionfilenames=[args.input]
	biomarkers=[args.biomarkers]
	if "," in args.input:
		predictionfilenames=args.input.split(",")
		if "," in args.biomarkers:
			biomarkers=args.biomarkers.split(",")
		else:
			biomarkers=[]
			for i in range(len(predictionfilenames)):
				biomarkers.append(args.biomarkers)
	
	ranklist=[classificationranks]
	if "," in classificationranks:
		ranklist=classificationranks.split(",")
	taxalist=[taxa]	
	if "," in taxa:
		taxalist=taxa.split(",")
		

	thresholdlist=[]
	intrathresholdlist=[]
	fmeasurelist=[]
	optthresholds=[]
	bestFmeasures=[]
	features=[]
	datasetnames=[]
	
	
	#load existing prediction for plotting
	i=0
	biomarkerlist=[]
	for predictionfilename in predictionfilenames:
		predictiondict=LoadPrediction(predictionfilename)
		biomarker=biomarkers[i]
		for rank in ranklist:
			prediction_datasets=predictiondict[rank]
			for datasetname in taxalist:
				prediction_datasetname=prediction_datasets[datasetname]
				thresholds,fmeasures,optthreshold,bestFmeasure,seqno,groupno,minalignmentlength,maxproportion=LoadPredictionForGivenRankAndDataset(prediction_datasetname)
				#if not (groupno < minGroupNo or seqno < minSeqNo or (minalignmentlength >0 and minalignmentlength < mincoverage)):	#for visualization
				if not (groupno < args.mingroupno or seqno < args.minseqno or maxproportion > args.maxproportion > args.maxproportion or optthreshold <args.mincutoff):	#for visualization
					thresholdlist.append(thresholds)
					fmeasurelist.append(fmeasures)
					optthresholds.append(optthreshold)
					bestFmeasures.append(bestFmeasure)
					features.append(rank)
					datasetnames.append(datasetname + "(" + str(seqno) + "," + str(groupno) + ")")
					biomarkerlist.append(biomarker)
		i=i+1			
	
				
	PlotPrediction(args.label,thresholdlist,fmeasurelist,optthresholds,bestFmeasures,features,datasetnames,biomarkerlist, figoutput)				
				
	