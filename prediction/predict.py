#!/usr/bin/env python
# -*- coding: utf-8 -*-
# FILE: predictlocally.py
# AUTHOR: Duong Vu
# CREATE DATE: 08 oct 2020

import os, argparse
import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
plt.rc('font',size=6)
#from matplotlib.patches import Polygon
import numpy as np
import multiprocessing
import json

parser=argparse.ArgumentParser(prog='predict.py', 
							   usage="%(prog)s [options] -i fastafile -c classificationfile -p classificationposition -st startingthreshold -et endthreshold -s step",
							   description='''Script that predicts an optimal threshold to separate the sequences based on the given classification''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-prefix','--prefix', help='the prefix of output filenames')
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-p','--classificationpositions', default="", help='the classification positions for the prediction, separated by ",".')
parser.add_argument('-st','--startingthreshold', type=float, default=0, help='starting threshold')
parser.add_argument('-et','--endthreshold', type=float, default=0, help='ending threshold')
parser.add_argument('-s','--step', type=float, default=0.001, help='the step to be increased for the threshold after each step of the prediction.')
parser.add_argument('-mc','--mincoverage', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, mc should probably be set to 100.')
parser.add_argument('-sim','--simfilename', help='The similarity matrix of the sequences if exists.')
parser.add_argument('-hp','--higherclassificationpositions', default="", help='The prediction is based on the whole dataset if hp="". Otherwise it will be predicted based on different datasets obtained at the higher classifications, separated by ",".')
parser.add_argument('-minGroupNo','--minimumgroupnumber', type=int, default=10, help='The minimum number of groups needed for prediction.')
parser.add_argument('-minSeqNo','--minimumsequencenumber', type=int, default=50, help='The minimum number of sequences needed for prediction.')
#parser.add_argument('-type','--predictiontype', default="global", help='The type of prediction. There are three options for prediction type: global, local and all.')
parser.add_argument('-redo','--redo', default="", help='Recompute F-measure for the current parameters.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
classificationpos=args.classificationpositions
threshold = args.startingthreshold
endthreshold=args.endthreshold
step=args.step
mincoverage=args.mincoverage
outputfolder=args.out
simfilename=args.simfilename
higherclassificationpos=args.higherclassificationpositions
minGroupNo=args.minimumgroupnumber
minSeqNo=args.minimumsequencenumber
#predictiontype=args.predictiontype
prefix=args.prefix
outputpath=args.out
redo=args.redo

if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)


nproc=multiprocessing.cpu_count()

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]  
	
def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	if "." in basename:
		basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + filename
	return path

class PointDef:
	def __init__(self, id,flag,neighbors,seqrecord):
		self.id=id
		self.flag=flag
		self.neighbors = neighbors
		self.seqrecord=seqrecord

class ClusterDef:
	def __init__(self, id, pointids):
		self.id=id
		self.pointids=pointids

def GetSeqIndex(seqname,seqlist):
	i=0
	for seq in seqlist:
		if (seqname == seq.id):
			return i
		i = i + 1
	return -1
	
def LoadSim(simfilename):
	simmatrix = {} #we use dictionary to reduce the memory constraints 
	simfile = open(simfilename)
	seqids=[]
	for line in simfile:
		numbers=line.rstrip().split(" ")
		i=numbers[0]
		j=numbers[1]
		seqids.append(i)
		seqids.append(j)
		if i not in simmatrix.keys():
			simmatrix.setdefault(i, {})
		simmatrix[i][j]=float(numbers[2])
	seqids=list(set(seqids))	
	for seqid1 in seqids:
		if not (seqid1 in simmatrix.keys()):
			simmatrix.setdefault(seqid1, {})
		simmatrix[seqid1][seqid1]=1	
		for seqid2 in seqids:
			if not seqid2 in simmatrix[seqid1].keys():
				simmatrix[seqid1][seqid2]=0
	simfile.close()		
	return simmatrix

def SaveSim(simmatrix,simfilename):
	simfile=open(simfilename,"w")
	for i in simmatrix.keys():
		for j in simmatrix[i].keys():
			simfile.write(str(i) + " " + str(j) + " " + str(simmatrix[i][j]) + "\n")
	simfile.close()
	
def ComputeSim(fastafilename,seqrecords,mincoverage):
	simmatrix={}
	for seqid1 in seqrecords.keys():
		simmatrix.setdefault(seqid1,{})
		for seqid2 in seqrecords.keys():
			if seqid1==seqid2:
				simmatrix[seqid1][seqid2]=1
			else:
				simmatrix[seqid1][seqid2]=0	
	#blast
	makedbcommand = "makeblastdb -in " + fastafilename + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + fastafilename + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	if mincoverage >=300:
		blastcommand = "blastn -query " + fastafilename + " -db  db -outfmt 6 -out out.txt -num_threads " + str(nproc)
	os.system(blastcommand)
	
	#read blast output
	blastoutputfile = open("out.txt")
	score=0
	for line in blastoutputfile:
		if line.rstrip()=="":
			continue
		words = line.split("\t")
		i = words[0].rstrip()
		j = words[1].rstrip()
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		score=sim
		if coverage < mincoverage:
			score=float(score * coverage)/mincoverage	
		if simmatrix[i][j] < score:
			simmatrix[i][j]=round(score,4)
			simmatrix[j][i]=round(score,4)
		#simmatrix[j][i]=score
	os.system("rm out.txt")
	return simmatrix

def LoadNeighbors(seqids,simmatrix,threshold):
	neighbordict={}
	for seqid in seqids:
		neighbordict.setdefault(seqid, [])
	for i in seqids:
		for j in seqids:
			if simmatrix[i][j] >= threshold:
				neighbordict[i].append(j)
				neighbordict[j].append(i)
	#os.system("rm out.txt")
	return neighbordict

def LoadPoints(neigbordict,seqrecords):
	points={}
	for seqid in  seqrecords.keys():
		point = PointDef(seqid,False,neigbordict[seqid],seqrecords[seqid])
		points[seqid]=point
	return points

def ExpandCluster(root,cluster,points):
	for i in root.neighbors:
		if points[i].flag==False:
			points[i].flag=True
			cluster.pointids.append(i)
			ExpandCluster(points[i],cluster,points)		

def Cluster(points,clusters):
	for pointid in points.keys():
		if points[pointid].flag==False:
			points[pointid].flag=True
			cluster = ClusterDef(len(clusters),[])
			cluster.pointids.append(pointid)
			ExpandCluster(points[pointid], cluster,points)
			clusters.append(cluster)
			
def ComputeFmeasure(classes,clusters):
	#compute F-measure
	f=0
	n=0
	for classname in classes.keys():
		group=classes[classname]
		m = 0
		for cluster in clusters:
			i = len(set(group) & set(cluster.pointids))
			v = float(2*i)/float((len(group) + len(cluster.pointids)))
			if m < v:
				m=v		
		n = n + len(group)
		f = f +	(len(group)*m)	
	return float(f)/float(n) 

def LoadClasses(seqids,classificationfilename,pos):
	classificationfile= open(classificationfilename)
#	records= open(classificationfilename,errors='ignore')
	classification={}
	classes={}
	for line in classificationfile:
		if line.startswith("#"):
			 continue
		words=line.split("\t")
		seqid= words[0].rstrip().replace(">","")
		if seqid in seqids:
			classname=""
			if pos < len(words):
				classname=words[pos].rstrip()
			classification[seqid]=classname
			if classname in classes.keys():
				classes[classname].append(seqid)
			else:
				classes.setdefault(classname,[seqid]) 
	classificationfile.close()
	return classes,classification

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False
	
def Predict(prediction_datasetname,records,classes):
	thresholds=[]
	fmeasures=[]	
	t=round(threshold,4)
	optthreshold=0
	bestFmeasure=0
	fmeasuredict={}
	if 'cut-off' in prediction_datasetname.keys():	
		optthreshold=prediction_datasetname['cut-off']
	if 'confidence' in prediction_datasetname.keys():		
		bestFmeasure=prediction_datasetname['confidence']
	if 'fmeasures' in prediction_datasetname.keys():			
		fmeasuredict=prediction_datasetname['fmeasures']
	#compute optimal threshold
	while t <= endthreshold:
		print("Computing F-measure for threshold " + str(t))
		fmeasure=0
		#if str(t) in fmeasuredict.keys() and datasetdict['fasta filename']==fastafilename and datasetdict['classification filename']==classificationfilename:
		if str(t) in fmeasuredict.keys() and args.redo=="":
			fmeasure=fmeasuredict[str(t)]
		else:	
			#compute fmeasure
			neighbordict = LoadNeighbors(records.keys(),simmatrix,t)	
			points=LoadPoints(neighbordict,records)
			clusters=[]
			Cluster(points,clusters)	
			fmeasure=round(ComputeFmeasure(classes,clusters),4)
			fmeasuredict[str(t)]=fmeasure
		if fmeasure > bestFmeasure or (fmeasure==bestFmeasure and optthreshold >t) :
			bestFmeasure=fmeasure
			optthreshold=t
		thresholds.append(t)
		fmeasures.append(fmeasure)
		t=round(t+step,4)
		print("F-measure: " + str(fmeasure))
	prediction_datasetname['cut-off']=optthreshold
	prediction_datasetname['confidence']=bestFmeasure
	prediction_datasetname['sequence number']=len(records)
	prediction_datasetname['group number']=len(classes)
	prediction_datasetname['fmeasures']=fmeasuredict
	return thresholds,fmeasures,optthreshold,bestFmeasure

#def GetHigherPositionList(classificationfilename,positionlist,ranklist,higherclassificationpos,predictiontype):
#	higherpositionlist=[]
#	if predictiontype=="global":
#		higherpositionlist=[0]*len(positionlist)
#	elif "," in higherclassificationpos:	
#		positions=higherclassificationpos.split(",")
#		for pos in positions:
#			higherpositionlist.append(int(pos))
#	elif higherclassificationpos !="":
#		higherpositionlist.append(int(higherclassificationpos))
#	elif predictiontype=="local":
#		classificationfile=open(classificationfilename)
#		header=classificationfile.readline()
#		classificationfile.close()
#		texts=header.rstrip().split("\t")
#		p_s=0
#		p_g=0
#		p_f=0
#		p_o=0
#		p_c=0
#		p_p=0
#		p_k=0
#		i=0
#		for text in texts:
#			if text.lower()=="species":
#				p_s=i
#			elif text.lower()=="genus":
#				p_g=i	
#			elif text.lower()=="family":
#				p_f=i	
#			elif text.lower()=="order":
#				p_o=i	
#			elif text.lower()=="class":
#				p_c=i	
#			elif text.lower()=="phylum":
#				p_p=i	
#			elif text.lower()=="kingdom":
#				p_k=i	
#			i=i+1	
#		for pos in positionlist:
#			higherpos=0
#			if pos==p_s:
#				higherpos=p_g
#			elif pos==p_g:
#				higherpos=p_f
#			elif pos==p_f:
#				higherpos=p_o
#			elif pos==p_o:
#				higherpos=p_c
#			elif pos==p_c:
#				higherpos=p_p
#			elif pos==p_p:
#				higherpos=p_k	
#			higherpositionlist.append(higherpos)	 
#	return higherpositionlist

#def GetHigherPositions(classificationfilename,pos):
#	higherpositionlist=[]
#	classificationfile=open(classificationfilename)
#	header=classificationfile.readline()
#	classificationfile.close()
#	texts=header.rstrip().split("\t")
#	p_s=-1
#	p_g=-1
#	p_f=-1
#	p_o=-1
#	p_c=-1
#	p_p=-1
#	p_k=-1
#	i=0
#	for text in texts:
#		if text.lower()=="species":
#			p_s=i
#		elif text.lower()=="genus":
#			p_g=i	
#		elif text.lower()=="family":
#			p_f=i	
#		elif text.lower()=="order":
#			p_o=i	
#		elif text.lower()=="class":
#			p_c=i	
#		elif text.lower()=="phylum":
#			p_p=i	
#		elif text.lower()=="kingdom":
#			p_k=i	
#		i=i+1	
#	level=0	
#	if pos==p_s:
#		level=5
#	elif pos==p_g:
#		level=4
#	elif pos==p_f:
#		level=3
#	elif pos==p_o:
#		level=2
#	elif pos==p_c:
#		level=1
#	if level >4:
#		higherpositionlist.append(p_g)	 
#	if level >3:
#		higherpositionlist.append(p_f)	 		
#	if level >2:
#		higherpositionlist.append(p_o)	 	
#	if level >1:
#		higherpositionlist.append(p_c)	 	
#	if level >0:
#		higherpositionlist.append(p_p)	 
#	return higherpositionlist

#def GenerateDatasetsForPrediction(seqrecords,classificationfilename,pos,higherpos):
#	datasets={}
#	#load classification
#	allseqids=[]
#	classificationfile= open(classificationfilename)
#	higherclassnames=[]
#	classnames=[]
#	for line in classificationfile:
#		if line.startswith("#"):
#			continue 		
#		texts=line.split("\t")
#		seqid=texts[0].replace(">","").rstrip()
#		classname=""
#		higherclassname=""
#		if pos < len(texts):
#			classname=texts[pos].rstrip()
#		if higherpos > 0:
#			if higherpos < len(texts):
#				higherclassname=texts[higherpos].rstrip()	
#		else:
#			higherclassname="All"
#		if classname != "" and higherclassname !="":
#			allseqids.append(seqid)
#			higherclassnames.append(higherclassname)
#		if not classname in classnames:
#			classnames.append(classname)
#	classificationfile.close()
#	for seqid in seqrecords.keys():
#		if seqid in allseqids:
#			higherclassname=higherclassnames[allseqids.index(seqid)]
#			if not higherclassname in datasets.keys():
#				datasets.setdefault(higherclassname,{})
#			datasets[higherclassname][seqid]=seqrecords[seqid]	
#	return datasets

#def GenerateDatasets(seqrecords,classificationfilename,pos,predictiontype):
#	alldatasets={}
#	if predictiontype=="local":		
#		higherpositionlist=GetHigherPositions(classificationfilename,pos)
#		for higherpos in higherpositionlist:
#			datasets=GenerateDatasetsForPrediction(seqrecords,classificationfilename,pos,higherpos)
#			alldatasets.update(datasets)
#		return alldatasets	
#	else:
#		alldatasets.setdefault("All",{})
#		classificationfile= open(classificationfilename)
#		for line in classificationfile:
#			texts=line.rstrip().split("\t")
#			seqid=texts[0].replace(">","")
#			classname=""
#			if pos < len(texts):
#				classname=texts[pos].rstrip()	
#			if classname != "" and seqid in seqrecords.keys():
#				alldatasets["All"][seqid]=seqrecords[seqid]
#		classificationfile.close()
#	return alldatasets

def GetPositionList(classificationfilename,classificationpos):
	positionlist=[]
	ranklist=[]
	if "," in classificationpos:
		positions=classificationpos.split(",")
		for pos in positions:
			positionlist.append(int(pos))
	elif classificationpos !="":
		positionlist.append(int(classificationpos))
	classificationfile=open(classificationfilename)
	header=classificationfile.readline()
	texts=header.split("\t")
	for pos in positionlist:
		if pos <len(texts):
			ranklist.append(texts[pos].rstrip())
		else:
			ranklist.append("position " + str(pos))
	classificationfile.close()
	return positionlist,ranklist
	
def GenerateDatasetsForPrediction(seqrecords,classificationfilename,pos,higherpos):
	datasets={}
	classificationfile= open(classificationfilename)
	for line in classificationfile:
		if line.startswith("#"):
			continue 		
		texts=line.rstrip().split("\t")
		seqid=texts[0].replace(">","")
		classname=""
		higherclassname=""
		if pos < len(texts):
			classname=texts[pos].rstrip()
		if higherpos > 0:
			if higherpos < len(texts):
				higherclassname=texts[higherpos].rstrip()	
		if classname != "" and higherclassname !="" and classname != "unidentified" and higherclassname !="unidentified" and seqid in seqrecords.keys():
			if not higherclassname in datasets.keys():
				datasets.setdefault(higherclassname,{})
			datasets[higherclassname][seqid]=seqrecords[seqid]	
	classificationfile.close()
	return datasets

def GenerateDatasets(seqrecords,classificationfilename,pos,higherclassificationpos):
	alldatasets={}
	if higherclassificationpos !="":		
		higherpositionlist,higherranklist=GetPositionList(classificationfilename,higherclassificationpos)
		for higherpos in higherpositionlist:
			datasets=GenerateDatasetsForPrediction(seqrecords,classificationfilename,pos,higherpos)
			alldatasets.update(datasets)
		return alldatasets	
	else:
		alldatasets.setdefault("All",{})
		classificationfile= open(classificationfilename)
		for line in classificationfile:
			texts=line.rstrip().split("\t")
			seqid=texts[0].replace(">","")
			classname=""
			if pos < len(texts):
				classname=texts[pos].rstrip()	
			if classname != "" and classname != "unidentified" and seqid in seqrecords.keys():
				alldatasets["All"][seqid]=seqrecords[seqid]
		classificationfile.close()
	return alldatasets

def LoadPrediction(predictionfilename):
	existingprediction={}
	#load classes
	with open(predictionfilename,encoding='latin1') as json_file:
		existingprediction = json.load(json_file)
	return existingprediction
			
def LoadPredictionAtPos(prediction_datasetname):
	thresholds=[]
	fmeasures=[]	
	optthreshold=prediction_datasetname['cut-off']
	bestFmeasure=prediction_datasetname['confidence']
	seqno=prediction_datasetname['sequence number']
	groupno=prediction_datasetname['group number']
	fmeasuredict=prediction_datasetname['fmeasures']
	for t in fmeasuredict.keys():
		thresholds.append(float(t))
		fmeasures.append(fmeasuredict[t])
	#sorting
	keydict = dict(zip(fmeasures,thresholds))
	fmeasures.sort(key=keydict.get)
	thresholds.sort()
	return thresholds,fmeasures,optthreshold,bestFmeasure,seqno,groupno

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
			del dataset["fmeasures"]
			seqno=dataset["sequence number"]
			groupno=dataset["group number"]
			if groupno < minGroupNo or seqno < minSeqNo:	#delete the cutoffs that dont have enough sequences and groups for prediction
				del datasets[datasetname]
	with open(outputnamewithoutfmeasures,"w") as json_file:
		if sys.version_info[0] >= 3:
			json.dump(finalresults,json_file,indent=2)
		else:
			json.dump(finalresults,json_file,encoding='latin1',indent=2)
	#save as tab. format
	textoutput=outputnamewithoutfmeasures+".txt"
	textfile=open(textoutput,"w")
	textfile.write("Rank\tDataset\tcut-off\tconfidence\tsequence number\tgroup number\tfasta filename\tclassification filename\tclassification position\n")
	cutoff=0
	confidence=0
	SeqNo=0
	GroupNo=0
	fastafilename=""
	classificationfilename=""
	classificationposition=-1
	for rank in finalresults.keys():
		datasets=finalresults[rank]
		for datasetname in datasets.keys():
			dataset=datasets[datasetname]
			cutoff=dataset["cut-off"]
			confidence=dataset["confidence"]
			SeqNo=dataset["sequence number"]
			GroupNo=dataset["group number"]
			fastafilename=dataset["fasta filename"]
			classificationfilename=dataset["classification filename"]
			classificationposition=dataset["classification position"]
			textfile.write(rank+"\t" + datasetname + "\t"+str(cutoff)+"\t"+str(confidence)+"\t"+str(SeqNo)+"\t"+str(GroupNo)+"\t"+fastafilename+"\t"+classificationfilename+"\t"+str(classificationposition)+"\n")
	textfile.close()
			

def GetTaxa(classificationfilename,pos,highertaxon):
	taxa=[]
	classificationfile=open(classificationfilename)
	for line in classificationfile:
		if line.startswith("#"):
			continue		
		texts=line.rstrip().split("\t")
		if not highertaxon in texts:
			continue
		if pos<len(texts):
			taxon=texts[pos]
			if not taxon in taxa:
				taxa.append(taxon)
	classificationfile.close()
	return taxa
	
def PlotPrediction(thresholdlist,fmeasurelist,optthresholds,bestFmeasures,features,datasetnames,figoutput):
	fig, ax = plt.subplots(figsize=(3,3)) 
	if len(thresholdlist) >1:
		fig, ax = plt.subplots(figsize=(6,3))
	ax.set_xlabel("Cut-off")
	ax.set_ylabel('F-measure')
	colors = plt.cm.rainbow(np.linspace(0, 1, len(thresholdlist)))
	labels=[]
	i=0
	for thresholds in thresholdlist:
		fmeasures=fmeasurelist[i]
		optthreshold=optthresholds[i]
		bestFmeasure=bestFmeasures[i]
		ax.plot(np.array(thresholds), np.array(fmeasures),color=colors[i])
		labels.append(features[i] + " cut-off for " + datasetnames[i] + ": "  + str(round(optthresholds[i],4)))
		ax.text(round(optthreshold,4), 0.97, round(bestFmeasure,4), transform=ax.get_xaxis_transform(), horizontalalignment='center', size='x-small', color=colors[i])
		i=i+1
	ax.set_title("Predicting similarity cut-offs for sequence identification")	
	plt.legend(labels, loc="top left", bbox_transform=plt.gcf().transFigure)
	plt.tight_layout()
	plt.savefig(figoutput, dpi = 500)
	plt.show()	
	
def PlotResults(optthresholds,bestFmeasures,features,datasetnames,localfigoutput):
	#sort all according to increasing order of optthresholds
	if len(optthresholds)==0:
		return
	optthresholds,bestFmeasures,datasetnames,features = zip(*sorted(zip(optthresholds,bestFmeasures,datasetnames,features)))
	labels=[]
	i=0
	for datasetname in datasetnames:
		if len(set(features))==1:
			labels.append(datasetnames[i] )
		else:
			labels.append(features[i] + " cut-off for " + datasetnames[i] )
		i=i+1
	x = np.arange(len(labels))  # the label locations
	#width = 0.35  # the width of the bars
	fig, ax = plt.subplots(figsize=(6,3)) 
	ax2=ax.twinx()
	#ax.bar(x,optthresholds,color='b')
	ax.plot(np.array(labels), np.array(optthresholds), color='b')
	ax2.plot(np.array(labels), np.array(bestFmeasures), color='r')
	ax.set_ylabel('Cut-off',color='b')
	ax2.set_ylabel("F-measure",color="r")
	if len(set(features))==1:
		ax.set_title(features[0] + " cut-offs and F-measures predicted for different groups")
	else:
		ax.set_title("Cut-offs and F-measures predicted for different groups")
	
	ax.set_xticks(x)
	ax.set_xticklabels(labels,rotation=90)
	ax.legend()
	plt.tight_layout()
	plt.savefig(localfigoutput, dpi = 500)
	plt.show()	
	
if __name__ == "__main__":
	if prefix=="" or prefix==None:
		basename=os.path.basename(fastafilename)
		prefix=basename[:-(len(basename)-basename.rindex("."))] 
	#load sequences
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	if len(seqrecords) >100:
		sys.setrecursionlimit(len(seqrecords)*2)
	if simfilename=="" or simfilename==None:
		simfilename=GetWorkingBase(prefix) + ".sim"
		
	#load similarity matrix
	simmatrix={}
	#load or compute simmatrix
	if endthreshold >=threshold and endthreshold >0:
		if os.path.exists(simfilename):
			print("Loading similarity matrix " + simfilename)
			simmatrix=LoadSim(simfilename)
		else:	
			print("Computing similarity matrix...")
			simmatrix=ComputeSim(fastafilename,seqrecords,mincoverage)
			print("Save similarity matrix " + simfilename)
			SaveSim(simmatrix,simfilename)	
	
	outputname=GetWorkingBase(prefix) + ".predicted"	
	predictiondict={}
	if os.path.exists(outputname):
		predictiondict=LoadPrediction(outputname)
	#load positions of classification for prediction
	positionlist,ranklist=GetPositionList(classificationfilename,classificationpos)
	if len(positionlist)==0:
		for rank in predictiondict.keys():
			prediction_datasets=predictiondict[rank]
			for datasetname in prediction_datasets.keys():
				datasetdict=prediction_datasets[datasetname]
				if "fasta filename" in datasetdict.keys() and "classification filename" in datasetdict.keys():
					if fastafilename== datasetdict["fasta filename"] and classificationfilename==datasetdict["classification filename"]:
						pos=datasetdict["classification position"]
						if not (pos in positionlist):
							positionlist.append(pos)
							ranklist.append(rank)
#	#load higher classification positions for loading the datasets for prediction
#	higherpositionlist=GetHigherPositionList(classificationfilename,positionlist,ranklist,higherclassificationpos,predictiontype)
#	if len(higherpositionlist)!=len(positionlist):
#		print("Please give a list of higher positions (-hp) for loading the datasets for the prediction.")
#		sys.exit()	
	thresholdlist=[]
	intrathresholdlist=[]
	fmeasurelist=[]
	optthresholds=[]
	bestFmeasures=[]
	features=[]
	datasetnames=[]
	i=0
	#predicting
	for pos in positionlist:
		rank=ranklist[i]
		prediction_datasets={}
		if rank in predictiondict.keys():
			prediction_datasets=predictiondict[rank]
		if step==0 or (endthreshold < threshold) or (endthreshold==0 and threshold==0):
			#load existing prediction
			for datasetname in prediction_datasets.keys():
				datasetdict=prediction_datasets[datasetname]
				if not ("fasta filename" in datasetdict.keys() and "classification filename" in datasetdict.keys()):
					continue
				elif fastafilename != datasetdict["fasta filename"] and classificationfilename!=datasetdict["classification filename"]:
					continue
#				if predictiontype=="global" and datasetname !="All":
#					continue
				seqno=datasetdict['sequence number']
				groupno=datasetdict['group number']			
				thresholds,fmeasures,optthreshold,bestFmeasure,seqno,groupno=LoadPredictionAtPos(datasetdict)
				if not (seqno < minGroupNo or seqno < minSeqNo):	#for visualization
					thresholdlist.append(thresholds)
					fmeasurelist.append(fmeasures)
					optthresholds.append(optthreshold)
					bestFmeasures.append(bestFmeasure)
					features.append(rank)
					datasetnames.append(datasetname + "(" + str(seqno) + "," + str(groupno) + ")")		
		else:
			#predicting
			#higherpos=higherpositionlist[i]
			#datasets=GenerateDatasetsForPrediction(seqrecords,classificationfilename,pos,higherpos)
			datasets=GenerateDatasets(seqrecords,classificationfilename,pos,higherclassificationpos)
			for datasetname in datasets.keys():
#				if predictiontype=="global" and datasetname !="All":
#					continue
				records=datasets[datasetname]
				seqno=len(records)
				#load classification at given positions 	
				classes,classification=LoadClasses(records.keys(),classificationfilename,pos)
#				#only predict when the numbers of the groups > 1
				if len(classes) < 2:
					continue
				datasetdict={}
				if datasetname in prediction_datasets.keys():
					datasetdict = prediction_datasets[datasetname]
				thresholds=[]
				fmeasures=[]		
				optthreshold=0
				bestFmeasure=0
				groupno=len(classes)
				print("Predicting optimal threshold to separate sequences at the " + rank + " level for " + datasetname)
				thresholds,fmeasures,optthreshold,bestFmeasure=Predict(datasetdict,records,classes)	
				datasetdict['fasta filename']=fastafilename
				datasetdict['classification filename']=classificationfilename
				datasetdict['classification position']=pos
				if not (datasetname in prediction_datasets.keys()):
					prediction_datasets[datasetname]=datasetdict		
				if not (groupno < minGroupNo or seqno < minSeqNo):	#for visualization
					thresholdlist.append(thresholds)
					fmeasurelist.append(fmeasures)
					optthresholds.append(optthreshold)
					bestFmeasures.append(bestFmeasure)
					features.append(rank)
					datasetnames.append(datasetname + "(" + str(seqno) + "," + str(groupno) + ")")
			if not rank in predictiondict.keys():	
				predictiondict[rank] = prediction_datasets 
		i=i+1
	#save prediction and cut-offs for classification and identification
	if len(positionlist) >0:	
		if len(predictiondict.keys())>0:
			outputwithoutfmeasures=GetBase(outputname) + ".cutoffs"	
			SavePrediction(predictiondict,outputname,outputwithoutfmeasures)
			#SaveCutoffs(predictiondict,outputcutoffs)
			print("The prediction and cutoffs are saved in the files " + outputname + " and " + outputwithoutfmeasures + ".")	
	else:
		#load existing prediction for plotting
		for rank in predictiondict.keys():
			prediction_datasets=predictiondict[rank]
			for datasetname in prediction_datasets.keys():
				prediction_datasetname=prediction_datasets[datasetname]
#				if predictiontype=="global" and datasetname !="All":
#					continue
				thresholds,fmeasures,optthreshold,bestFmeasure,feature,seqno,groupno=LoadPredictionAtPos(prediction_datasetname)
				if not (groupno < minGroupNo or seqno < minSeqNo):	#for visualization
					thresholdlist.append(thresholds)
					fmeasurelist.append(fmeasures)
					optthresholds.append(optthreshold)
					bestFmeasures.append(bestFmeasure)
					features.append(feature)
					datasetnames.append(datasetname + "(" + str(seqno) + "," + str(groupno) + ")")
	#Plot the prediction results
	globalfigoutput=""
	localfigoutput=""
	if len(positionlist) == 1:
		globalfigoutput=GetBase(outputname) + "." + str(positionlist[0]) + ".global.png"
		localfigoutput=GetBase(outputname) + "." + str(positionlist[0]) + ".local.png"
	elif len(positionlist) > 1 or len(positionlist) == 1: 	
		globalfigoutput=GetBase(outputname) + ".global.png"
		barplotfigoutput=GetBase(outputname) + ".local.png"
	if higherclassificationpos=="":	
		if len(thresholdlist) >0:
			#plot all predictions			
			PlotPrediction(thresholdlist,fmeasurelist,optthresholds,bestFmeasures,features,datasetnames,globalfigoutput)	
			print("The prediction plot is saved in the file " + globalfigoutput + ".")
		else:
			print("Please check the parameters.")
	else:	
		if len(optthresholds) >0:
			#barplot the prediction results only
			PlotResults(optthresholds,bestFmeasures,features,datasetnames,localfigoutput)
			print("The prediction plot is saved in the file " + localfigoutput + ".")	
		else:
			print("Please check the parameters.")		
	
