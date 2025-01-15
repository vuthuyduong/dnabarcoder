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
import numpy as np
import multiprocessing
import json
import random

nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='predict.py', 
							   usage="%(prog)s [options] -i fastafile -c classificationfile -p classificationposition -st startingthreshold -et endthreshold -s step -ml minalignmentlength",
							   description='''Script that predicts an optimal threshold to separate the sequences based on the given classification''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-prefix','--prefix', default="", help='the prefix of output filenames.')
parser.add_argument('-label','--label',default="", help='The label to display in the figure.')
parser.add_argument('-labelstyle','--labelstyle', default='normal', help='The label style to be displayed: normal, italic, or bold.')
parser.add_argument('-c','--classification', default="", help='the classification file in tab. format.')
parser.add_argument('-rank','--classificationranks', default="", help='the classification ranks for the prediction, separated by ",".')
parser.add_argument('-st','--startingthreshold', type=float, default=0, help='starting threshold')
parser.add_argument('-et','--endthreshold', type=float, default=0, help='ending threshold')
parser.add_argument('-s','--step', type=float, default=0.001, help='the step to be increased for the threshold after each step of the prediction.')
parser.add_argument('-ml','--minalignmentlength', type=int,default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should probably be set to smaller, 50 for instance.')
parser.add_argument('-sim','--simfilename', help='The similarity matrix of the sequences if exists.')
#parser.add_argument('-hp','--higherclassificationpositions', default="", help='The prediction is based on the whole dataset if hp="". Otherwise it will be predicted based on different datasets obtained at the higher classifications, separated by ",".')
parser.add_argument('-higherrank','--higherclassificationranks', default="", help='The prediction is done on the whole dataset if higherranks="". Otherwise it will be predicted for different datasets obtained at the higher classifications, separated by ",".')
parser.add_argument('-mingroupno','--mingroupno', type=int, default=10, help='The minimum number of groups needed for prediction.')
parser.add_argument('-minseqno','--minseqno', type=int, default=30, help='The minimum number of sequences needed for prediction.')
parser.add_argument('-maxseqno','--maxseqno', type=int, default=20000, help='Maximum number of the sequences of the predicted taxon name from the classification file will be selected for the comparison to find the best match. If it is not given, all the sequences will be selected.')
parser.add_argument('-maxproportion','--maxproportion', type=float, default=1, help='Only predict when the proportion of the sequences the largest group of the dataset is less than maxproportion. This is to avoid the problem of inaccurate prediction due to imbalanced data.')
parser.add_argument('-mincutoff','--mincutoff', type=float, default=0, help='The minimum cutoff for selection.')
parser.add_argument('-taxa','--taxa', default="", help='The selected taxa separated by commas for local prediction. If taxa=="", all the clades at the given higher positions are selected for prediction.')
parser.add_argument('-removecomplexes','--removecomplexes',default="", help='If removecomplexes="yes", indistinguishable groups will be removed before the prediction.')
parser.add_argument('-redo','--redo', default="", help='Recompute F-measure for the current parameters.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')
parser.add_argument('-display','--display',default="", help='If display=="yes" then the plot figure is displayed.')
parser.add_argument('-ncpus','--ncpus', type=int, default=nproc, help='The number of CPUs used for searching. The default value is the total number of CPUs.')


args=parser.parse_args()
fastafilename= args.input
classificationfilename=args.classification
#classificationpos=args.classificationpositions
#higherclassificationpos=args.higherclassificationpositions
classificationranks=args.classificationranks
higherclassificationranks=args.higherclassificationranks
threshold = args.startingthreshold
endthreshold=args.endthreshold
step=args.step
mincoverage=args.minalignmentlength
outputfolder=args.out
simfilename=args.simfilename
minGroupNo=args.mingroupno
minSeqNo=args.minseqno
prefix=args.prefix
label=args.label
outputpath=args.out
redo=args.redo
nproc=args.ncpus

if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)




def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]  
	
#def GetWorkingBase(filename):
#	basename=os.path.basename(filename)
#	if "." in basename:
#		basename=basename[:-(len(basename)-basename.rindex("."))] 
#	path=outputpath + "/" + basename
#	return path
	
def GetWorkingBase(filename): 
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
		#if len(seqids) < args.maxSimMatrixSize: #load full matrix
		for seqid2 in seqids:
			if not seqid2 in simmatrix[seqid1].keys():
				#simmatrix[seqid1].setdefault(seqid2,0)
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
	blastoutput = fastafilename + ".blast.out"		
	blastdb=fastafilename + ".db"		
	#blast
	print("Comparing the sequences of " + fastafilename + " using Blast...")
	makedbcommand = "makeblastdb -in " + fastafilename + " -dbtype \'nucl\' " +  " -out " + blastdb
	print(makedbcommand)
	os.system(makedbcommand)
	blastcommand = "blastn -query " + fastafilename + " -db  " + blastdb + " -task blastn-short -outfmt 6 -out " + blastoutput + " -num_threads " + str(nproc)
	if mincoverage >=400:
		blastcommand = "blastn -query " + fastafilename + " -db " + blastdb + " -outfmt 6 -out " + blastoutput + " -num_threads " + str(nproc)
	print(blastcommand)
	os.system(blastcommand)
	if not os.path.exists(blastoutput):
		print("Cannot compare the sequences of " + fastafilename + " using Blast...")
		logfile=open(GetWorkingBase((os.path.basename(args.input))) + ".predict.log","w")
		logfile.write("Cannot compare the sequences of " + fastafilename + " using Blast...")
		logfile.write("Make BLAST database command: " + makedbcommand + "\n")
		logfile.write("No output for the BLAST command: " + blastcommand + "\n")
		logfile.write("Please rerun prediction for " + os.path.basename(fastafilename) + ".")
		logfile.close()
		return {}
	print("Reading Blast results of " + fastafilename + "...")
	simmatrix={}
	for seqid1 in seqrecords.keys():
		simmatrix.setdefault(seqid1,{})
		for seqid2 in seqrecords.keys():
			if seqid1==seqid2:
				simmatrix[seqid1][seqid2]=1
			#elif len(seqrecords.keys()) < args.maxSimMatrixSize: #load full matrix
			else:
				simmatrix[seqid1][seqid2]=0
	#read blast output
	blastoutputfile = open(blastoutput)
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
			simmatrix[i][j] = round(score, 4)
			simmatrix[j][i] = round(score, 4)
		# if len(seqrecords.keys()) < args.maxSimMatrixSize: #the full sim matrix has been loaded
		# 	if simmatrix[i][j] < score:
		# 		simmatrix[i][j]=round(score,4)
		# 		simmatrix[j][i]=round(score,4)
		# else:
		# 	if j in simmatrix[i].keys():
		# 		if simmatrix[i][j] < score:
		# 			simmatrix[i][j]=round(score,4)
		# 			simmatrix[j][i]=round(score,4)
		# 	else:
		# 		simmatrix[i][j]=round(score,4)
		# 		simmatrix[j][i]=round(score,4)
	os.system("rm " + blastoutput)
	os.system("rm " + blastdb + "*")
	#os.system("rm " + blastdb + ".*")
	return simmatrix

def LoadNeighbors(seqids,subsimmatrix,threshold):
	neighbordict={}
	for seqid in seqids:
		neighbordict.setdefault(seqid, [])
	for i in seqids:
		for j in seqids:
			if subsimmatrix[i][j] >= threshold:
				neighbordict[i].append(j)
				neighbordict[j].append(i)
	# if len(subsimmatrix.keys()) < args.maxSimMatrixSize:	 #the full matrix has been loaded
	# 	for i in seqids:
	# 		for j in seqids:
	# 			if subsimmatrix[i][j] >= threshold:
	# 				neighbordict[i].append(j)
	# 				neighbordict[j].append(i)
	# else:
	# 	for i in seqids:
	# 		for j in seqids:
	# 			if j in subsimmatrix[i].keys():
	# 				if subsimmatrix[i][j] >= threshold:
	# 					neighbordict[i].append(j)
	# 					neighbordict[j].append(i)
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

def LoadClassification(classificationfilename,rank,higherranklist):
	allclassification={}
	seqidpos,positionlist,higherpositionlist,isError = GetPositionList(classificationfilename,[rank],higherranklist)
	pos=positionlist[0]
	if isError == True:
		sys.exit()
	classificationfile = open(classificationfilename,errors='ignore')
	for line in classificationfile:
		texts = line.rstrip().split("\t")
		seqid = texts[seqidpos].rstrip()
		classname = ""
		higherclassname = ""
		if pos >-1 and pos < len(texts):
			classname = texts[pos].rstrip()
		if classname=="" or classname=="unidentified" or classname[-15:].lower()==" incertae sedis" or classname[-15:].lower()=="_incertae_sedis" or classname[-3]==" sp" or classname[-3:]=="_sp" or ("unculture" in classname):
			continue
		allclassification.setdefault(seqid,{})
		allclassification[seqid][rank]=classname
		i=0
		for higherpos in higherpositionlist:
			higherclassname=""
			if higherpos > -1 and higherpos < len(texts):
			 	higherclassname = texts[higherpos].rstrip()
			allclassification[seqid][higherranklist[i]] = higherclassname
			i = i + 1
	classificationfile.close()
	return allclassification

def LoadClasses(records,rank,allclassification):
	classification={}
	classes={}
	for seqid in records.keys():
		description=seqrecords[seqid].description
		classname=""
		if allclassification=={}:
			classname=GetTaxonName(description,rank)
		else:
			try:
				classname = allclassification[seqid][rank]
			except KeyError:
				continue
		classification[seqid]=classname
		if classname in classes.keys():
			classes[classname].append(seqid)
		else:
			classes.setdefault(classname,[seqid]) 
	return classes,classification

def LoadClassesFromClassification(records,classification):
	classes={}
	for seqid in records.keys():
		try:
			classname=classification[seqid]
		except KeyError:
			continue
		if classname in classes.keys():
			classes[classname].append(seqid)
		else:
			classes.setdefault(classname,[seqid])
	return classes

def ComputeMaxProportion(classes,seqno):
	maxproportion=0
	if seqno==0:
		for classname in classes.keys():
			seqno=seqno + len(classes[classname])
	for classname in classes.keys():
		n=len(classes[classname])
		if maxproportion < float(n/seqno):
			maxproportion= round(float(n/seqno),4)
	return maxproportion

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

def ComputeSubSim(datasetname,records,simmatrix):
#	if datasetname=="All":
#		return simmatrix
#	subsimmatrix={}
	if simmatrix!={}:
#		for seqid1 in records.keys():
#			subsimmatrix.setdefault(seqid1,{})
#			for seqid2 in records.keys():
#				subsimmatrix[seqid1][seqid2]=simmatrix[seqid1][seqid2]	
		return simmatrix
	else:
		#save sequence records to a fasta file
		subfastafilename=GetWorkingBase(datasetname) + ".fasta"
		SeqIO.write(records.values(), subfastafilename, "fasta")			
		subsimmatrix=ComputeSim(subfastafilename,records,mincoverage)
		if subsimmatrix!={}:
			os.system("rm " + subfastafilename)
	return subsimmatrix

def MergeComplexes(complexes,classes):
	removednames=[]
	complexnames=[]
	for comp in complexes:
		for classname in comp:
			if not (classname in complexnames):
				complexnames.append(classname)
	namepoints={}
	for classname in complexnames:
		namepoint = PointDef(classname,False,[],None)
		namepoints.setdefault(classname,namepoint)	
	#Load neighbor for points
	for comp in complexes:
		for classname in comp:
			for classname2 in comp:
				if not classname2 in namepoints[classname].neighbors:
					namepoints[classname].neighbors.append(classname2)			
	nameclusters=[]	
	Cluster(namepoints,nameclusters)	
	for namecluster in nameclusters:
		maxseqnumber=0
		ref=""
		for complexname in namecluster.pointids:
			seqnumber=len(classes[complexname])
			if seqnumber>=maxseqnumber:
				maxseqnumber=seqnumber
				ref=complexname
		for complexname in namecluster.pointids:
			if complexname!=ref:
				removednames.append(complexname)
	return removednames	

def RemoveComplexes(records,classification,subsimmatrix):
	distinguishablerecords={}
	newclasses={}
	#load neighbors
	neighbordict = LoadNeighbors(records.keys(),subsimmatrix,1)	
	#cluster
	points=LoadPoints(neighbordict,records)
	clusters=[]
	Cluster(points,clusters)	
	#compute complexes
	complexes=[]
	for cluster in clusters:
		speciescomplex=[]
		for seqid in cluster.pointids:
			seqrecord = records[seqid]
			classname=classification[seqid]
			if not classname in speciescomplex:
				speciescomplex.append(classname)
		complexes.append(speciescomplex)
	#merge complexes
	complexnames=MergeComplexes(complexes,classes)
	#remove complexes
	for seqid in records.keys():
		seqrecord=records[seqid]
		if not classification[seqid] in complexnames:
			distinguishablerecords.setdefault(seqid,seqrecord)
	newclasses=LoadClassesFromClassification(distinguishablerecords,classification)		
	return distinguishablerecords,newclasses

def Predict(datasetname,prediction_datasetname,records,classes,classification,simmatrix):
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
	if args.redo == "yes" and optthreshold>=t and optthreshold<=endthreshold:
		optthreshold = 0
		bestFmeasure = 0
	isError=False		
	subsimmatrix={}	
	print("Number of sequences for prediction: " + str(len(records)))
	#compute optimal threshold
	while t <= endthreshold:
		print("Computing F-measure for threshold " + str(t))
		fmeasure=0
		#if str(t) in fmeasuredict.keys() and datasetdict['fasta filename']==fastafilename and datasetdict['classification filename']==classificationfilename:
		if str(t) in fmeasuredict.keys() and args.redo=="":
			fmeasure=fmeasuredict[str(t)]
		else:
			if subsimmatrix=={}:
				#remove complexes if required	
				if args.removecomplexes=="yes":
					#compute sub simmatrix
					subsimmatrix=ComputeSubSim(datasetname,records,simmatrix)
					if subsimmatrix=={}:
						print("Cannot compute the similarity matrix for " + datasetname + ".")
						sys.exit()
					records,classes=RemoveComplexes(records,classification,subsimmatrix)
				else:	
					#compute sub simmatrix
					subsimmatrix=ComputeSubSim(datasetname,records,simmatrix)
					if subsimmatrix=={}:
						isError=True	
						break
			#compute fmeasure
			neighbordict = LoadNeighbors(records.keys(),subsimmatrix,t)	
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
	print(datasetname + ": similarity cutoff: " + str(optthreshold) + " with best F-measure: " + str(bestFmeasure))
	if isError==False:	
		prediction_datasetname['cut-off']=optthreshold
		prediction_datasetname['confidence']=bestFmeasure
		#prediction_datasetname['min alignment length']=mincoverage
		prediction_datasetname['sequence number']=len(records)
		prediction_datasetname['group number']=len(classes)
		prediction_datasetname['fmeasures']=fmeasuredict
	return thresholds,fmeasures,optthreshold,bestFmeasure,isError

def GetPositionList(classificationfilename,ranklist,higherranklist):
	positionlist=[]
	higherpositionlist=[]
	classificationfile=open(classificationfilename, errors='ignore')
	header=classificationfile.readline()
	header=header.rstrip()
	classificationfile.close()
	texts=header.split("\t")
	isError=False
	i=0
	for text in texts:
		if text.lower()==args.idcolumnname.lower():
			seqidpos=i
		i=i+1
	if 	seqidpos==-1:
		print("Please specify the sequence id columnname by using -idcolumnname.")
		isError=True
	for rank in ranklist:
		if rank in texts:
			pos=texts.index(rank)
			positionlist.append(pos)
		else:
			print("The rank " + rank + " is not given in the classification." )
			isError=True
	if isError==False:
		for rank in higherranklist:
			if rank in texts:
				pos=texts.index(rank)
				higherpositionlist.append(pos)
			else:
				print("The higher rank " + rank + " is not given in the classification." )
				isError=True
	return seqidpos,positionlist,higherpositionlist,isError

def GetTaxonName(description,rank):
	taxonname=""
	species=""
	genus=""
	family=""
	order=""
	bioclass=""
	phylum=""
	kingdom=""
	if " " in description:
		#description=description.split(" ")[1]
		description=description[description.index(" ") + 1:]	
	texts=description.split("|")
	for text in texts:
		text=text.rstrip()
		taxa=text.split(";")	
		for taxon in taxa:
			if taxon[-2:]=="__" or taxon[-15:].lower()==" incertae sedis" or taxon[-15:].lower()=="_incertae_sedis" or taxon[-3]==" sp" or taxon[-3:]=="_sp" or ("unculture" in taxon):
				continue
			if taxon.startswith("k__"):
				kingdom=taxon.replace("k__","")
			elif taxon.startswith("p__"):
				phylum=taxon.replace("p__","")
			elif taxon.startswith("c__"):
				bioclass=taxon.replace("c__","")	
			elif taxon.startswith("o__"):
				order=taxon.replace("o__","")
			elif taxon.startswith("f__"):
				family=taxon.replace("f__","")	
			elif taxon.startswith("g__"):
				genus=taxon.replace("g__","")
			elif taxon.startswith("s__") and (" " in taxon.replace("s__","") or "_" in taxon.replace("s__","")):
				species=taxon.replace("s__","")
				species=species.replace("_"," ")
	if rank.lower()=="species":
		taxonname=species
	elif rank.lower()=="genus":
		taxonname=genus
	elif rank.lower()=="family":
		taxonname=family
	elif rank.lower()=="order":
		taxonname=order
	elif rank.lower()=="class":
		taxonname=bioclass
	elif rank.lower()=="phylum":
		taxonname=phylum
	elif rank.lower()=="kingdom":
		taxonname=kingdom		
	return taxonname

def SelectList(higherclasses,maxseqno):
	##we make sure that most classes having at least one sequence selected
	selectedids = []
	for higherclassname in higherclasses.keys():
		seqids=[]
		classes=higherclasses[higherclassname]
		for classname in classes.keys():
			for seqid in classes[classname]:
				seqids.append(seqid)
		tmpids = []
		if maxseqno > 0 and maxseqno < len(seqids):
			# if maxseqno < len(classes.keys()):
			# 	selectedclassnames = random.sample(list(classes.keys()), k=maxseqno)
			# 	for classname in selectedclassnames:
			# 		tmp=random.sample(classes[classname], k=1)
			# 		tmpids.append(tmp[0])
			# else:
			# 	m=int(maxseqno/len(classes.keys()))
			# 	for classname in classes.keys():
			# 		tmp=classes[classname]
			# 		if m < len(classes[classname]):
			# 			tmp = random.sample(classes[classname], k=m)
			# 		tmpids=tmpids + tmp
			# 	notyetselected = list(set(seqids) - set(tmpids))
			# 	furtherselected = random.sample(notyetselected, k=maxseqno - len(tmpids))
			# 	tmpids = tmpids + furtherselected
			tmpids = random.sample(seqids, k=maxseqno)
			selectedids = selectedids + tmpids
		else:
			selectedids=selectedids + seqids
	return selectedids

def GenerateDatasetsForPredictionFromDescription(seqrecords,rank,higherrank,taxa,maxseqno):
	taxalist=[]
	if "," in taxa:
		taxalist=taxa.split(",")
	elif taxa!="":
		taxalist.append(taxa)
	datasets={}
	higherclasses={}
	for seqid in seqrecords.keys():
		description=seqrecords[seqid].description
		classname=GetTaxonName(description,rank)
		higherclassname=GetTaxonName(description,higherrank)
		if (len(taxalist) >0) and higherclassname!="" and not (higherclassname in taxalist) :
			continue
		if classname != "" and higherclassname !="" and classname != "unidentified" and higherclassname !="unidentified":
			if not higherclassname in higherclasses.keys():
				higherclasses.setdefault(higherclassname, {})
			if not (classname in higherclasses[higherclassname].keys()):
				higherclasses[higherclassname].setdefault(classname,[])
			higherclasses[higherclassname][classname].append(seqid)
	selectedseqids=SelectList(higherclasses,maxseqno)
	for seqid in selectedseqids:
		description = seqrecords[seqid].description
		higherclassname = GetTaxonName(description, higherrank)
		if not higherclassname in datasets.keys():
			datasets.setdefault(higherclassname, {})
		datasets[higherclassname][seqid] = seqrecords[seqid]
	return datasets

def GenerateDatasetsFromDescription(seqrecords,rank,higherranks,taxa,maxseqno):
	alldatasets={}
	if len(higherranks) !=0:		
		for higherrank in higherranks:
			datasets=GenerateDatasetsForPredictionFromDescription(seqrecords,rank,higherrank,taxa,maxseqno)
			alldatasets.update(datasets)
		return alldatasets	
	else:
		alldatasets.setdefault("All",{})
		classes={}
		for seqid in seqrecords.keys():
			description=seqrecords[seqid].description
			classname=GetTaxonName(description,rank)
			if classname != "" and classname != "unidentified":
				#alldatasets["All"][seqid]=seqrecords[seqid]
				if not (classname in classes.keys()):
					classes.setdefault(classname, [])
				classes[classname].append(seqid)
		higherclasses = {}
		higherclasses.setdefault("All", classes)
		selectedseqids=SelectList(higherclasses,maxseqno)
		for seqid in selectedseqids:
			alldatasets["All"][seqid] = seqrecords[seqid]
	return alldatasets

def GenerateDatasetsForHigherTaxa(seqrecords,allclassification,higherrank,taxa,maxseqno):
	taxalist=[]
	if "," in taxa:
		taxalist=taxa.split(",")
	elif taxa!="":
		taxalist.append(taxa)
	datasets={}
	higherclasses={}
	for seqid in seqrecords.keys():
		try:
			higherclassname=allclassification[seqid][higherrank]
		except KeyError:
			continue
		else:
			if higherclassname=="" or higherclassname=="unidentified":
				continue
			if (len(taxalist) > 0) and not (higherclassname in taxalist):
				continue
			if not higherclassname in higherclasses.keys():
				higherclasses.setdefault(higherclassname, {})
			classname = allclassification[seqid][rank]
			if not (classname in higherclasses[higherclassname].keys()):
				higherclasses[higherclassname].setdefault(classname,[])
			higherclasses[higherclassname][classname].append(seqid)
	selectedseqids=SelectList(higherclasses,maxseqno)
	for seqid in selectedseqids:
		higherclassname = allclassification[seqid][higherrank]
		if not higherclassname in datasets.keys():
			datasets.setdefault(higherclassname, {})
		datasets[higherclassname][seqid] = seqrecords[seqid]
	return datasets
	
def GenerateDatasets(seqrecords,allclassification,higherranklist,taxa,maxseqno):
	alldatasets={}
	if len(higherranklist) !=0:
		for higherrank in higherranklist:
			datasets=GenerateDatasetsForHigherTaxa(seqrecords,allclassification,higherrank,taxa,maxseqno)
			alldatasets.update(datasets)
		return alldatasets
	else:
		alldatasets.setdefault("All",{})
		classes={}
		for seqid in seqrecords.keys():
			classname=""
			try:
				classname=allclassification[seqid][rank]
			except KeyError:
				continue
			else:
				#alldatasets["All"][seqid] = seqrec
				if not (classname in classes.keys()):
					classes.setdefault(classname, [])
				classes[classname].append(seqid)
		higherclasses={}
		higherclasses.setdefault("All",classes)
		selectedseqids=SelectList(higherclasses,maxseqno)
		for seqid in selectedseqids:
			alldatasets["All"][seqid] = seqrecords[seqid]
	return alldatasets

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
			cutoff=0 
			if "cut-off" in dataset.keys():
				cutoff=dataset["cut-off"]	
			maxproportion =0
			if "max proportion" in dataset.keys():
				maxproportion=dataset["max proportion"]	
			#if groupno < minGroupNo or seqno < minSeqNo or (minalignmentlength>0 and minalignmentlength <mincoverage):	#delete the cutoffs that dont have enough sequences and groups for prediction
			if groupno < minGroupNo or seqno < minSeqNo or maxproportion > args.maxproportion or  cutoff < args.mincutoff:	#remove the cutoffs that dont have enough sequences and groups for prediction, or having a dominant group with its proportion greater than maxpropostion, or having a predicted cutoff lower than the boundary mincutoff
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
	
def PlotPrediction(datasetname,thresholdlist,fmeasurelist,optthresholds,bestFmeasures,features,datasetnames,figoutput):
	if len(thresholdlist) >5:
		fig, ax = plt.subplots(figsize=(6,3))
	else:	
		fig, ax = plt.subplots(figsize=(4,3)) 
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
	ax.set_title(datasetname + ": predicting similarity cut-offs for sequence identification")	
	plt.legend(labels, loc="lower left", bbox_transform=plt.gcf().transFigure,prop={'style': args.labelstyle})
	plt.tight_layout()
	plt.savefig(figoutput, dpi = 500)
	print("The prediction plot is saved in file " + figoutput + ".")
	if args.display=="yes":
		plt.show()	
	
def PlotResults(prefix,optthresholds,bestFmeasures,features,datasetnames,localfigoutput):
# 	#sort all according to increasing order of optthresholds
# 	if len(optthresholds)==0:
# 		return
# 	optthresholds,bestFmeasures,datasetnames,features = zip(*sorted(zip(optthresholds,bestFmeasures,datasetnames,features)))
	#sort all according to increasing order of bestFmeasures
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
	if prefix=="" or prefix==None:
		basename=os.path.basename(fastafilename)
		prefix=basename[:-(len(basename)-basename.rindex("."))] 
	outputname=GetWorkingBase(prefix) + ".predicted"	
	predictiondict={}
	if os.path.exists(outputname):
		predictiondict=LoadPrediction(outputname)
	ranklist=[]
	if classificationranks!="":
		ranklist=classificationranks.split(",")
	higherranklist=[]
	if higherclassificationranks!="":
		higherranklist=higherclassificationranks.split(",")
	if len(ranklist)==0:
		for rank in predictiondict.keys():
			prediction_datasets=predictiondict[rank]
			for datasetname in prediction_datasets.keys():
				datasetdict=prediction_datasets[datasetname]
				if "fasta filename" in datasetdict.keys() and "classification filename" in datasetdict.keys():
					if fastafilename== datasetdict["fasta filename"] and classificationfilename==datasetdict["classification filename"]:
						if not (rank in ranklist):
							ranklist.append(rank)	
	if len(ranklist)==0:
		print("Please specify the ranks for similarity cut-offs prediction by using -ranks.")	
		sys.exit()		
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
			#question=input("Do you want to load the similarity matrix from the file " + simfilename + " (yes/no)?")
			#if question=="yes":
			print("Loading similarity matrix " + simfilename)
			simmatrix=LoadSim(simfilename)
		elif higherclassificationranks=="":	#predict globally
			if len(seqrecords.keys()) <= args.maxseqno:
				print("Computing similarity matrix...")
				simmatrix=ComputeSim(fastafilename,seqrecords,mincoverage)
				print("Save similarity matrix " + simfilename)
				SaveSim(simmatrix,simfilename)
	thresholdlist=[]
	intrathresholdlist=[]
	fmeasurelist=[]
	optthresholds=[]
	bestFmeasures=[]
	features=[]
	datasetnames=[]
	#load positions of classification for prediction
	seqidpos=-1
	positionlist=[]
	higherpositionlist=[]
	isError=False
	if classificationfilename!="":
		seqidpos,positionlist,higherpositionlist,isError=GetPositionList(classificationfilename,ranklist,higherranklist)
		if isError==True:
			sys.exit()
	#predicting
	i=0
	for rank in ranklist:
		prediction_datasets={}
		if rank in predictiondict.keys():
			prediction_datasets=predictiondict[rank]
		#load datasets for prediction
		datasets={}
		allclassification={}
		if classificationfilename!="":
			#pos=positionlist[i]
			#load classification
			allclassification=LoadClassification(classificationfilename,rank,higherranklist)
			datasets=GenerateDatasets(seqrecords,allclassification,higherranklist,args.taxa,args.maxseqno)
		else:
			datasets=GenerateDatasetsFromDescription(seqrecords,rank,higherranklist,args.taxa,args.maxseqno)
		if datasets=={}:
			print("Please provide classification for the rank " + rank + ".")	
			sys.exit()	
		if step==0 or (endthreshold < threshold) or (endthreshold==0 and threshold==0):
			#load existing prediction
			for datasetname in prediction_datasets.keys():
				if not datasetname in datasets.keys():
					continue 
				datasetdict=prediction_datasets[datasetname]
				if not ("fasta filename" in datasetdict.keys() and "classification filename" in datasetdict.keys()):
					continue
				elif fastafilename != datasetdict["fasta filename"] and classificationfilename!=datasetdict["classification filename"]:
					continue
				thresholds,fmeasures,optthreshold,bestFmeasure,seqno,groupno,minalignmentlength,maxproportion=LoadPredictionForGivenRankAndDataset(datasetdict)																													   
				#if not (groupno < minGroupNo or (seqno < minSeqNo) or (minalignmentlength >0 and minalignmentlength < mincoverage)):	#for visualization				
				if not (groupno < minGroupNo or seqno < minSeqNo or maxproportion > args.maxproportion or optthreshold <args.mincutoff):	#for visualization
					if (len(higherpositionlist)==0 and datasetname == "All") or (len(higherpositionlist)!=0):	
						thresholdlist.append(thresholds)
						fmeasurelist.append(fmeasures)
						optthresholds.append(optthreshold)
						bestFmeasures.append(bestFmeasure)
						features.append(rank)
						datasetnames.append(datasetname + "(" + str(seqno) + "," + str(groupno) + ")")			
		else:
			if os.path.exists(GetWorkingBase((os.path.basename(args.input))) + ".predict.log"):		
				os.system("rm " + GetWorkingBase((os.path.basename(args.input))) + ".predict.log")
			for datasetname in datasets.keys():
				records=datasets[datasetname]
				seqno=len(records)
				#load classification at the given rank
				classes={}
				classification={}
				classes, classification=LoadClasses(records,rank,allclassification)
				maxproportion=ComputeMaxProportion(classes, seqno)
				#only predict when the numbers of the groups > 1
				if len(classes) < 2:
					continue
				#only proportion of the largest group is less than maxproportion
				if maxproportion >= args.maxproportion:
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
				thresholds,fmeasures,optthreshold,bestFmeasure,isError=Predict(datasetname,datasetdict,records,classes,classification,simmatrix)	
				if isError==False:
					datasetdict['min alignment length']=mincoverage
					datasetdict['fasta filename']=fastafilename
					datasetdict['classification filename']=classificationfilename
					datasetdict['max proportion']=maxproportion
					if not (datasetname in prediction_datasets.keys()):
						prediction_datasets[datasetname]=datasetdict		
					if not (groupno < minGroupNo or seqno < minSeqNo or maxproportion > args.maxproportion or optthreshold <args.mincutoff):	#for visualization
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
	if len(ranklist) >0:		
		if len(predictiondict.keys())>0:
			outputwithoutfmeasures=GetBase(outputname) + ".cutoffs.json"	
			SavePrediction(predictiondict,outputname,outputwithoutfmeasures)
			#SaveCutoffs(predictiondict,outputcutoffs)
			print("Only cut-offs for the clades with the numbers of sequences and subclades greater than  " + str(minSeqNo) + " and " + str(minGroupNo) + ", with the proportion of the largest group less than " + str(args.maxproportion) + ", and the predicted cutoff greater than mincutoff are saved. If you wish to save cut-offs for clades with less numbers of sequences and groups, please reset minseqno, mingroupno, maxproportion, mincutoff with -minseqno, -mingroupno,-maxproportion, and -mincutoff.")
			print("The prediction and cut-offs are saved in the files " + outputname + ", " + outputwithoutfmeasures + " and " + outputwithoutfmeasures + ".txt.")
	else:
		#load existing prediction for plotting
		for rank in predictiondict.keys():
			prediction_datasets=predictiondict[rank]
			for datasetname in prediction_datasets.keys():
				prediction_datasetname=prediction_datasets[datasetname]
				thresholds,fmeasures,optthreshold,bestFmeasure,seqno,groupno,minalignmentlength,maxproportion=LoadPredictionForGivenRankAndDataset(prediction_datasetname)
				#if not (groupno < minGroupNo or seqno < minSeqNo or (minalignmentlength >0 and minalignmentlength < mincoverage)):	#for visualization
				if not (groupno < minGroupNo or seqno < minSeqNo or maxproportion > args.maxproportion > args.maxproportion or optthreshold <args.mincutoff):	#for visualization
					thresholdlist.append(thresholds)
					fmeasurelist.append(fmeasures)
					optthresholds.append(optthreshold)
					bestFmeasures.append(bestFmeasure)
					features.append(rank)
					datasetnames.append(datasetname + "(" + str(seqno) + "," + str(groupno) + ")")
	#Plot the prediction results
	globalfigoutput=""
	localfigoutput=""
	if len(ranklist) == 1:
		globalfigoutput=GetBase(outputname) + "." + ranklist[0].replace(" ","_") + ".global.png"
		localfigoutput=GetBase(outputname) + "." + ranklist[0].replace(" ","_") + ".local.png"
	else: 	
		globalfigoutput=GetBase(outputname) + ".global.png"
		barplotfigoutput=GetBase(outputname) + ".local.png"
	if label=="":
		label=prefix	
	if len(higherranklist)==0 or len(thresholdlist)==1:
		if len(higherranklist)==0:
			#plot all predictions		
			PlotPrediction(label,thresholdlist,fmeasurelist,optthresholds,bestFmeasures,features,datasetnames,globalfigoutput)	
		elif len(thresholdlist) >0:
			#plot all predictions		
			PlotPrediction(label,thresholdlist,fmeasurelist,optthresholds,bestFmeasures,features,datasetnames,localfigoutput)	
		else:
			print("Please check the parameters.")
	else:	
		
		if len(optthresholds) >0:
			#barplot the prediction results only
			PlotResults(label,optthresholds,bestFmeasures,features,datasetnames,localfigoutput)
		else:
			print("Please check the parameters.")	
	if os.path.exists(GetWorkingBase((os.path.basename(args.input))) + ".predict.log"):		
		print("Please check the file " + GetWorkingBase((os.path.basename(args.input))) + ".predict.log for the prediction.")		
	
	
