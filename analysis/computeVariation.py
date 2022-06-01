#!/usr/bin/env python
# FILE: computeVariation.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys, argparse
if sys.version_info[0] >= 3:
	unicode = str
import numpy as np
import os
from Bio import SeqIO
import json
import random
import matplotlib.pyplot as plt
plt.rc('font',size=6)
from matplotlib.patches import Polygon
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='computeVariation.py',  
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -ranks classificationranks -ml minalignmentlength  -o output",
							   description='''Script that computes the median and minimum similarity scores within the groups. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be clustered.')
parser.add_argument('-ml','--minalignmentlength', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should be set to smaller, 50 for instance.')
parser.add_argument('-o','--out',default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', default="", help='the classification file in tab. format.')
parser.add_argument('-rank','--classificationranks', default="", help='the classification ranks to compute variation, separated by ",".')
parser.add_argument('-m','--maxSeqNo', type=int, default=0, help='The maximum number of randomly selected sequences of each class to be computed in the case the groups are too big.')
parser.add_argument('-plt','--plottype', default="boxplot", help='The type of plots. There are two options: boxplot and plot.')
parser.add_argument('-sim','--simfilename', default="", help='The similarity matrix of the sequences if exists.')
parser.add_argument('-prefix','--prefix',default="", help='The prefix of the output files.')
parser.add_argument('-label','--label',default="", help='The label to display in the figure.')
parser.add_argument('-maxSimMatrixSize','--maxSimMatrixSize', type=int, default=20000, help='The maximum number of sequences to load or compute a full similarity matrix. In case the number of sequences is greater than this number, only similarity values greater than 0 will be loaded to avoid memory problems.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')


args=parser.parse_args()
referencename= args.input
mincoverage = args.minalignmentlength
classificationfilename=args.classification
jsonvariationfilename =args.out
plottype=args.plottype
simfilename=args.simfilename
prefix=args.prefix
label=args.label


maxSeqNo=0
if args.maxSeqNo !=None:
	maxSeqNo=args.maxSeqNo
outputpath=args.out
if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)	

def GetBase(filename):
	if not ("." in filename):
		return filename
	return filename[:-(len(filename)-filename.rindex("."))] 

#def GetWorkingBase(filename):
#	basename=os.path.basename(filename)
#	if "." in basename:
#		basename=basename[:-(len(basename)-basename.rindex("."))] 
#	path=outputpath + "/" + basename
#	return path

def GetWorkingBase(basename):
	path=outputpath + "/" + basename
	return path

def LoadClassificationFromDescription(seqrecords,rank):
	for seqid in seqrecords.keys():
		description=seqrecords[seqid].description
		species=""
		genus=""
		family=""
		order=""
		bioclass=""
		phylum=""
		kingdom=""
		if " " in description:
			description=description.split(" ")[1]
		texts=description.split("|")
		for text in texts:
			text=text.rstrip()
			taxa=text.split(";")	
			for taxon in taxa:
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
			classname=species
		elif rank.lower()=="genus":
			classname=genus
		elif rank.lower()=="family":
			classname=family
		elif rank.lower()=="order":
			classname=order
		elif rank.lower()=="class":
			classname=bioclass
		elif rank.lower()=="phylum":
			classname=phylum
		elif rank.lower()=="kingdom":
			classname=kingdom
		if classname=="" or ("unidentified" in classname):
			continue 
		if seqid in seqrecords.keys():
			if not (classname in classes.keys()):
				classes.setdefault(classname,[])	
			classes[classname].append(seqrecords[seqid])				
	return classes

def LoadClassification(seqrecords,classificationfilename,pos,seqidpos):
	classes={}
	if classificationfilename == "":
		return classes
	records= open(classificationfilename)
	next(records)
	for line in records:
		elements=line.split("\t")
		seqid = elements[seqidpos].replace(">","").rstrip()
		classname=""
		if pos < len(elements):
			 classname=elements[pos].rstrip()
		if classname=="" or classname=="unidentified":
			continue 
		if seqid in seqrecords.keys():
			if not (classname in classes.keys()):
				classes.setdefault(classname,[])	
			classes[classname].append(seqrecords[seqid])
	records.close()			
	return classes

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
		if len(seqids) < args.maxSimMatrixSize: #load full matrix	
			for seqid2 in seqids:
				if not seqid2 in simmatrix[seqid1].keys():
					#simmatrix[seqid1].setdefault(seqid2,0)
					simmatrix[seqid1][seqid2]=0
	simfile.close()		
	return simmatrix

def GetSeqIndex(seqname,seqrecords):
	i=0
	for seqrecord in seqrecords:
		if (seqname == seqrecord.id):
			return i
		i = i + 1
	return -1

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
			elif len(seqrecords.keys()) < args.maxSimMatrixSize: #load full matrix	
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
		if len(seqrecords.keys()) < args.maxSimMatrixSize: #the full sim matrix has been loaded 
			if simmatrix[i][j] < score:
				simmatrix[i][j]=round(score,4)
				simmatrix[j][i]=round(score,4)
		else:		
			if j in simmatrix[i].keys():
				if simmatrix[i][j] < score:
					simmatrix[i][j]=round(score,4)
					simmatrix[j][i]=round(score,4)
			else:
				simmatrix[i][j]=round(score,4)
				simmatrix[j][i]=round(score,4)	
		#simmatrix[j][i]=score
	os.system("rm " + blastoutput)
	os.system("rm " + blastdb + "*")
	#os.system("rm " + blastdb + ".*")
	return simmatrix

def ComputeVariation(reffilename,mincoverage,simmatrix):
	#load sequeces from the fasta files
	records = SeqIO.to_dict(SeqIO.parse(reffilename, "fasta"))
	keys=list(records.keys())
	scorelist=[]
	check=False
	if simmatrix!={}:
		check=True
		for i in range(0,len(keys)-2):
			if not keys[i] in simmatrix.keys():
				check==False
				break
			for j in range(i+1,len(keys)-1):
				if i!=j:
					if not keys[j] in simmatrix[keys[i]].keys():
						check==False
						break
					scorelist.append(simmatrix[keys[i]][keys[j]])	
	if check==False:
		scorematrix=ComputeSim(reffilename,records,mincoverage)
		for i in range(0,len(keys)-2):
			for j in range(i+1,len(keys)-1):
				if i!=j:
					scorelist.append(scorematrix[keys[i]][keys[j]])	
	threshold=1
	minthreshold=1		
	if len(scorelist) >0:
		x=np.array(scorelist)
		minthreshold=round(float(np.min(x)),4)
		threshold=round(float(np.median(x)),4)
	return threshold,minthreshold

def ComputeVariations(variationfilename,classes,mincoverage,simmatrix):
	#create json dict
	variations={}
	i=0
	for taxonname in classes.keys():
		threshold=0
		minthreshold=0
		sequences=classes[taxonname]
		if len(sequences) >0:
			if maxSeqNo==0 or (len(sequences) < maxSeqNo):
				fastafilename=taxonname.replace(" ","_") + ".fasta"
				SeqIO.write(sequences,fastafilename,"fasta")
				threshold,minthreshold=ComputeVariation(fastafilename,mincoverage,simmatrix)
				os.system("rm " + fastafilename)
			else:
				threshold,minthreshold=EvaluateVariation(taxonname,sequences,mincoverage)
			currentvariation=[threshold,minthreshold,len(sequences)]
			variations[taxonname]=currentvariation
		i=i+1	
	#write to file
	with open(variationfilename,"w") as json_file:
		if sys.version_info[0] >= 3:
			json.dump(variations,json_file,indent=2)	
		else:
			json.dump(variations,json_file,encoding='latin1',indent=2)	
	return variations

def EvaluateVariation(taxonname,sequences,mincoverage,simmatrix):
	selectedindexes=random.sample(range(0, len(sequences)), k=maxSeqNo)
	selectedsequences=[]
	for index in selectedindexes:
		selectedsequences.append(sequences[index])
	fastafilename=taxonname.replace(" ","_") + ".fasta"
	SeqIO.write(selectedsequences,fastafilename,"fasta")
	threshold,minthreshold=ComputeVariation(fastafilename,mincoverage,simmatrix)
	return threshold,minthreshold
		
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

def SaveVariationInTabFormat(output,variation):
	outputfile=open(output,"w")
	outputfile.write("Taxonname\tMedian similarity score\tMin similarity score\tNumber of sequences\n")	
	for classname in variation.keys():
		threshold=variation[classname][0]
		minthreshold=variation[classname][1]
		seqno=variation[classname][2]
		outputfile.write(classname + "\t" + str(threshold) + "\t" + str(minthreshold) + "\t" + str(seqno) + "\n")
	outputfile.close()
	
def Plot(datasetname,figoutput,variations,rank,displayed):
	#sort variations based on median thresholds with decreasing order
	sorted_variations = sorted(variations.items(), key=lambda x: x[1][0], reverse=True)
	thresholds=[]
	minthresholds=[]
	seqnos=[]
	for item in sorted_variations:
		threshold=item[1][0]
		minthreshold=item[1][1]
		seqno=item[1][2]
		thresholds.append(threshold)
		minthresholds.append(minthreshold)
		seqnos.append(seqno)
	median_median=round(np.median(np.array(thresholds))	,4)
	median_min=round(np.median(np.array(minthresholds))	,4)
	x = np.arange(len(seqnos))
	fig, ax2 = plt.subplots(figsize=(3,3))
	ax = ax2.twinx()  # instantiate a second axes that shares the same x-axi
	ax2.set_ylabel('Number of sequences')  # we already handled the x-label with ax1
	ax2.plot(x, np.array(seqnos), color='g')
	ax.set_title(datasetname + ": median and minimum similarity scores of the " + rank.lower() )
	ax2.set_xlabel("Group index")
	ax.set_ylabel('Similarity score')
	#plt.plot(x, np.array(thresholds), 'r--', x,minthresholds, 'bs') #green 'g^'
	ax.plot(x, np.array(thresholds), 'b--', label=  'Median. Median score' + str(median_median))
	ax.plot(x, np.array(minthresholds), 'rs', label='Min. Median score: ' + str(median_min))
	plt.legend()
	plt.tight_layout()
	plt.rcParams['font.size'] = 6.0
	plt.savefig(figoutput, dpi = 500)
	if displayed==True:
		plt.show()

def PlotAll(datasetname,figoutput,variationlist,labels):
	data=[]
	for variations in variationlist:
		sorted_variations = sorted(variations.items(), key=lambda x: x[1][0], reverse=True)
		thresholds=[]
		for item in sorted_variations:
			threshold=item[1][0]
			#minthreshold=item[1][1]
			#seqno=item[1][2]
			thresholds.append(threshold)
			#minthresholds.append(minthreshold)
			#seqnos.append(seqno)
		data.append(thresholds)
	colors = plt.cm.Set1(np.linspace(0, 1,len(data)))	
	fig, ax = plt.subplots(figsize=(3,3))
	if datasetname=="":
		ax.set_title("Median similarity scores of all groups")
	else:
		ax.set_title(datasetname + ": median similarity scores of all groups")
	ax.set_xlabel("Group index")
	ax.set_ylabel('Median similarity score')
	#plt.plot(x, np.array(thresholds), 'r--', x,minthresholds, 'bs') #green 'g^'
	k=0
	for thresholds in data:
		x = np.arange(len(thresholds))
		median_median=np.median(np.array(thresholds))
		ax.plot(x, np.array(thresholds), color=colors[k], label=  labels[k] + '. Median ' + str(median_median))
		k=k+1
	plt.legend()
	plt.tight_layout()
	plt.rcParams['font.size'] = 6.0
	plt.savefig(figoutput, dpi = 500)
	plt.show()
	
def BoxPlot(datasetname,figoutput,variations,rank,displayed):
	#sort variations based on median thresholds with decreasing order
	sorted_variations = sorted(variations.items(), key=lambda x: x[1][0], reverse=True)
	thresholds=[]
	minthresholds=[]
	seqnos=[]
	for item in sorted_variations:
		threshold=item[1][0]
		minthreshold=item[1][1]
		seqno=item[1][2]
		thresholds.append(threshold)
		minthresholds.append(minthreshold)
		seqnos.append(seqno)
	labels=["Median", "Min"]	
	data=[np.array(thresholds),np.array(minthresholds)]
#	fig, ax = plt.subplots(figsize=(10, 6))
#	#fig.canvas.set_window_title('Variation')
#	fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
	fig, ax = plt.subplots(figsize=(3,3))
	box_colors = ['b','r']#['darkkhaki', 'royalblue']
	bp = ax.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
	plt.setp(bp['boxes'], color='black')
	plt.setp(bp['whiskers'], color='black')
	plt.setp(bp['fliers'], color='red', marker='+')
	
	ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)

	# Hide these grid behind plot objects
	ax.set_axisbelow(True)
	if datasetname=="":
		ax.set_title('Median and min. similarity scores of ' + rank)
	else:
		ax.set_title(datasetname + ': median and min. similarity scores of ' + rank)
	#ax.set_xlabel('')
	ax.set_ylabel('Similarity score')
	num_boxes=len(data)
	medians=np.empty(num_boxes)
	for i in range(num_boxes):
		box=bp['boxes'][i]
		boxX=[]
		boxY=[]
		for j in range(5):
			boxX.append(box.get_xdata()[j])
			boxY.append(box.get_ydata()[j])
		box_coords = np.column_stack([boxX, boxY])	
		# Fill in the color
		ax.add_patch(Polygon(box_coords, facecolor=box_colors[i % 2]))
		med = bp['medians'][i]
		medianX=[]
		medianY=[]
		for j in range(2):
			medianX.append(med.get_xdata()[j])
			medianY.append(med.get_ydata()[j])
		medians[i] = medianY[0]
		#plot the average value
		ax.plot(np.average(med.get_xdata()), np.average(data[i]),color='w', marker='*', markeredgecolor='k')
	#add labels	
	ax.set_xticklabels(np.array(labels))	
	#add median values 
	upper_labels = [str(np.round(s, 4)) for s in medians]
	pos = np.arange(num_boxes) + 1
	k=0
	for tick, label in zip(range(num_boxes), ax.get_xticklabels()):
		ax.text(pos[tick], 0.97, upper_labels[tick], transform=ax.get_xaxis_transform(), horizontalalignment='center', size='x-small', color=box_colors[k])
		k=k+1
	#plt.legend()
	plt.tight_layout()
	plt.rcParams['font.size'] = 6.0
	plt.savefig(figoutput, dpi = 500)
	if displayed==True:
		plt.show()
		
def BoxPlotAll(datasetname,figoutput,variationlist,labels):
	data=[]
	labels2=[]
	colors=[]
	i=0
	for variations in variationlist:
		sorted_variations = sorted(variations.items(), key=lambda x: x[1][0], reverse=True)
		thresholds=[]
		minthresholds=[]
		for item in sorted_variations:
			threshold=item[1][0]
			minthreshold=item[1][1]
			#seqno=item[1][2]
			thresholds.append(threshold)
			minthresholds.append(minthreshold)
			#seqnos.append(seqno)
		data.append(thresholds)
		data.append(minthresholds)
		labels2.append("Median_" + labels[i])
		labels2.append("Min_" + labels[i])
		colors.append('b')
		colors.append('r')
		i=i+1
#	fig, ax = plt.subplots(figsize=(10, 6))
#	#fig.canvas.set_window_title('Variation')
#	fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
	#box_colors = ['r','b']#['darkkhaki', 'royalblue']
	fig, ax = plt.subplots(figsize=(4,3))
	bp = ax.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
	plt.setp(bp['boxes'], color='black')
	plt.setp(bp['whiskers'], color='black')
	plt.setp(bp['fliers'], color='red', marker='+')
	
	ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
	# Hide these grid behind plot objects
	ax.set_axisbelow(True)
	ax.set_title(datasetname + ': median and min. similarity scores of all groups')
	#ax.set_xlabel('')
	ax.set_ylabel('Similarity score')
	num_boxes=len(data)
	#colors = plt.cm.Set1(np.linspace(0, 1,num_boxes))
	medians=np.empty(num_boxes)
	for i in range(num_boxes):
		box=bp['boxes'][i]
		boxX=[]
		boxY=[]
		for j in range(5):
			boxX.append(box.get_xdata()[j])
			boxY.append(box.get_ydata()[j])
		box_coords = np.column_stack([boxX, boxY])	
		# Fill in the color
		ax.add_patch(Polygon(box_coords, facecolor=colors[i]))
		med = bp['medians'][i]
		medianX=[]
		medianY=[]
		for j in range(2):
			medianX.append(med.get_xdata()[j])
			medianY.append(med.get_ydata()[j])
		medians[i] = medianY[0]
		#plot the average value
		ax.plot(np.average(med.get_xdata()), np.average(data[i]),color='w', marker='*', markeredgecolor='k')	
	#add labels	
	ax.set_xticklabels(np.array(labels2), rotation=90)	
	#add median values 
	upper_labels = [str(np.round(s, 4)) for s in medians]
	pos = np.arange(num_boxes) + 1
	k=0
	for tick, label in zip(range(num_boxes), ax.get_xticklabels()):
		ax.text(pos[tick], 0.97, upper_labels[tick], transform=ax.get_xaxis_transform(), horizontalalignment='center', size='x-small', color=colors[k])
		k=k+1
	#plt.legend()
	plt.tight_layout()
	plt.rcParams['font.size'] = 6.0
	plt.savefig(figoutput, dpi = 500)
	plt.show()		

def GetPositionList(classificationfilename,ranks):
	ranklist=[]	
	if "," in ranks:
		ranklist=ranks.split(",")
	elif ranks !="":
		ranklist.append(ranks)
	positionlist=[]	
	isError=False	
	seqidpos=-1	
	classificationfile=open(classificationfilename)
	header=classificationfile.readline()
	header=header.rstrip()
	classificationfile.close()
	texts=header.rstrip().split("\t")
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
	return seqidpos,positionlist,ranklist,isError

##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]

if prefix=="":
	prefix=GetBase(os.path.basename(referencename))

seqidpos,positionlist,ranklist,isError=GetPositionList(classificationfilename,args.classificationranks)	
if isError==True :
	sys.exit()
displayed=False
if len(positionlist)==1:
	displayed=True
#load similarity matrix
simmatrix={}	
if os.path.exists(simfilename):
	print("Loading similarity matrix " + simfilename)
	simmatrix=LoadSim(simfilename)	
#load reference seq records
referencerecords =  SeqIO.to_dict(SeqIO.parse(referencename, "fasta"))
variationlist=[]
labels=[]
i=0
jsonvariationfilename=""
figoutput=""
for rank in ranklist:
	rank=rank.lower()
	jsonvariationfilename = GetWorkingBase(prefix) + "." + rank + ".variation"
	figoutput=GetBase(jsonvariationfilename) + ".variation.png" 
	#Load classes, classification:
	classes={}
	if classificationfilename !="":
		classificationposition=positionlist[i]
		classes=LoadClassification(referencerecords,classificationfilename, classificationposition,seqidpos)
	else:
		classes=LoadClassificationFromDescription(referencerecords,rank)
	variations={}
	if not os.path.exists(jsonvariationfilename):
		variations=ComputeVariations(jsonvariationfilename,classes,mincoverage,simmatrix)
	else:
		print("The variation file " + jsonvariationfilename + " exists. Please delete the file if you wish to recalculate the variation.")
		with open(jsonvariationfilename) as variation_file:
			variations = json.load(variation_file)
		SaveVariationInTabFormat(jsonvariationfilename + ".txt",variations)
		print("The variations are saved in the json file  " + jsonvariationfilename + " and tab file " + jsonvariationfilename + ".txt. The figure is saved in " + figoutput + "."  )
	variationlist.append(variations)
	labels.append(rank)	
	i=i+1	
if label=="":
	label=prefix	
if len(positionlist)>1:
	jsonvariationfilename = GetWorkingBase(prefix) + ".variation"
	figoutput=jsonvariationfilename + ".png" 
if plottype=="plot":
	PlotAll(label,figoutput,variationlist,labels)
else:	
	BoxPlotAll(label,figoutput,variationlist,labels)
print("All variations and theirs figure are saved in file " + jsonvariationfilename + " and " + figoutput + ".")
			


