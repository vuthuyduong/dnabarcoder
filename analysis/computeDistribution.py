#!/usr/bin/env python
# FILE: computeDistribution.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys, argparse
if sys.version_info[0] >= 3:
	unicode = str
import numpy as np
import os
from Bio import SeqIO
#import json
#import random
import matplotlib.pyplot as plt
plt.rc('font',size=6)

import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='computeDistribution.py',  
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -ranks classificationranks  -o output",
							   description='''Script that computes the distribution of the sequences based on given classification. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be clustered.')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', default='', help='the classification file in tab. format.')
parser.add_argument('-rank','--classificationranks', default="", help='the classification ranks to compute distribution, separated by ",".')
parser.add_argument('-n','--numberofdisplayedlabels', type=int, default=5, help='The number of labels to be displayed.')
parser.add_argument('-labelstyle','--labelstyle', default='normal', help='The label style to be displayed: normal, italic, or bold.')
parser.add_argument('-method','--visualizationmethod', default="plot", help='The visualization method. There are two methods to be selected: krona and plot.')
parser.add_argument('-prefix','--prefix',default="", help='The prefix of the output files.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')

args=parser.parse_args()
referencename= args.input
classificationfilename=args.classification
jsonvariationfilename =args.out
labelno=args.numberofdisplayedlabels
method=args.visualizationmethod
outputpath=args.out
prefix=args.prefix

if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)

def GetBase(filename):
	if not ("." in filename):
		return filename
	return filename[:-(len(filename)-filename.rindex("."))] 

def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	if "." in basename:
		basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + basename
	return path

def LoadClassificationFromDescription(seqrecords,ranklist):
	classificationdict={}
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
		classification=[]
		if "kingdom" in ranklist:
			classification.append(kingdom)			
		if "phylum" in ranklist:
			classification.append(phylum)		
		if "class" in ranklist:
			classification.append(bioclass)	
		if "order" in ranklist:
			classification.append(order)
		if "family" in ranklist:
			classification.append(family)	
		if "genus" in ranklist:
			classification.append(genus)
		if "species" in ranklist:
			classification.append(species)	
		classname=""
		i=0
		for taxonname in classification:
			if taxonname=="" or ("unidentified" in taxonname):
				taxonname="unidentified"
			if i==0:
				classname=taxonname
			else:
				classname=classname + "\t" + taxonname
			i=i+1			
		if classname in classificationdict.keys():
			classificationdict[classname]=classificationdict[classname] + 1
		else:
			classificationdict.setdefault(classname, 1)
			
	return classificationdict	

def LoadClassification(seqrecords,classificationfilename,poslist,seqidpos):
	classificationdict={}
	classificationfile= open(classificationfilename)
	next(classificationfile)
	for line in classificationfile:
		elements=line.split("\t")
		seqid = elements[seqidpos].rstrip()
		if seqid in seqrecords.keys():
			classname=""
			i=0
			for pos in poslist:
				taxonname=elements[pos].rstrip()
				if taxonname=="":
					taxonname="unidentified"
				if i==0:
					classname=taxonname
				else:
					classname=classname + "\t" + taxonname
				i=i+1	
			if classname in classificationdict.keys():
				classificationdict[classname]=classificationdict[classname] + 1
			else:
				classificationdict.setdefault(classname, 1)
	return classificationdict

def SaveDistributionInTabFormat(output,classification):
	outputfile=open(output,"w")
	outputfile.write("Taxonname\tNumber of sequences\t%\n")	
	values = classification.values() 
	total= sum(values)
	sorted_classification = sorted(classification.items(), key=lambda x: x[1], reverse=True)
	i=0
	count=0
	for item in sorted_classification:
		classname=item[0]
		seqno=item[1]
		percentage=round(seqno*100/total,2)
		outputfile.write(classname + "\t" + str(seqno) + "\t" + str(percentage) + "\n")
		if i<args.numberofdisplayedlabels:
			count=count+seqno
		i=i+1
	percentage=round(count*100/total,2)
	print("The first " + str(args.numberofdisplayedlabels) + " groups of " + output + " consist of " + str(count) + "/" + str(total) + " (" +  str(percentage) + "%) sequences")	
	outputfile.close()
	
def PlotPieChart(figoutput,title,classification,displayed):
	values = classification.values() 
	total= sum(values)
	unseqno=0
	if "unidentified" in classification.keys():
		unseqno=classification["unidentified"]
		del classification['unidentified']
	#sort variations based on median thresholds with decreasing order
	sorted_classification = sorted(classification.items(), key=lambda x: x[1], reverse=True)
	classnames=[]
	seqnos=[]
	for item in sorted_classification:
		seqno=item[1]
		p=round(seqno*100/total,2)
		classname=item[0] + " " + str(p) + "%" 
		classnames.append(classname)
		seqnos.append(seqno)
	if unseqno > 0:
		p=round(unseqno*100/total,2)
		classnames.append("unidentified" + " " + str(p) + "%" )	
		seqnos.append(unseqno)
	#select labels to show
	n=len(classnames)
	if labelno > 0:
		n=min(labelno,len(classnames))
	displayedlabels=classnames[0:n]
	cmap = plt.get_cmap("tab20c")
	colors = plt.cm.Set1(np.linspace(0, 1, labelno)) #prism
	if labelno < len(classnames):
		colors1=cmap(np.arange(len(classnames)-labelno)) 
		#colors1 = plt.cm.rainbow(np.linspace(0, 1, (len(classnames)-labelno)))
		colors=np.concatenate((colors,colors1), axis=0)
	if unseqno >0:	
		colors[len(classnames)-1]= np.array([0,0,0,1])#black	
	#pie chart
	plt.close("all")
	fig, ax = plt.subplots(figsize=(3,3))
	#ax.pie(seqnos, startangle = 90, wedgeprops = { 'linewidth': 2, "edgecolor" :"k" }) #autopct='%.2f'
	ax.set_title
	ax.pie(seqnos,colors=colors,startangle = 90) #autopct='%.2f'
	#plt.title(title , horizontalalignment='center')
	plt.gca().axis("equal")
	plt.legend(displayedlabels,  loc="lower right", bbox_transform=plt.gcf().transFigure,prop={'style': args.labelstyle})#bbox_to_anchor=(1,1)
	plt.rcParams['font.size'] = 6.0
	plt.tight_layout()
	#plt.subplots_adjust(left=0.1, bottom=0.1, right=0.5)
	plt.savefig(figoutput, dpi = 500)
	if displayed==True:
		plt.show()

def PlotNestedPieCharts(figoutput,title,classificationlist,labels):
	#colors = plt.cm.Set1(np.linspace(0, 1,len(data)))	
	# create a figure with two subplots
	#pie chart
	plt.close("all")
	fig, ax = plt.subplots(figsize=(3,3))
	ax.set_title(title)
	size = 0.1
	r=1
	i=0
	sorted_classificationlist=sorted(classificationlist, key=lambda x: len(x))
	##print(sorted_classificationlist)
	displayedlabels=[]
	cmap = plt.get_cmap("tab20c")
	for classification in sorted_classificationlist:
		values = classification.values() 
		total= sum(values)
		unseqno=0
		if "unidentified" in classification.keys():
			unseqno=classification["unidentified"]
			del classification['unidentified']
		sorted_classification = sorted(classification.items(), key=lambda x: x[1], reverse=True)
		classnames=[]
		seqnos=[]
		for item in sorted_classification:
			seqno=item[1]	
			p=round(seqno*100/total,2)
			classname=item[0] + " " + str(p) + "%" 
			classnames.append(classname)
			seqnos.append(seqno)
		if unseqno > 0:
			p=round(unseqno*100/total,2)
			classnames.append("unidentified" + " " + str(p) + "%" )	
			seqnos.append(unseqno)	
		#pie chart
		#select labels to show
		if i==0:
		#select labels to show
			n=len(classnames)
			if labelno > 0:
				n=min(labelno,len(classnames))
				displayedlabels=classnames[0:n]
		#plt.pie(seqnos, labels = displayedlabels, startangle = 90)
		#colors=cmap(np.arange(len(classnames))) 
		colors = plt.cm.Set1(np.linspace(0, 1, labelno)) #prism
		if labelno < len(classnames):
			colors1=cmap(np.arange(len(classnames)-labelno)) 
			#colors1 = plt.cm.rainbow(np.linspace(0, 1, (len(classnames)-labelno)))
			colors=np.concatenate((colors,colors1), axis=0)
		if unseqno >0:	
			colors[len(classnames)-1]= np.array([0,0,0,1])#black			
		ax.pie(seqnos, radius=r, startangle = 90, colors=colors, wedgeprops=dict(width=size, edgecolor='w'))
		r=r-size
		i=i+1	
	#plt.title(title, horizontalalignment='center')
	plt.gca().axis("equal")
	plt.legend(displayedlabels,  loc="lower right", bbox_transform=plt.gcf().transFigure,prop={'style': args.labelstyle})#bbox_to_anchor=(1,1)
	plt.rcParams['font.size'] = 6.0
	plt.tight_layout()
	plt.savefig(figoutput, dpi = 500)
	plt.show()
	
def KronaPieCharts(classification,kronareport,kronahtml):
	kronareportfile=open(kronareport,"w")
	for classname in classification.keys():
		kronareportfile.write(str(classification[classname]) + "\t" + classname + "\n")
	kronareportfile.close()	
	#create kronahtml
	command="ImportText.pl " + kronareport + " -o " + kronahtml
	print(command)
	os.system(command)
	os.system("firefox " + kronahtml) 
	
def GetPositionList(classificationfilename,ranklist):	
	positionlist=[]
	classificationfile=open(classificationfilename)
	header=classificationfile.readline()
	header=header.rstrip()
	classificationfile.close()
	texts=header.rstrip().split("\t")
	isError=False
	seqidpos=-1
	i=0
	for text in texts:
		if text.lower()==args.idcolumnname.lower():
			seqidpos=i
		i=i+1
	if 	seqidpos==-1:
		print("Please specify ID columnname by using --idcolumnname.")
		isError=True
	for rank in ranklist:
		if rank in texts:
			pos=texts.index(rank)
			positionlist.append(pos)
		else:
			print("The rank " + rank + " is not given in the classification." )
			isError=True
	return positionlist,seqidpos,isError
##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]
displayed=True
poslist=[]

if prefix=="":
	prefix=GetBase(os.path.basename(referencename))
			
ranklist=[]
if "," in args.classificationranks:
	ranklist=args.classificationranks.split(",")
else:
	ranklist=[args.classificationranks]	

poslist=[]
seqidpos=0
if classificationfilename!="":
	poslist,seqidpos,isError=GetPositionList(classificationfilename,ranklist)
	if isError==True:
		sys.exit()		

displayed=False
if len(ranklist)==1:
	displayed=True
	
#load train seq records
referencerecords = SeqIO.to_dict(SeqIO.parse(referencename, "fasta"))
if method!="krona":
	classificationlist=[]
	labels=[]
	for rank in ranklist:
		jsonfilename = GetWorkingBase(prefix) + "." + rank + ".distribution"
		figoutput=GetBase(jsonfilename) + ".distribution.png" 
		#Load classes, classification:
		classificationdict={}
		if classificationfilename!="":
			pos=poslist[ranklist.index(rank)]
			classificationdict=LoadClassification(referencerecords,classificationfilename,[pos],seqidpos)
		else:
			classificationdict=LoadClassificationFromDescription(referencerecords,[rank])
		title=""
		if rank.lower()== "species":
			title=prefix + ": the distribution of the sequences at the species level"		
		elif rank.lower()== "genus":
			title=prefix + ": the distribution of the sequences at the genus level"	
		elif rank.lower()== "family":
			title=prefix + ": the distribution of the sequences at the family level"	
		elif rank.lower()== "order":
			title=prefix + ": the distribution of the sequences at the order level"	
		elif rank.lower()== "class":
			title=prefix + ": the distribution of the sequences at the class level"	
		elif rank.lower()== "phylum":
			title=prefix + ": the distribution of the sequences at the phylum level"	
		elif rank.lower()== "kingdom":
			title=prefix + ": the distribution of the sequences at the kingdom level"	
		else:
			title= prefix + ": the distribution of the the sequences of the groups at the columnname " + rank
		#save classification	
		SaveDistributionInTabFormat(jsonfilename + ".txt",classificationdict)	
		newclassification=classificationdict.copy()
		classificationlist.append(newclassification)
		labels.append(rank)
		#plot
		PlotPieChart(figoutput,title,classificationdict,displayed)
		print("The results are saved in the json file  " + jsonfilename + " and tab file " + jsonfilename + ".txt. The figure is saved in " + figoutput + "."  )	
	if len(ranklist)>1:
		jsonfilename=""
		jsonfilename = GetWorkingBase(prefix) + ".distribution"
		figoutput=jsonfilename + ".png" 
		print("The figure of the variations of all groups are saved in file " + figoutput + ".")
		title=prefix + ": the distribution of the sequences"
		PlotNestedPieCharts(figoutput,title,classificationlist,labels)
else:
	kronareport = GetWorkingBase(prefix) + ".krona.report"
	kronahtml=GetBase(kronareport) + ".html"
	classificationdict={}
	if classificationfilename!="":
		classificationdict=LoadClassification(referencerecords,classificationfilename,poslist,seqidpos)
	else:
		classificationdict=LoadClassificationFromDescription(referencerecords,ranklist)
	KronaPieCharts(classificationdict,kronareport,kronahtml)
	print("The krona report and html are saved in iles " + kronareport + " and " + kronahtml + ".") 
		
	
			


