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
import json
import random
import matplotlib.pyplot as plt
plt.rc('font',size=6)
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
import numpy as np
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='computeDistribution.py',  
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -p classificationposition  -o output",
							   description='''Script that computes the distribution of the sequences based on given classification. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be clustered.')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', help='the classification file in tab. format.')
parser.add_argument('-p','--classificationpos', help='the classification positions for computing sequence variation.')
parser.add_argument('-n','--numberofdisplayedlabels', type=int, default=5, help='The number of labels to be displayed.')
parser.add_argument('-method','--visualizationmethod', default="plot", help='The visualization method. There are two methods to be selected: krona and plot.')


args=parser.parse_args()
referencename= args.input
classificationfilename=args.classification
jsonvariationfilename =args.out
labelno=args.numberofdisplayedlabels
method=args.visualizationmethod
outputpath=args.out
if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + basename
	return path

def LoadClassification(seqIDs,seqrecords,classificationfilename,poslist):
	classification={}
	if classificationfilename == "":
		classification["unidentified"]=len(seqrecords)
		return classification
	records= list(open(classificationfilename, "r"))
	feature=""
	for line in records:
		if line.startswith("#"):
			i=0 
			for pos in poslist:
				if i==0:
					feature=line.rstrip().split("\t")[pos]
				else:
					feature=feature + line.rstrip().split("\t")[pos]
				i=i+1
			continue 
		elements=line.split("\t")
		seqid = elements[0].replace(">","").rstrip()
		if seqid in seqIDs:
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
			if classname in classification.keys():
				classification[classname]=classification[classname] + 1
			else:
				classification.setdefault(classname, 1)
	return classification,feature

def SaveDistributionInTabFormat(output,classification):
	outputfile=open(output,"w")
	outputfile.write("Taxonname\tNumber of sequences\t%\n")	
	values = classification.values() 
	total= sum(values)
	sorted_classification = sorted(classification.items(), key=lambda x: x[1], reverse=True)
	for item in sorted_classification:
		classname=item[0]
		seqno=item[1]
		percentage=round(seqno*100/total,2)
		outputfile.write(classname + "\t" + str(seqno) + "\t" + str(percentage) + "\n")
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
	fig, ax = plt.subplots(figsize=(3,3))
	#ax.pie(seqnos, startangle = 90, wedgeprops = { 'linewidth': 2, "edgecolor" :"k" }) #autopct='%.2f'
	#ax.set_title
	ax.pie(seqnos,colors=colors,startangle = 90) #autopct='%.2f'
	#plt.title(title , horizontalalignment='center')
	plt.gca().axis("equal")
	plt.legend(displayedlabels,  loc="lower right", bbox_transform=plt.gcf().transFigure)#bbox_to_anchor=(1,1)
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
	fig, ax = plt.subplots(figsize=(3,3))
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
	plt.legend(displayedlabels,  loc="lower right", bbox_transform=plt.gcf().transFigure)#bbox_to_anchor=(1,1)
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
##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]
displayed=True
poslist=[]

if args.classificationpos==None or args.classificationpos=="":
	classificationfile=open(classificationfilename)
	firstline=classificationfile.readline()
	classificationfile.close()
	featureNo=len(firstline.split("\t"))
	for i in range(1,featureNo):
		poslist.append(i)
	displayed=False
elif "," in args.classificationpos:
	texts=args.classificationpos.split(",")
	for t in texts:
		poslist.append(int(t))
	displayed=False	
else:
	poslist.append(int(args.classificationpos))
	
#load train seq records
referencerecords =  list(SeqIO.parse(referencename, "fasta"))
referenceIDs=[]
for seq in referencerecords:
	referenceIDs.append(seq.id)
if method!="krona":
	classificationlist=[]
	labels=[]
	for classificationposition in poslist:
		jsonfilename = GetWorkingBase(referencename) + "." + str(classificationposition) + ".distribution"
		figoutput=GetBase(jsonfilename) + ".distribution.png" 
		#Load classes, classification:
		classification,feature= LoadClassification(referenceIDs,referencerecords,classificationfilename, [classificationposition])
		title=""
		if feature=="":
			title="The distribution of the the sequences of the groups at position " + str(classificationposition)
		elif feature.lower()== "species":
			title="The distribution of the sequences at the species level"		
		elif feature.lower()== "genus":
			title="The distribution of the sequences at the genus level"	
		elif feature.lower()== "family":
			title="The distribution of the sequences at the family level"	
		elif feature.lower()== "order":
			title="The distribution of the sequences at the order level"	
		elif feature.lower()== "class":
			title="The distribution of the sequences at the class level"	
		elif feature.lower()== "phylum":
			title="The distribution of the sequences at the phylum level"	
		elif feature.lower()== "kingdom":
			title="The distribution of the sequences at the kingdom level"	
		#save classification	
		SaveDistributionInTabFormat(jsonfilename + ".txt",classification)	
		newclassification=classification.copy()
		classificationlist.append(newclassification)
		labels.append(feature)
		#plot
		PlotPieChart(figoutput,title,classification,displayed)
		print("The results are saved in the json file  " + jsonfilename + " and tab file " + jsonfilename + ".txt. The figure is saved in " + figoutput + "."  )	
	if len(poslist)>1:
		jsonfilename=""
		jsonfilename = GetWorkingBase(referencename) + ".distribution"
		figoutput=jsonfilename + ".png" 
		print("The figure of the variations of all groups are saved in file " + figoutput + ".")
		title="The distribution of the sequences"
		PlotNestedPieCharts(figoutput,title,classificationlist,labels)
else:
	kronareport = GetWorkingBase(referencename) + ".krona.report"
	kronahtml=GetBase(kronareport) + ".html"
	classification,feature= LoadClassification(referenceIDs,referencerecords,classificationfilename, poslist)
	KronaPieCharts(classification,kronareport,kronahtml)
	print("The krona report and html are saved in iles " + kronareport + " and " + kronahtml + ".") 
		
	
			


