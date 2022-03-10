#!/usr/bin/env python
# FILE: dnabarcodes_visualize.py
# AUTHOR: Duong Vu
# CREATE DATE: 20 oct 2020

import sys, argparse
import os
if sys.version_info[0] >= 3:
	unicode = str
import numpy
import matplotlib.pyplot as plt
plt.rc('font',size=6)
from mpl_toolkits import mplot3d	
from Bio import SeqIO
import json
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='visualize.py',  
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -rank classificationprank",
							   description='''Script that visualizes the dna sequences based on a given classification. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out',default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', help='the classification file in tab. format.')
#parser.add_argument('-p','--classificationpos', type=int, default=0, help='the classification positions for the prediction.')
parser.add_argument('-rank','--classificationrank', default="", help='the classification rank for coloring the sequences.')
parser.add_argument('-sim','--simfilename', help='The similarity file if exists, othter it will be computed.')
parser.add_argument('-ml','--minalignmentlength', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should be set to smaller, 50 for instance.')
parser.add_argument('-coord','--coordinates', help='A file containing coordinates of the sequences, computed by LargeVis. If these coordinates will be computed if this file is not given.')
parser.add_argument('-dim','--dimension', type=int, default=3,help='The dimension 2D or 3D for visualization.')
parser.add_argument('-kneigh','--kneigbors', type=int, default=150,help='The k-neighbors number for visualization.')
parser.add_argument('-ms','--minsim', type=float, default=0, help='The minimum similarity for visualization.')
parser.add_argument('-method','--visualizationmethod', default="", help='The visualization method. There are two methods to be selected: dive and plot.')
#parser.add_argument('-lim','--lim', type=int, default=50, help='The lim for visualization.')
parser.add_argument('-size','--size', type=float, default=1, help='The size of the dot.')
parser.add_argument('-n','--numberofdisplayedlabels', type=int, default=5, help='The size of the dot.')
parser.add_argument('-prefix','--prefix',default="", help='The prefix of the output files.')
parser.add_argument('-label','--label',default="", help='The label to display in the figure.')

args=parser.parse_args()
fastafilename= args.input
simfilename=args.simfilename
coordfilename=args.coordinates
classificationfilename=args.classification
#classificationpos=args.classificationpos
rank=args.classificationrank
mincoverage=args.minalignmentlength
minsim=args.minsim
dim=args.dimension
kneigh=args.kneigbors
method=args.visualizationmethod
size=args.size
numberofdisplayedlabels=args.numberofdisplayedlabels
outputpath=args.out
prefix=args.prefix
label=args.label

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

def LoadSim(simfilename,n):
	simmatrix = {} #we use dictionary to reduce the memory constraints 
	simfile = open(simfilename)
	for line in simfile:
		texts=line.rstrip().split(" ")
		i=texts[0]
		j=texts[1]
		if i not in simmatrix.keys():
			simmatrix.setdefault(i, {})
		if j not in simmatrix[i].keys():
			simmatrix[i].setdefault(j, 0)
		if float(texts[2]) > simmatrix[i][j]:
			simmatrix[i][j]=float(texts[2])
	simfile.close()		
	return simmatrix

def SaveSim(simmatrix,simfilename,ms):
	simfile=open(simfilename,"w")
	for i in simmatrix.keys():
		for j in simmatrix[i].keys():
			if simmatrix[i][j] >= ms:
				simfile.write(str(i) + " " + str(j) + " " + str(simmatrix[i][j]) + "\n")
	simfile.close()

def ComputeSim(fastafilename,seqids,mincoverage,minsim):
	#scorematrix = [[0 for x in range(len(seqrecords))] for y in range(len(seqrecords))] 
	simmatrix={}
	for seqid in seqids:
		simmatrix.setdefault(seqid,{})
		simmatrix[seqid][seqid]=1
#	simfilename=GetBase(fastafilename) + ".sim"
#	simfile=open(simfilename,"w")
	#blast
	makedbcommand = "makeblastdb -in " + fastafilename + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + fastafilename + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	if mincoverage >=300:
		blastcommand = "blastn -query " + fastafilename + " -db  db -outfmt 6 -out out.txt -num_threads " + str(nproc)
	os.system(blastcommand)
	
	#read blast output
	blastoutputfile = open("out.txt")
	for line in blastoutputfile:
		if line.rstrip()=="":
			continue
		words = line.split("\t")
#		queryid = words[0].rstrip()
#		refid = words[1].rstrip()
		i = words[0].rstrip()
		j = words[1].rstrip()
#		i = seqids.index(queryid)
#		j = seqids.index(refid)
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		score=sim
		if coverage < mincoverage:
			score=float(score * coverage)/mincoverage
		if j in simmatrix[i].keys():
			if simmatrix[i][j] < score:
				simmatrix[i][j]=round(score,4)
				simmatrix[j][i]=round(score,4)
		else:
			simmatrix[i][j]=round(score,4)
			simmatrix[j][i]=round(score,4)
#		if scorematrix[i][j] < score:
#			scorematrix[i][j]=score
#			scorematrix[j][i]=score
	os.system("rm out.txt")
	return simmatrix

def ComputeCoordinates(base,simfilename,dim,edgeNo,kneigh):
	coordfilename==GetBase(fastafilename) + ".coord"
	largeViscommand = "LargeVis -fea 0 -input " + simfilename + " -output " + coordfilename + " -outdim " + str(dim) +  " -threads "  + str(nproc) + " -log 1 -samples " + str(edgeNo) + " -neigh " + str(kneigh);
	print(largeViscommand)
	os.system(largeViscommand)
	return coordfilename

def LoadFullClassification(seqids,classificationfilename):
	features=""
	classification=[""]*len(seqids)
	if classificationfilename!=None:
		classificationfile=open(classificationfilename)
		firstline=classificationfile.readline()
		features=firstline.replace("\n","").split("\t")
		if not firstline.startswith("#"):
			seqid = firstline.split("\t")[0].replace(">","")
			if seqid in seqids:  
				classification[seqids.index(seqid)]=firstline.replace("\n","")	        
		for line in classificationfile:
			seqid = line.split("\t")[0].replace(">","")
			if seqid in seqids:
				classification[seqids.index(seqid)]=line.replace("\n","").split("\t")
		classificationfile.close()				
	else:
		features=["seqindex"]
		for i in range(0,len(seqids)):
			classification[i]=[str(i)]
	return features,classification

def LoadClassification(seqids,classificationfilename,classificationpos):
	labels=[""]*len(seqids)
	rank=""
	if classificationfilename!=None:
		classificationfile=open(classificationfilename)
		i=0
		for line in classificationfile:
			seqid = line.split("\t")[0].replace(">","")
			if i==0:
				rank=line.replace("\n","").split("\t")[classificationpos]
			if seqid in seqids:
				labels[seqids.index(seqid)]=line.replace("\n","").split("\t")[classificationpos]
			i=i+1	
		classificationfile.close()				
	else:
		for i in range(0,len(seqids)):
			labels[i]=str(i)
	return labels, rank

def LoadCoordinates(coorfilename):
	coordfile=open(coordfilename)
	firstline=coordfile.readline().rstrip()
	seqNo=int(firstline.split(" ")[0])
	dim=int(firstline.split(" ")[1])
	coordinates={}
	for line in coordfile:
		line=line.rstrip()
		seqid=line.split(" ")[0]
		line=line[line.index(" ")+1:]
		coordinate=line.split(" ")
		coordinates[seqid]=coordinate
	coordfile.close()	
	return seqNo,dim,coordinates

def CreateDiVEFiles(base,seqids,classificationfilename,coorfilename):
	jsonfilename=GetWorkingBase(fastafilename) + ".dive.json"
	divefilename=base + "DiVE/data/data.js"
	features,classification=LoadFullClassification(seqids,classificationfilename)
	seqNo,dim,coordinates=LoadCoordinates(coorfilename)
	coorddict={}
	coorddict["NamesOfProperties"]=features
	i=0
	for seqid in seqids:
		item={}
		item["Coordinates"]=coordinates[seqid]
		item["Properties"]=classification[i]
		coorddict[seqid]=item
		#coorddict[str(i)]=item
		i=i+1
	#write to file
	with open(jsonfilename,"w") as json_file:
		if sys.version_info[0] >= 3:
			json.dump(coorddict,json_file,indent=2)
		else:
			json.dump(coorddict,json_file,encoding='latin1',indent=2)
	f = open(jsonfilename, "r")
	contents = f.readlines()
	f.close()
	contents.insert(0, "const data_all = ")
	f = open(divefilename, "w")
	contents = "".join(contents)
	f.write(contents)
	f.close()
	return divefilename,jsonfilename

def VisualizeUsingDiVE(base,seqids,classificationfilename,coorfilename):
	divefilename,jsonfilename=CreateDiVEFiles(base,seqids,classificationfilename,coorfilename)
	print("The DiVE file is created as " + jsonfilename + ".")
	os.system("firefox " + base + "DiVE/index.html");
	
def Plot(prefix,seqids,coordfilename,classificationfilename,classificationpos,size,output):
	labels,rank = LoadClassification(seqids,classificationfilename,classificationpos)
	seqNo,dim,coordinates=LoadCoordinates(coordfilename)
	all_data = {}
	cmap = plt.get_cmap("tab20c")
	for i in range(0,seqNo):
		all_data.setdefault(labels[i], [0,[]])
	i=0	
	for seqid in seqids:
		coord=coordinates[seqid]
		all_data[labels[i]][0]=all_data[labels[i]][0] + 1
		if dim==2:
			#all_data.setdefault(labels[i], []).append((float(coord[-2]), float(coord[-1])))
			all_data[labels[i]][1].append((float(coord[-2]), float(coord[-1])))
		else:
			all_data[labels[i]][1].append((float(coord[-3]),float(coord[-2]),float(coord[-1])))
		i=i+1	
	#sort all_data based on number of points with a decreasing order
	sorted_data = sorted(all_data.items(), key=lambda x: x[1], reverse=True)
	#colors = plt.cm.rainbow(numpy.linspace(0, 1, len(all_data)))
	n1=min(10,len(all_data))
	colors = plt.cm.tab10(numpy.linspace(0, 1, n1)) #prism
#	if len(colors) >=5:
#		colors=numpy.delete(colors,5,0)#delete the yellow color from the list
#		if len(colors) >=7:
#			colors=numpy.delete(colors,7,0)#delete the yellow color from the list		
#	n1=len(colors)
#	print(n1)
	if n1 < len(all_data):
		#colors1 = plt.cm.rainbow(numpy.linspace(0, 1, (len(all_data)-n1)))
		colors1=cmap(numpy.arange(len(all_data)-n1)) 
		colors=numpy.concatenate((colors,colors1), axis=0)
	 #Set1, primse 
	fig, ax = plt.subplots(figsize=(4,3)) 
	if prefix=="":
		ax.set_title("Sequence distribution \n at the " + rank + " level.",loc='left')
	else:
		ax.set_title(prefix + ": sequence distribution \n at the " + rank + " level.",loc='left')
	labels=[]
	if dim==2:
		# Hide grid lines
		#ax.grid(False)
		# Hide axes ticks
		#ax.set_xticks([])
		#ax.set_yticks([])
		for color, item in zip(colors, sorted_data):
			x = [t[0] for t in item[1][1]]
			y = [t[1] for t in item[1][1]]
			taxon=item[0]
			if taxon=="":
				taxon="unidentified"
			if taxon=="unidentified":
				color="black"	
			label=taxon + " " + str(item[1][0])
			ax.plot(x, y, '.', color = color, markersize = size)	
			labels.append(label)	
	elif dim==3:
		ax = plt.axes(projection='3d')
		ax.set_title(prefix + ": visualization of the sequences \n at the " + rank + " level.", loc='left')
		# Hide grid lines
		#ax.grid(False)
		# Hide axes ticks
		#ax.set_xticks([])
		#ax.set_yticks([])
		#ax.set_zticks([])
		for color, item in zip(colors, sorted_data):
			x = [t[0] for t in item[1][1]]
			y = [t[1] for t in item[1][1]]
			z = [t[2] for t in item[1][1]]
			taxon=item[0]
			if taxon=="":
				taxon="unidentified"
			if taxon=="unidentified":
				color="black"
			label=taxon + " " + str(item[1][0])
			ax.plot3D(x, y, z, '.', color = color, markersize = size)
			labels.append(label)	
			#ax.scatter3D(x, y, z, '.', color = color)
			#ax.axis('off')
	leg=fig.legend(labels[0:numberofdisplayedlabels],loc='best')
	#leg=fig.legend(labels[0:numberofdisplayedlabels],loc='lower left')
	i=0
	for text in leg.get_texts():
		color= colors[i]
		label=text.get_text()
		if label.startswith("unidentified"):
			color="black"
		plt.setp(text, color = color)
		i=i+1
	#plt.xlim(-lim, lim)
    #plt.ylim(-lim, lim)
	plt.rcParams['font.size'] = 6.0
	plt.savefig(output, dpi = 500)
	plt.show()		
	
def GetPosition(classificationfilename,rank):
	pos=-1
	classificationfile=open(classificationfilename)
	header=classificationfile.readline()
	header=header.rstrip()
	classificationfile.close()
	texts=header.split("\t")
	isError=False
	if rank in texts:
		pos=texts.index(rank)
	else:
		print("The rank " + rank + " is not given in the classification." )
		isError=True
	return pos,isError
	
if __name__ == "__main__":
	#load sequences
	if prefix=="":
		prefix=GetBase(os.path.basename(fastafilename))
	base = sys.argv[0][0: sys.argv[0].rindex("/")+1]
	seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
	seqids=[]
	for rec in seqrecords:
		seqids.append(rec.id)
	if coordfilename==None:
		if simfilename==None or simfilename=="":			
			simfilename=GetWorkingBase(prefix)  + ".sim"
			simfilename_minsim=""
			if minsim >0:
				simfilename_minsim=GetWorkingBase(prefix)  + "." + str(int(round(minsim*100,0))) + ".sim"
			else:
				simfilename_minsim=simfilename
		else:
			simfilename_minsim=simfilename
		coordfilename=GetBase(simfilename_minsim)  + "." + str(dim) + ".coord"	
		if not os.path.isfile(coordfilename):
			if not os.path.isfile(simfilename_minsim):
				if not os.path.isfile(simfilename):
					#compute similarity matrix
					print("Computing similarity matrix..")
					simmatrix=ComputeSim(fastafilename,seqids,mincoverage,minsim)
					print("Saving similarity matrix " + simfilename )
					SaveSim(simmatrix,simfilename,0)
					if minsim >0:
						print("Saving similarity matrix " + simfilename_minsim )
						SaveSim(simmatrix,simfilename_minsim,minsim)
				else:		
					print("Loading the exising similarities " + simfilename)
					simmatrix=LoadSim(simfilename,minsim)
					if minsim>0:
						print("Saving similarity matrix " + simfilename_minsim )
						SaveSim(simmatrix,simfilename_minsim,minsim)
			else:
				#use the existing simfilename.
				print("The existing file " + simfilename_minsim + " is used for loading similarity matrix. If you wish to recalculate the similarity matrix, please delete the existing file.")
			#compute coordinates	
			edgeNo=int(len(seqrecords)/100)
			if edgeNo <50:
				edgeNo=50			
			coordfilename=ComputeCoordinates(base,simfilename_minsim,dim,edgeNo,kneigh)
			print("The coordinates are saved in file " + coordfilename)
		else:
			#use the existing simfilename.
			print("The existing file " + coordfilename + " is used for visualization. If you wish to recalculate the coordinates, please delete the existing file.")
	if method.lower()=="dive":
		VisualizeUsingDiVE(base,seqids,classificationfilename,coordfilename)		
	else:
		classificationpos,isError=GetPosition(classificationfilename,rank)
		output=GetWorkingBase(prefix)  + "." + str(classificationpos) + ".visualization.png"
		if label=="":
			label=prefix
		Plot(label,seqids,coordfilename,classificationfilename,classificationpos,size,output)
		print("The visualization is saved in " + output + ".")
			

