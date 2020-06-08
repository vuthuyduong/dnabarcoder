#!/usr/bin/env python
# FILE: clusterlocally.py
# AUTHOR: Duong Vu
# CREATE DATE: 08 oct 2019
import os
import sys
from Bio import SeqIO
import multiprocessing

fastafilename = sys.argv[1] # fasta file containing DNA sequences
threshold = 0.97
if len(sys.argv)>2:
	threshold = float(sys.argv[2])
mincoverage = 100
if len(sys.argv)>3:
	mincoverage=int(sys.argv[3])
classificationfilename=""
if len(sys.argv) >4:	
	classificationfilename = sys.argv[4] # the file containing taxonomic classification: class order family genus species
classificationpos=0
if len(sys.argv) >5:
	classificationpos = int(sys.argv[5]) # the position of the feature in the header of the sequences used for comparisons
outputname= ""
if len(sys.argv) > 6:
	outputname  =sys.argv[6]

nproc=multiprocessing.cpu_count()

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))] 

class PointDef:
	def __init__(self, id,name,flag,neighbors,seqrecord):
		self.id=id
		self.name=name
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

def ComputeBLASTscore(fastafilename,seqrecords,mincoverage):
	scorematrix = [[0 for x in range(len(seqrecords))] for y in range(len(seqrecords))] 
	seqids=[]
	for rec in seqrecords:
		seqids.append(rec.id)
	#blast
	makedbcommand = "makeblastdb -in " + fastafilename + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + fastafilename + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	os.system(blastcommand)
	
	#read blast output
	blastoutputfile = open("out.txt")
	refid = ""
	score=0
	queryid=""
	for line in blastoutputfile:
		if line.rstrip()=="":
			continue
		words = line.split("\t")
		queryid = words[0].rstrip()
		refid = words[1].rstrip()
		i = seqids.index(queryid)
		j = seqids.index(refid)
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		score=sim
		if coverage < mincoverage:
			score=float(score * coverage)/mincoverage
		if scorematrix[i][j] < score:
			scorematrix[i][j]=score
			scorematrix[j][i]=score
	os.system("rm out.txt")
	return scorematrix

def LoadNeighbors(seqrecords,scorematrix,threshold):
	neighborlist=[]
	for i in range(len(seqrecords)):
		neighborlist.append([])
	for i in range(0,len(seqrecords)-1):
		for j in range(i+1,len(seqrecords)):
				if scorematrix[i][j] >= threshold:
					if (j not in neighborlist[i]):
						neighborlist[i].append(j)
					if (i not in neighborlist[j]):
						neighborlist[j].append(i)
	#os.system("rm out.txt")
	return neighborlist

def LoadPoints(neigborlist,seqrecords):
	points=[]
	i=0
	for seqrecord in  seqrecords:
		point = PointDef(i,seqrecord.id,False,neigborlist[i],seqrecord)
		points.append(point)
		i=i+1
	return points

def ExpandCluster(root,cluster,points):
	for i in root.neighbors:
		if points[i].flag==False:
			points[i].flag=True
			cluster.pointids.append(i)
			ExpandCluster(points[i],cluster,points)		

def Cluster(points,clusters):
	for point in points:
		if point.flag==False:
			point.flag=True
			cluster = ClusterDef(len(clusters),[])
			cluster.pointids.append(point.id)
			ExpandCluster(point, cluster,points)
			clusters.append(cluster)
		
def ComputeFmeasure(classes,clusters):
	#compute F-measure
	f=0
	n=0
	for group in classes:
		m = 0
		for cluster in clusters:
			i = len(set(group) & set(cluster.pointids))
			v = float(2*i)/float((len(group) + len(cluster.pointids)))
			if m < v:
				m=v		
		n = n + len(group)
		f = f +	(len(group)*m)	
	return float(f)/float(n) 

def LoadClasses(seqrecords,classificationfilename,pos):
	records= open(classificationfilename)
#	records= open(classificationfilename,errors='ignore')
	allseqids=[]
	allclassification=[]
	for line in records:
		words=line.split("\t")
		seqid = words[0].rstrip().replace(".1","").replace(">","")
		allseqids.append(seqid)
		allclassification.append(words[pos].rstrip())
	classes = []
	classnames = []
	classification=[]
	i=0
	for seqrecord in seqrecords:
		classname=""
		if seqrecord.id in allseqids:
			index=allseqids.index(seqrecord.id)
			classname=allclassification[index]
		if classname in classnames:
			classid=classnames.index(classname)
			classes[classid].append(i)
		else:
			classnames.append(classname)
			refclass=[]
			refclass.append(i)
			classes.append(refclass)
		classification.append(classname)
		i=i+1
	return classes,classnames,classification

def SaveClusters(clusters,seqrecords,classification,output):
	outputfile=open(output,"w")
	outputfile.write("ClusterID\tSequenceID\tClassification\tPrediction\n")
	for cluster in clusters:
		classnames=[]
		seqnumbers=[]
		for id in cluster.pointids:
			seqrecord = seqrecords[id]
			classname=classification[id]
			if classname in classnames:
				i=classnames.index(classname)
				seqnumbers[i]=seqnumbers[i] +1 
			else:
				classnames.append(classname)
				seqnumbers.append(1)
		maxindex=seqnumbers.index(max(seqnumbers))
		classname=classnames[maxindex]
		for id in cluster.pointids:
			seqrecord = seqrecords[id]
			outputfile.write(str(cluster.id) + "\t" + seqrecord.description + "\t" + classification[id] + "\t" + classname + "\n")		
	outputfile.close()

if __name__ == "__main__":
	if outputname=="":
		outputname=GetBase(fastafilename) + ".locally.clustered"
	seqrecords=list(SeqIO.parse(fastafilename, "fasta"))
	sys.setrecursionlimit(len(seqrecords)*2)
	#seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	classes,classnames,classification=LoadClasses(seqrecords,classificationfilename,classificationpos)
	#load similarity matrix
	scorematrix=ComputeBLASTscore(fastafilename,seqrecords,mincoverage)
	#load neighbors
	neighborlist = LoadNeighbors(seqrecords,scorematrix,threshold)
	print("Threshold\tFmeasure")
	points=LoadPoints(neighborlist,seqrecords)
	clusters=[]
	Cluster(points,clusters)
	fmeasure=ComputeFmeasure(classes,clusters)
	print(str(threshold) + "\t" + str(fmeasure))
	output=GetBase(fastafilename) + ".classified"
	SaveClusters(clusters,seqrecords,classification,outputname)


