#!/usr/bin/env python
# FILE: predictlocally.py
# AUTHOR: Duong Vu
# CREATE DATE: 08 oct 2019
import os
import sys
from Bio import SeqIO
import multiprocessing

fastafilename = sys.argv[1] # fasta file containing DNA sequences
classificationfilename = sys.argv[2] # the file containing taxonomic classification: class order family genus species
classificationpos = sys.argv[3] # the position of the feature in the header of the sequences used for comparisons
threshold = 0.97
if len(sys.argv)>4:
	threshold = float(sys.argv[4])
endthreshold = 1
if len(sys.argv)>5:
	endthreshold = float(sys.argv[5])
step = 0.001
if len(sys.argv)>6:
	step = float(sys.argv[6])
mincoverage = 300
if len(sys.argv)>7:
	mincoverage=int(sys.argv[7])	
outputname= ""
if len(sys.argv) > 8:
	outputname  =sys.argv[8]

nproc=multiprocessing.cpu_count()

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))] 

def FilterSequencesWithoutClassification(fastafilename,classificationfilename,positionlist):
	#load classification
	allseqids=[]
	records= open(classificationfilename)
	rank=str(positionlist[0])
	for record in records:
		texts=record.split("\t")
		if record.startswith("#"):
			for pos in positionlist:					        
				if pos < len(texts):
					rank=rank + texts[pos].rstrip() + "."
			continue 		
		seqid=texts[0].replace(">","").replace(".1","").rstrip()
		isclassified=True
		for pos in positionlist:
			if pos < len(texts):
				if texts[pos].rstrip()=="":   
					isclassified=False
		if isclassified == True:
			allseqids.append(seqid)
	records.close()
	
	fastafile=open(fastafilename)
	newfastafilename=GetBase(fastafilename) + "." + rank + "fasta"
	newfastafile=open(newfastafilename,"w")
	writetofile=False
	for line in fastafile:
		if line.startswith(">"):
			writetofile=False
			seqid=line.split("|")[0].replace(">","").rstrip()	
			if seqid in allseqids:
				#write to file
				newfastafile.write(">" + seqid + "\n")
				writetofile=True
		else:
			if writetofile==True:
				newfastafile.write(line)
	fastafile.close()
	newfastafile.close()
	return newfastafilename, rank

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
			#print("cluster id:" + str(cluster.id) + "\t" + str(len(cluster.pointids)))
			ExpandCluster(points[i],cluster,points)		

def Cluster(points,clusters):
	for point in points:
		if point.flag==False:
			point.flag=True
			cluster = ClusterDef(len(clusters),[])
			cluster.pointids.append(point.id)
			#print("cluster id:" + str(cluster.id) + "\t" + str(len(cluster.pointids)))
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
	records= open(classificationfilename,errors='ignore')
	allseqids=[]
	allclassification=[]
	for line in records:
		words=line.split("\t")
		seqid = words[0].rstrip().replace(".1","").replace(">","")
		if pos < len(words):
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

if __name__ == "__main__":
	positionlist=[]
	if "," in classificationpos:
		positions=classificationpos.split(",")
		for pos in positions:
			positionlist.append(int(pos))
	else:
		positionlist.append(int(classificationpos))
	newfastafilename,rank=FilterSequencesWithoutClassification(fastafilename,classificationfilename,positionlist)
	if outputname=="":
		outputname=GetBase(fastafilename) + "." + rank + "locally.predicted"
	#load sequences
	seqrecords=list(SeqIO.parse(newfastafilename, "fasta"))
	print(len(seqrecords))
	sys.setrecursionlimit(len(seqrecords)*2)
	#seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	#load similarity matrix
	scorematrix=ComputeBLASTscore(newfastafilename,seqrecords,mincoverage)	
	for pos in positionlist:
		#load classification at given positions 
		classes,classnames,classification=LoadClasses(seqrecords,classificationfilename,pos)
		if not os.path.isfile(outputname):
			outputfile=open(outputname,"w")
		else:
			outputfile=open(outputname,"a")
		outputfile.write("Threshold\tFmeasure at the classification position " + str(pos) + " \n")
		outputfile.close()
		t=threshold
		optthreshold=threshold
		bestFmeasure=0
		#compute optimal threshold
		while t < endthreshold:
			#load neighbors
			neighborlist = LoadNeighbors(seqrecords,scorematrix,t)	
			points=LoadPoints(neighborlist,seqrecords)
			clusters=[]
			Cluster(points,clusters)	
			fmeasure=ComputeFmeasure(classes,clusters)
			if fmeasure > bestFmeasure:
				bestFmeasure=fmeasure
				optthreshold=t
			#print(str(t) + "\t" + str(fmeasure))
			outputfile=open(outputname,"a")
			outputfile.write(str(t) + "\t" + str(fmeasure)+"\n")
			outputfile.close()
			t=t+step
		outputfile=open(outputname,"a")	
		outputfile.write("Optimal threshold at the classification position " + str(pos) + ": " + str(optthreshold) + "\tFmeasure: " +   str(bestFmeasure) + "\n")
		outputfile.close()
		print("Optimal threshold at position " + str(pos) + "\tFmeasure")
		print(str(optthreshold) + "\t" + str(bestFmeasure))
	
	
	

