#!/usr/bin/env python
# FILE: removesimilarsequences.py
# AUTHOR: Duong Vu
# CREATE DATE: 02 july 2020
import os
import sys, argparse
from Bio import SeqIO
import multiprocessing
nproc=multiprocessing.cpu_count()
parser=argparse.ArgumentParser(prog='removeComplexes.py',  
							   usage="%(prog)s [options] -i fastafile -t threshold -c classification -p position -out outputname",
							   description='''Script that removes similar sequences of the same complexes based on a given threshold. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )
parser.add_argument('-i','--input', required=True, help='the fasta file to be clustered.')
parser.add_argument('-t','--cutoff', type=float, default=1, help='The threshold for the classification.')
parser.add_argument('-ml','--minalignmentlength', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should probably be set to smaller, 50 for instance.')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', default="", help='the classification file in tab. format.')
#parser.add_argument('-p','--classificationpos', type=int, default=0, help='the classification position to load the classification.')
parser.add_argument('-rank','--classificationrank', default="", help='the classification rank for loading the complexes.')
parser.add_argument('-sim','--simfilename', help='The similarity matrix of the sequences if exists.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')

args=parser.parse_args()
fastafilename= args.input
threshold=args.cutoff
mincoverage = args.minalignmentlength
classificationfilename=args.classification
#classificationpos=args.classificationpos
rank=args.classificationrank
outputpath=args.out
simfilename=args.simfilename
if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))] 
	
def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + basename
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
	for line in simfile:
		numbers=line.rstrip().split(" ")
		i=numbers[0]
		j=numbers[1]
		if i not in simmatrix.keys():
			simmatrix.setdefault(i, {})
		if j not in simmatrix[i].keys():
			simmatrix[i].setdefault(j, 0)
		if float(numbers[2]) > simmatrix[i][j]:
			simmatrix[i][j]=float(numbers[2])
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
	for seqid in seqrecords.keys():
		simmatrix.setdefault(seqid,{})
		simmatrix[seqid][seqid]=1
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
		if j in simmatrix[i].keys():
			if simmatrix[i][j] < score:
				simmatrix[i][j]=round(score,4)
				simmatrix[j][i]=round(score,4)
		else:
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
			if not j in simmatrix[i].keys():
				continue
			if simmatrix[i][j] >= threshold:
				if (j not in neighbordict[i]):
					neighbordict[i].append(j)
				if (i not in neighbordict[j]):
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
	if n==0:
		return f
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
			elif taxon.startswith("s__") and (" " in taxon or "_" in taxon):
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

def LoadClasses(allseqrecords,classificationfilename,pos,seqidpos):
	classificationfile= open(classificationfilename)
#	records= open(classificationfilename,errors='ignore')
	classification={}
	classes={}
	seqrecords={}
	for line in classificationfile:
		words=line.split("\t")
		seqid= words[seqidpos].rstrip()
		if seqid in allseqrecords.keys():
			classname=""
			if pos < len(words):
				classname=words[pos].rstrip()
			if classname=="" or classname=="unidentified":	
				continue
			classification[seqid]=classname
			seqrecords[seqid]=allseqrecords[seqid]
			if classname in classes.keys():
				classes[classname].append(seqid)
			else:
				classes.setdefault(classname,[seqid]) 
	classificationfile.close()
	return seqrecords,classes,classification

def LoadClassesFromDescription(allseqrecords,rank):
	classification={}
	classes={}
	seqrecords={}
	for seqid in allseqrecords.keys():
		description=allseqrecords[seqid].description
		classname=GetTaxonName(description, rank)
		if classname=="" or classname=="unidentified":	
			continue
		classification[seqid]=classname
		seqrecords[seqid]=allseqrecords[seqid]
		if classname in classes.keys():
			classes[classname].append(seqid)
		else:
			classes.setdefault(classname,[seqid]) 
	return seqrecords,classes,classification

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

def SaveClusters(clusters,seqrecords,classes,classification,output,outputfastafilename):
	outputfile=open(output,"w")
	outputfastafile=open(outputfastafilename,"w")
	outputfile.write("ClusterID\tSequenceID\tClassification\tPrediction\n")
	complexes=[]
	numberlist=[]
	for cluster in clusters:
		speciescomplex=[]
		seqnumbers=[]
		for seqid in cluster.pointids:
			seqrecord = seqrecords[seqid]
			classname=classification[seqid]
			if classname in speciescomplex:
				i=speciescomplex.index(classname)
				seqnumbers[i]=seqnumbers[i] +1 
			else:
				speciescomplex.append(classname)
				seqnumbers.append(1)
		complexes.append(speciescomplex)
		numberlist.append(seqnumbers)			
		maxindex=seqnumbers.index(max(seqnumbers))
		classname=speciescomplex[maxindex]
		for id in cluster.pointids:
			seqrecord = seqrecords[id]
			outputfile.write(str(cluster.id) + "\t" + seqrecord.description + "\t" + classification[id] + "\t" + classname + "\n")		
	outputfile.close()
	#remove complexes
	complexnames=MergeComplexes(complexes,classes)
	for seqid in seqrecords.keys():
		seqrecord=seqrecords[seqid]
		if not classification[seqid] in complexnames:
			outputfastafile.write(">"+seqrecord.description+"\n")
			outputfastafile.write(str(seqrecord.seq)+"\n")
	outputfastafile.close()	

def GetPosition(classificationfilename,rank):
	pos=-1
	classificationfile=open(classificationfilename)
	header=classificationfile.readline()
	header=header.rstrip()
	classificationfile.close()
	texts=header.split("\t")
	isError=False
	seqidpos=-1
	i=0
	for text in texts:
		if text.lower()==args.idcolumnname.lower():
			seqidpos=i
		i=i+1
	if 	seqidpos==-1:
		print("Please specify the sequence id columnname by using -idcolumnname.")
		isError=True
	if rank in texts:
		pos=texts.index(rank)
	else:
		print("The rank " + rank + " is not given in the classification." )
		isError=True
	return seqidpos,pos,isError
	
if __name__ == "__main__":
	outputfastafilename=GetWorkingBase(fastafilename) + ".diff.fasta"	
	outputname=GetWorkingBase(fastafilename) + ".similar"
	allseqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	sys.setrecursionlimit(len(allseqrecords)*2)
	if classificationfilename!="":
		seqidpos,classificationpos,isError=GetPosition(classificationfilename,rank)
		if isError==True:
			sys.exit()
		seqrecords,classes,classification=LoadClasses(allseqrecords,classificationfilename,classificationpos,seqidpos)
	else:
		seqrecords,classes,classification=LoadClassesFromDescription(allseqrecords,rank)
	if seqrecords=={}:
		print("No classification names are available for the sequences at the rank " + rank + ".")
		os.sys.exit()
	#load similarity matrix
	if simfilename=="" or simfilename==None:
		simfilename=GetWorkingBase(fastafilename) + ".sim"
	simmatrix={}
	if os.path.exists(simfilename):
		print("Loading similarity matrix " + simfilename)
		simmatrix=LoadSim(simfilename)
	else:	
		print("Computing similarity matrix...")
		simmatrix=ComputeSim(fastafilename,seqrecords,mincoverage)
		print("Save similarity matrix " + simfilename)
		if len(seqrecords.keys())!=len(allseqrecords.keys()):
			simfilename=GetWorkingBase(simfilename) + "." + str(classificationpos) + ".sim"
		SaveSim(simmatrix,simfilename)	
	
	#load neighbors
	neighborlist = LoadNeighbors(seqrecords,simmatrix,threshold)
	print("Threshold\tFmeasure")
	points=LoadPoints(neighborlist,seqrecords)
	#cluster the sequences
	clusters=[]
	Cluster(points,clusters)
	fmeasure=ComputeFmeasure(classes,clusters)
	print(str(threshold) + "\t" + str(fmeasure))
	output=GetWorkingBase(fastafilename) + ".classified"
	SaveClusters(clusters,seqrecords,classes,classification,outputname,outputfastafilename)
	print("The remained sequences are saved in file: " + outputfastafilename )
	print("The clusters are saved in file: " + outputname )


