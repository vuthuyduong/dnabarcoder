#!/usr/bin/env python
# FILE: blastlocally.py, resulted as similarity matrix
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import os
import multiprocessing
from Bio import SeqIO


queryname=sys.argv[1] #folder containing the classifier name and config file
refname=sys.argv[2]
mincoverage=100
if len(sys.argv) >3:
	mincoverage=int(sys.argv[3]) #300 for ITS sequences
nproc=multiprocessing.cpu_count()
if len(sys.argv) >4:
	nproc=int(sys.argv[4])
shortsequences=0
if len(sys.argv) >5:
	shortsequences=int(sys.argv[5])
outputname=""
if len(sys.argv) >6:
	outputname=sys.argv[6]

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))] 

def ComputeBLASTscore(queryname,refname,mincoverage):
	queryrecords = list(SeqIO.parse(queryname, "fasta"))
	refrecords = list(SeqIO.parse(refname, "fasta"))
	scorematrix = [[0 for x in range(len(refrecords))] for y in range(len(queryrecords))] 
	simmatrix = [[0 for x in range(len(refrecords))] for y in range(len(queryrecords))] 
	coveragematrix = [[0 for x in range(len(refrecords))] for y in range(len(queryrecords))] 
	queryname= queryname
	refname=  refname
	queryids =[]
	for rec in queryrecords:
		queryids.append(rec.id)
	refids =[]
	for rec in refrecords:
		refids.append(rec.id)
	
	#blast
	makedbcommand = "makeblastdb -in " + refname + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + queryname + " -db  db -outfmt 6 -out blastout.txt -num_threads " + str(nproc)
	if shortsequences==1:
		blastcommand = "blastn -query " + queryname + " -db db -task blastn-short -outfmt 6 -out blastout.txt -num_threads " + str(nproc)
	
	os.system(blastcommand)
	#read blast output
	blastoutputfile = open("blastout.txt")
	refid = ""
	score=0
	queryid=""
	for line in blastoutputfile:
		if line.rstrip()=="":
			continue
		words = line.split("\t")
		queryid = words[0].rstrip()
		refid = words[1].rstrip()
		i = queryids.index(queryid)
		j = refids.index(refid)
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
			simmatrix[i][j]=sim
			simmatrix[j][i]=sim
			coveragematrix[i][j]= coverage
			coveragematrix[j][i]= coverage
	os.system("rm blastout.txt")
	return queryrecords,refrecords,scorematrix,simmatrix,coveragematrix


if __name__ == "__main__":
	if outputname=="":
		outputname=GetBase(queryname) + ".locally.blastresult"
	output=open(outputname,"w")
	queryrecords,refrecords,scorematrix,localsimmatrix,localcoveragematrix=ComputeBLASTscore(queryname,refname,mincoverage)
	for i in range(0,len(queryrecords)):
		for j in range(0,len(refrecords)):
			output.write(queryrecords[i].id + "\t" + refrecords[j].id + "\t" + str(scorematrix[i][j]) + "\t" + str(localsimmatrix[i][j]) + "\t" + str(localcoveragematrix[i][j]) + "\n")
			print(queryrecords[i].id + "\t" + refrecords[j].id + "\t" + str(scorematrix[i][j]) + "\t"  + str(localsimmatrix[i][j]) + "\t" + str(localcoveragematrix[i][j]) )
		
	output.close()
