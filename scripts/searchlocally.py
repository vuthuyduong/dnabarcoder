#!/usr/bin/env python
# FILE: searchlocally.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019

import os
import sys
from Bio import SeqIO
import multiprocessing

query = sys.argv[1] # fasta file containing genes/proteins sequences top50_UNITE.fas
reference = sys.argv[2] # fasta file containing reference sequenc  ITSBarcodesforblast.fas
mincoverage=100
if len(sys.argv) >3:
	mincoverage = int(sys.argv[3]) # the minimum coverage of the sequences in comparisons, for ITS, it is between 300-600
outputname = ""
if len(sys.argv) >4:
	outputname = sys.argv[4]

nproc=multiprocessing.cpu_count()

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

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


def GetSeqIndex(seqid,seqrecords):
	i=0
	for seqrecord in seqrecords:
		if (seqid == seqrecord.id):
			return i
		i = i + 1
	return -1

def ComputeBLASTscore(query,reference,mincoverage):
	indexed_query= IndexSequences(query)

	#load sequeces from the fasta files
	queryrecords = list(SeqIO.parse(indexed_query, "fasta"))
	#refrecords = list(SeqIO.parse(reference, "fasta"))

	bestscorelist =[0] * len(queryrecords)
	bestsimlist = [0] * len(queryrecords)
	bestoverlaplist = [0] * len(queryrecords)
	bestrefidlist = [""] * len(queryrecords)
	bestminposlist = [0] * len(queryrecords)
	bestmaxposlist =[0] * len(queryrecords)
	bestrefminposlist = [0] * len(queryrecords)
	bestrefmaxposlist = [0] * len(queryrecords)

	#blast
	makedbcommand = "makeblastdb -in " + reference + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + indexed_query + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	#blastcommand = "blastn -query " + indexed_query + " -subject " + reference + " -outfmt 6 -out out.txt"
	os.system(blastcommand)
	
	#read blast output
	blastoutputfile = open("out.txt")
	refid = ""
	queryid=""
	for line in blastoutputfile:
		words = line.split("\t")
		queryid = words[0]
		refid = words[1]
		pos1 = int(words[6])
		pos2 = int(words[7])
		refpos1 =  int(words[8])
		refpos2 = int(words[9])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		minpos=min(pos1,pos2)
		maxpos=max(pos1,pos2)
		refminpos=min(refpos1,refpos2)
		refmaxpos=max(refpos1,refpos2)
		i = int(queryid.split("|")[0])	
		score=sim
		if coverage < mincoverage:
			score=float(score * coverage) /mincoverage
		if score > bestscorelist[i]:
			bestscorelist[i]= score
			bestsimlist[i]= sim
			bestoverlaplist[i]= coverage
			bestminposlist[i]=minpos
			bestmaxposlist[i]=maxpos
			bestrefminposlist[i]=refminpos
			bestrefmaxposlist[i]=refmaxpos
			bestrefidlist[i]=refid			
	os.system("rm " + indexed_query)		
	os.system("rm out.txt")
	return queryrecords,bestrefidlist,bestscorelist,bestsimlist,bestoverlaplist,bestminposlist,bestmaxposlist,bestrefminposlist,bestrefmaxposlist
#MAIN
if __name__ == "__main__":
	if outputname=="":
		outputname=GetBase(query) + ".locally.searched"
	queryrecords,bestrefidlist,bestscorelist,bestsimlist,bestoverlaplist,bestminposlist,bestmaxposlist,bestrefminposlist,bestrefmaxposlist=ComputeBLASTscore(query,reference,mincoverage)
	outputfile = open(outputname,"w") 
	i=0
	outputfile.write("Sequence id\tReference sequence id\tBest score\tIdentity\tCoverage\tBest local sim\tBest local coverage\tQuery first position \tQuery last position\tRef first position\tRef last position\n")
	for seqrecord in queryrecords:
		description=seqrecord.description
		description=description[description.index("|")+1:]
		outputfile.write(description + "\t"+ bestrefidlist[i] + "\t" + str(bestscorelist[i]) + "\t" + str(bestsimlist[i]) + "\t" + str(bestoverlaplist[i]) + "\t" + str(bestminposlist[i]) + "\t" + str(bestmaxposlist[i]) + "\t" + str(bestrefminposlist[i]) + "\t" + str(bestrefmaxposlist[i]) + "\n")
		i=i+1
	outputfile.close()
