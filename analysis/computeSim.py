#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 12:28:09 2021

@author: Duong Vu
"""
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
import multiprocessing
from Bio import SeqIO

nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='computeSim.py',  
							   usage="%(prog)s [options] -i fastafile -o output",
							   description='''Script that computes BLAST similarity matrix for sequences of a fasta file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out',default="dnabarcoder", help='The output folder.') 
parser.add_argument('-ml','--minalignmentlength', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should be set to smaller, 50 for instance.')
parser.add_argument('-ms','--minsim', type=float, default=0, help='The minimum similarity score that will be saved for the output.')

args=parser.parse_args()
fastafilename= args.input
mincoverage=args.minalignmentlength
minsim=args.minsim
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

def SaveSim(simmatrix,simfilename,ms):
	simfile=open(simfilename,"w")
	for i in simmatrix.keys():
		for j in simmatrix[i].keys():
			if simmatrix[i][j] >= ms:
				simfile.write(str(i) + " " + str(j) + " " + str(simmatrix[i][j]) + "\n")
	simfile.close()

def ComputeSim(fastafilename,mincoverage,minsim):
	#scorematrix = [[0 for x in range(len(seqrecords))] for y in range(len(seqrecords))] 
#	simfilename=GetBase(fastafilename) + ".sim"
#	simfile=open(simfilename,"w")
	#blast
	makedbcommand = "makeblastdb -in " + fastafilename + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + fastafilename + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	if mincoverage >=400:
		blastcommand = "blastn -query " + fastafilename + " -db  db -outfmt 6 -out out.txt -num_threads " + str(nproc)
	os.system(blastcommand)
	
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	simmatrix={}
	for seqid in seqrecords.keys():
		simmatrix.setdefault(seqid,{})
		simmatrix[seqid][seqid]=1
	
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

if __name__ == "__main__":
	#load sequences
	base = sys.argv[0][0: sys.argv[0].rindex("/")+1]
	simmatrix={}
	print("Computing similarity matrix using BLAST..")
	simmatrix=ComputeSim(fastafilename,mincoverage,minsim)
	output=GetWorkingBase(fastafilename) + ".sim"
	#save the simmatrix
	SaveSim(simmatrix,output,minsim)
	print("The similarity file is saved in " + output + ".")		
			
		
