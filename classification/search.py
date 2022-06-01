#!/usr/bin/env python
# FILE: search.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
import sys
#from sklearn.metrics import precision_recall_fscore_support
#from sklearn.metrics import cohen_kappa_score
#from sklearn.metrics import matthews_corrcoef
#from sklearn.metrics import confusion_matrix
#from sklearn.metrics import accuracy_score
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
from Bio import SeqIO
#import json
import multiprocessing

nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='search.py',  
							   usage="%(prog)s [options] -i fastafile -r referencefastafile -ml minalignmentlength ",
							   description='''Script that classifies the sequences of the fasta files using BLAST with a cut-off value or the cut-off values given in the cutoffs file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be classified.')
parser.add_argument('-r','--reference', required=True, help='the reference fasta file.')
parser.add_argument('-ml','--minalignmentlength', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should be set to smaller, 50 for instance.')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-prefix','--prefix', help='the prefix of output filenames')

args=parser.parse_args()
testdataset= args.input
traindataset = args.reference
mincoverage = args.minalignmentlength

prefix=args.prefix
outputpath=args.out
if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]
	
def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	if "." in basename:
		basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + filename
	return path

def GetSeqIndex(seqname,seqrecords):
	i=0
	for seqrecord in seqrecords:
		if (seqname == seqrecord.id):
			return i
		i = i + 1
	return -1

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

def ComputeBestBLASTscore(query,reference,mincoverage):
	indexed_query= IndexSequences(query)

	#load sequeces from the fasta files
	queryrecords = list(SeqIO.parse(indexed_query, "fasta"))
	#refrecords = list(SeqIO.parse(reference, "fasta"))

	bestscorelist =[0] * len(queryrecords)
	bestsimlist =[0] * len(queryrecords)
	bestcoveragelist =[0] * len(queryrecords)
	bestrefidlist = [""] * len(queryrecords)

	#blast
	makedbcommand = "makeblastdb -in " + reference + " -dbtype \'nucl\' " +  " -out db"
	print(makedbcommand)
	os.system(makedbcommand)
	#for short read
	blastcommand = "blastn -query " + indexed_query + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	#for long read
	if mincoverage >=400:
		blastcommand = "blastn -query " + indexed_query + " -db  db -outfmt 6 -out out.txt -num_threads " + str(nproc)
	print(blastcommand)	
	os.system(blastcommand)
	
	#read blast output
	blastoutputfile = open("out.txt")
	refid = ""
	score=0
	queryid=""
	for line in blastoutputfile:
		words = line.split("\t")
		queryid=words[0]
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		refid=words[1]
		score=sim
		if coverage < mincoverage:
				score=float(score * coverage)/mincoverage
		i = int(queryid.split("|")[0])		
		#if score > bestscorelist[i]:
		if score > bestscorelist[i] or (score == bestscorelist[i] and coverage > bestcoveragelist[i]):	
			bestscorelist[i]= score
			bestrefidlist[i]=refid
			bestsimlist[i]=sim
			bestcoveragelist[i]=coverage
	os.system("rm " + indexed_query)		
	os.system("rm out.txt")
	return bestrefidlist,bestscorelist,bestsimlist,bestcoveragelist

def SavePrediction(testseqIDs,bestscorelist,bestsimlist,bestcoveragelist,bestrefidlist,outputname):
	output=open(outputname,"w")
	output.write("ID\tReferenceID\tBLAST score\tBLAST sim\tBLAST coverage\n")
	i=0
	for seqid in testseqIDs:
		output.write(seqid + "\t"  + bestrefidlist[i] + "\t" +  str(bestscorelist[i]) + "\t" + str(bestsimlist[i]) + "\t" + str(bestcoveragelist[i]) +"\n")
		i=i+1
	output.close()
	
##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]

#load ref seq records
refseqrecords = SeqIO.to_dict(SeqIO.parse(traindataset, "fasta"))

#load test seq records
testseqrecords = SeqIO.to_dict(SeqIO.parse(testdataset, "fasta"))

#search for a best match of a test sequence in a train dataset
bestmatchlist,bestscorelist,bestsimlist,bestcoveragelist=ComputeBestBLASTscore(testdataset,traindataset,mincoverage)

#Save prediction by searching 
if prefix=="" or prefix==None:
	prefix=GetBase(testdataset)
	if "/" in prefix:
		prefix=prefix[prefix.rindex("/")+1:]	
basename=GetBase(traindataset)
if "/" in basename:
	basename=basename[basename.rindex("/")+1:]		
reportfilename=GetWorkingBase(prefix) + "." + basename + "_BLAST.bestmatch"
SavePrediction(testseqrecords.keys(),bestscorelist,bestsimlist,bestcoveragelist,bestmatchlist,reportfilename)
print("The results are saved in file  " + reportfilename)

