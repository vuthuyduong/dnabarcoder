#!/usr/bin/env python
# FILE: classify.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
#from sklearn.model_selection import StratifiedKFold
#from sklearn.preprocessing import StandardScaler
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse

from dnabarcoder.classification.search import ComputeBestBLASTscore
from dnabarcoder.classification.classify import GetAssignment
from dnabarcoder.classification.classify import LoadClassification
from dnabarcoder.classification.classify import LoadClassificationFromDescription
from dnabarcoder.classification.classify import AddCutoffsToTaxonomy

import json
from Bio import SeqIO
import multiprocessing

def main():
	parser=argparse.ArgumentParser(prog='classifyOne.py',  
								   usage="%(prog)s [options] -i sequences -r referencepath -o output",
								   description='''Script that assigns the classified sequences of the prediction file to their BLAST best match based on the given cutoffs.''',
								   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
	   )

	parser.add_argument('-i','--input', help='The input can be a fasta file name or a text of sequences for classification')
	parser.add_argument('-r','--referencepath', default="references", help='The path to the folder containing reference datasets and their associated similarity cutoffs files.')
	#parser.add_argument('-f','--fasta', default="", help='The fasta file of the sequences for saving unidentified sequences. Optional.')
	#parser.add_argument('-c','--classification', default="", help='the classification file in tab. format.')
	#parser.add_argument('-r','--reference', default="", help='the reference fasta file, in case the classification of the sequences is given in the sequence headers.')
	parser.add_argument('-o','--out', default="dnabarcoder_output", help='The output folder.')
	#parser.add_argument('-fmt','--inputformat', default="tab delimited", help='the format of the classified file. The inputfmt can have two values "tab delimited" and "blast". The value "tab delimited" is given as default, and the "blast" fmt is the format of the BLAST output with outfmt=6.')
	#parser.add_argument('-cutoff','--globalcutoff', type=float, default=-1,help='The global cutoff to assign the sequences to predicted taxa. If the cutoffs file is not given, this value will be taken for sequence assignment.')
	#parser.add_argument('-confidence','--globalconfidence', type=float,default=-1,help='The global confidence to assign the sequences to predicted taxa')
	#parser.add_argument('-rank','--classificationrank', default="", help='the classification rank')
	#parser.add_argument('-prefix','--prefix', help='the prefix of output filenames')
	#parser.add_argument('-cutoffs','--cutoffs', help='The json file containing the local cutoffs to assign the sequences to the predicted taxa.')
	#parser.add_argument('-minseqno','--minseqno', type=int, default=0, help='the minimum number of sequences for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
	#parser.add_argument('-mingroupno','--mingroupno', type=int, default=0, help='the minimum number of groups for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
	#parser.add_argument('-ml','--minalignmentlength', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should probably be set to smaller, 50 for instance.')
	#parser.add_argument('-saveclassifiedonly','--saveclassifiedonly',default=False, help='The option to save all (False) or only classified sequences (True) in the classification output.')
	parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')
	#parser.add_argument('-display','--display',default="", help='If display=="yes" then the krona html is displayed.')
	parser.add_argument('-ncpus','--ncpus', type=int, default=0, help='The number of CPUs used for searching. The default value is the total number of CPUs.')

	args=parser.parse_args()
	sequences=args.input
	referencepath=args.referencepath
	#globalcutoff=args.globalcutoff
	#globalconfidence=args.globalconfidence
	#cutoffsfilename=args.cutoffs
	#classificationfilename=args.classification
	#classificationrank=args.classificationrank
	#fastafilename= args.fasta
	#referencefastafilename= args.reference
	#mincoverage = args.minalignmentlength
	#prefix=args.prefix
	nproc=args.ncpus
	if nproc==0:
		nproc=multiprocessing.cpu_count()
	outputpath=args.out
	
	finalbestmatchdict=classifyAgainstReferenceSet(sequences,referencepath,args.idcolumnname,outputpath,nproc)
	classification_result=Dict2Tab(finalbestmatchdict)
	outputfilename=outputpath + "/query.classification"
	if os.path.exists(sequences):
		outputfilename=outputpath + "/" + os.path.basename(sequences)[0:os.path.basename(sequences).index(".")] + ".classification"
	with open(outputfilename, "w") as f:
		f.write(classification_result)
	print("The classification output is saved in " + outputfilename + ".")	
	
	#classifyOneSequence(sequence,referencefastafilename,classificationfilename,args.idcolumnname,cutoffsfilename,globalcutoff,globalconfidence,mincoverage,args.minseqno,args.mingroupno,classificationrank,outputpath,nproc)
	
def Assign(refclassificationdict,taxonomy,bestmatchdict,classificationrank):
	#classificationlevel=GetLevel(classificationrank)
	#classification_result="ID\tReferenceID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\trank\tscore\tcutoff\tconfidence\n"
	classificationdict={}
	
	for seqid in bestmatchdict.keys():
		classificationdict.setdefault(seqid,{})
		rank=""
		level=-1
		refid=bestmatchdict[seqid]["refid"]
		bestscore=bestmatchdict[seqid]["score"]
		#sim=bestmatchdict[seqid]["sim"]
		#coverage=bestmatchdict[seqid]["alignmentlength"]
		confidence=-1
		cutoff=-1
		classification="k__unidentified;p__unidentified;c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified"
		if refid!="":
			classification,predictedname,rank,level,cutoff,confidence=GetAssignment(refid,refclassificationdict,bestscore,taxonomy,classificationrank)		
		bestmatchdict[seqid]["cutoff"]=cutoff
		bestmatchdict[seqid]["confidence"]=confidence
		if sys.version_info[0] < 3:
			predictedname=unicode(predictedname,'latin1')
		cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","").replace("_"," ")
		cleanclassification=cleanclassification.replace(";","\t")
		#classification_result=classification_result + seqid + "\t" + refid + "\t" + cleanclassification + "\t" + rank + "\t" + str(bestscore) + "\t" + cutoff_str + "\t" + confidence_str + "\n"
		taxa=cleanclassification.split("\t")
		bestmatchdict[seqid]["kingdom"]=taxa[0]
		bestmatchdict[seqid]["phylum"]=taxa[1]
		bestmatchdict[seqid]["class"]=taxa[2]
		bestmatchdict[seqid]["order"]=taxa[3]
		bestmatchdict[seqid]["family"]=taxa[4]
		bestmatchdict[seqid]["genus"]=taxa[5]
		bestmatchdict[seqid]["species"]=taxa[6]

def classifySequences(sequences,referencefastafilename,classificationfilename,idcolumnname,cutoffsfilename,globalcutoff,globalconfidence,mincoverage,minseqno,mingroupno,classificationrank,outputpath,nproc):
	if not os.path.exists(outputpath):
		os.system("mkdir " + outputpath)
	
	query = outputpath + "/query.fasta"
	if os.path.exists(sequences): #input is a fasta file
		os.system("cp " + sequences + " " + query)
	else: #the input is a text of sequences	
		if not (">" in sequences):
			command = "echo \">query\n" + sequences + "\" > " + query
		else:
			command = "echo \"" + sequences + "\" > " + query
		os.system(command)
	queryrecords=SeqIO.to_dict(SeqIO.parse(query, "fasta"))	
	#search for a best match of a test sequence in a train dataset
	bestmatchlist,bestscorelist,bestsimlist,bestcoveragelist=ComputeBestBLASTscore(query,referencefastafilename,mincoverage,outputpath,nproc)
	bestmatchdict={}
	for seqid in queryrecords.keys():
		bestmatchdict.setdefault(seqid,{})
		bestmatchdict[seqid]["refid"]=bestmatchlist[0]
		bestmatchdict[seqid]["score"]=bestscorelist[0]
		bestmatchdict[seqid]["sim"]=bestsimlist[0]
		bestmatchdict[seqid]["coverage"]=bestcoveragelist[0]
	
	refclassificationdict={}
	#load classification for the sequences
	if classificationfilename!="":
		refclassificationdict,taxonomy,isError = LoadClassification(classificationfilename,idcolumnname)
		if isError==True:
			sys.exit()
	else:
		refclassificationdict,taxonomy = LoadClassificationFromDescription(classificationfilename)
	
	cutoffs={}
	if cutoffsfilename!="" and cutoffsfilename!=None:
		with open(cutoffsfilename) as cutoffsfile:
			cutoffs = json.load(cutoffsfile)	
	#add cutoffs to taxa for sequence identification		
	AddCutoffsToTaxonomy(taxonomy,globalcutoff,globalconfidence,cutoffs,minseqno,mingroupno)
	Assign(refclassificationdict,taxonomy,bestmatchdict,classificationrank)

	return bestmatchdict
	
def getReferences(referencepath):
	referencedict={}
	filenames = os.listdir(referencepath)

	for filename in filenames:
		reference=filename.split(".")[0]
		if ".json" in filename:
			referencedict.setdefault(reference,{})
	for filename in filenames:
		reference=filename.split(".")[0]
		if not reference in referencedict.keys():
			continue
		if ".json" in filename:
			referencedict[reference]["cutoffs"]=referencepath + "/" + filename
		if ".classification" in filename:
			referencedict[reference]["classificationfilename"]=referencepath + "/" + filename	
		if ".fas" in filename:
			referencedict[reference]["referencefilename"]=referencepath + "/" + filename
		if ("ITS1" in filename) or ("ITS2" in filename):
			referencedict[reference]["mincoverage"]=50
		else:
			referencedict[reference]["mincoverage"]=400
		if ("speciescomplexes" in filename):
			referencedict[reference]["speciescomplexes"]=referencepath + "/" + filename	
		referencedict[reference]["minseqno"]=0
		referencedict[reference]["mingroupno"]=0
	return referencedict

def Dict2Tab(finalbestmatchdict):
	#output 				
	classification_result="ID\tReferenceID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\trank\tscore\tsim\tcoverage\tcutoff\tconfidence\n"	
	for seqid in finalbestmatchdict.keys():
		classification_result=classification_result + seqid
		classification_result=classification_result + "\t" + finalbestmatchdict[seqid]["refid"]
		classification_result=classification_result + "\t" + finalbestmatchdict[seqid]["kingdom"]
		classification_result=classification_result + "\t" + finalbestmatchdict[seqid]["phylum"]
		classification_result=classification_result + "\t" + finalbestmatchdict[seqid]["order"]
		classification_result=classification_result + "\t" + finalbestmatchdict[seqid]["family"]
		classification_result=classification_result + "\t" + finalbestmatchdict[seqid]["genus"]
		classification_result=classification_result + "\t" + finalbestmatchdict[seqid]["species"]
		classification_result=classification_result + "\t" + str(finalbestmatchdict[seqid]["score"]) + "\t" + str(finalbestmatchdict[seqid]["sim"]) + "\t" +  str(finalbestmatchdict[seqid]["coverage"]) + "\t"+ str(finalbestmatchdict[seqid]["cutoff"]) + "\t" + str(finalbestmatchdict[seqid]["confidence"]) + "\n"
	return classification_result

def classifyAgainstReferenceSet(sequences,referencepath,idcolumnname,outputpath,nproc):
	referencedict=getReferences(referencepath)
	finalbestmatchdict={}
	for reference in referencedict.keys():
		referencefilename=referencedict[reference]["referencefilename"]
		classificationfilename=referencedict[reference]["classificationfilename"]
		cutoffsfilename=referencedict[reference]["cutoffs"]
		mincoverage=referencedict[reference]["mincoverage"]
		minseqno=referencedict[reference]["minseqno"]
		mingroupno=referencedict[reference]["mingroupno"]
		globalcutoff=0
		globalconfidence=0
		classificationrank=""
		bestmatchdict=classifySequences(sequences,referencefilename,classificationfilename,idcolumnname,cutoffsfilename,globalcutoff,globalconfidence,mincoverage,minseqno,mingroupno,classificationrank,outputpath,nproc)
		if finalbestmatchdict=={}:
			finalbestmatchdict=bestmatchdict.copy()
		else:
			for seqid in bestmatchdict.keys():
				if seqid in finalbestmatchdict.keys():
					if finalbestmatchdict[seqid]["score"] < bestmatchdict[seqid]["score"]:
						finalbestmatchdict[seqid]=bestmatchdict[seqid].copy()
					elif (finalbestmatchdict[seqid]["score"] == bestmatchdict[seqid]["score"]) and (finalbestmatchdict[seqid]["coverage"] < bestmatchdict[seqid]["coverage"]):
						finalbestmatchdict[seqid]=bestmatchdict[seqid].copy()
				else:
					finalbestmatchdict.setdefault(seqid,bestmatchdict[seqid])
		#add indistinguishable species to finalbestmatchdict
		for seqid in finalbestmatchdict.keys():
			species=finalbestmatchdict[seqid]["species"]
			if species=="" or species=="unidentified":
				continue
			if "speciescomplexes" in referencedict[reference].keys():
				speciescomplexes={}
				with open(referencedict[reference]["speciescomplexes"], 'r') as f:
					speciescomplexes=json.load(f)
				try:
					finalbestmatchdict[seqid]["indistinguishable species"]=speciescomplexes[species]
					finalbestmatchdict[seqid]["indistinguishable species"]
				except KeyError:
					pass
	
	#output 				
	#classification_result=Dict2Tab(finalbestmatchdict)	
	return finalbestmatchdict
if __name__ == "__main__":
	main()
	
	