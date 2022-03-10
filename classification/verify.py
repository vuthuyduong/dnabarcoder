#!/usr/bin/env python
# FILE: verify.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
#from sklearn.model_selection import StratifiedKFold
#from sklearn.preprocessing import StandardScaler
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
#from sklearn.metrics import precision_recall_fscore_support
#from sklearn.metrics import cohen_kappa_score
#from sklearn.metrics import matthews_corrcoef
#from sklearn.metrics import confusion_matrix
#from sklearn.metrics import accuracy_score
#import json
from Bio import SeqIO
from Bio import Phylo
import pylab
import random
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='verify.py',  
							   usage="%(prog)s [options] -i classified file -f the fasta file -r referencefastafilename -c classificationfile -o output",
							   description='''Script that assigns the classified sequences of the prediction file to their BLAST best match based on the given cutoffs.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the classified file')
parser.add_argument('-f','--fasta', required=True, help='the fasta file')
parser.add_argument('-r','--reference', required=True, help='the reference fasta file')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', default="", help='the classification file in tab. format.')
parser.add_argument('-m','--maxseqno', type=int, default=50, help='Maximum number of the sequences of the predicted taxon name from the classification file will be selected for the comparison to find the best match. If it is not given, all the sequences will be selected.')
parser.add_argument('-rank','--classificationrank', default="", help='the classification rank')
parser.add_argument('-redo','--redo', default="", help='the classification rank')
parser.add_argument('-prefix','--prefix', help='the prefix of output filenames')

args=parser.parse_args()
predictionfilename=args.input
fastafilename= args.fasta
referencefastafilename= args.reference
classificationfilename=args.classification
maxseqno=args.maxseqno
prefix=args.prefix
verifyingrank=args.classificationrank
redo=args.redo
outputpath=args.out

if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)
nproc=multiprocessing.cpu_count()

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]
	
def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	if "." in basename:
		basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + filename
	return path

def GetTaxonomicClassification(level,header,texts):
	classification=""
	p_s=len(texts)
	p_g=len(texts)
	p_f=len(texts)
	p_o=len(texts)
	p_c=len(texts)
	p_p=len(texts)
	p_k=len(texts)
	i=0
	for text in header.split("\t"):
		text=text.rstrip()
		if text.lower()=="species":
			p_s=i
		elif text.lower()=="genus":
			p_g=i	
		elif text.lower()=="family":
			p_f=i	
		elif text.lower()=="order":
			p_o=i	
		elif text.lower()=="class":
			p_c=i	
		elif text.lower()=="phylum":
			p_p=i	
		elif text.lower()=="kingdom":
			p_k=i	
		i=i+1 
	species="unidentified"
	genus="unidentified"
	family="unidentified"
	order="unidentified" 
	bioclass="unidentified"
	phylum="unidentified"
	kingdom="unidentified"
	if p_s< len(texts):
		species=texts[p_s].rstrip()
		if species=="" or ("unculture" in species) or ("_sp_" in species) or ("unidentified" in species) or ("environmental" in species.lower()) or ("(" in species):
			species="unidentified"
		elif not (" " in species or "_" in species):
			species="unidentified"
		species=species.replace(" ","_")
	if p_g< len(texts):
		genus=texts[p_g].rstrip()
		if genus=="" or ("unculture" in genus) or ("unidentified" in genus):
			genus="unidentified"
	if p_f< len(texts):
		family=texts[p_f].rstrip()
		if family=="" or ("unculture" in family) or ("unidentified" in family):
			family="unidentified"
	if p_o< len(texts):
		order=texts[p_o].rstrip()
		if order=="" or ("unculture" in order) or ("unidentified" in order):
			order="unidentified"
	if p_c< len(texts):
		bioclass=texts[p_c].rstrip()
		if bioclass=="" or ("unculture" in bioclass) or ("unidentified" in bioclass):
			bioclass="unidentified"
	if p_p< len(texts):
		phylum=texts[p_p].rstrip()
		if phylum=="" or ("unculture" in phylum) or ("unidentified" in phylum):
			phylum="unidentified"
	if p_k< len(texts):
		kingdom=texts[p_k].rstrip()
		if kingdom=="" or ("unculture" in kingdom) or ("unidentified" in kingdom):
			kingdom="unidentified"
	taxonname=""
	rank=""
	if level <7 and kingdom!="unidentified":
		taxonname=kingdom
		rank="kingdom"
		classification="k__" + kingdom +";p__unidentified;c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified"	
	if level <6 and phylum!="unidentified":
		taxonname=phylum
		rank="phylum"
		classification="k__" + kingdom +";p__"+phylum +";c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified"
	if level <5 and bioclass!="unidentified":
		taxonname=bioclass
		rank="class"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass+";o__unidentified;f__unidentified;g__unidentified;s__unidentified"
	if level <4 and order!="unidentified":
		taxonname=order
		rank="order"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order + ";f__unidentified;g__unidentified;s__unidentified"
	if level <3 and family!="unidentified":
		taxonname=family
		rank="family"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family +";g__unidentified;s__unidentified"
	if level <2 and genus!="unidentified":
		taxonname=genus
		rank="genus"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family + ";g__"+ genus + ";s__unidentified"
	if level <1 and species!="unidentified":
		taxonname=species
		rank="species"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family + ";g__"+ genus+";s__"+species 	
	return classification,taxonname,rank

def LoadClassification(seqrecords,classificationfilename):
	classificationdict={}
	classes={}
	if classificationfilename == "":
		return {},{}
	classificationfile= open(classificationfilename)
	header=next(classificationfile)
	for line in classificationfile:
		texts=line.split("\t")
		seqid = texts[0].replace(">","").rstrip()
		classification,taxonname,rank=GetTaxonomicClassification(0,header,texts)
		if classification!="":
			classificationdict.setdefault(seqid,{})
			classificationdict[seqid]["classification"]=classification
			classificationdict[seqid]["taxonname"]=taxonname
			classificationdict[seqid]["rank"]=rank	
		if seqid in seqrecords.keys():
			taxonnames=classification.split(";")
			for taxonname in taxonnames:
				if taxonname=="":
					continue
				taxonname=taxonname.split("__")[1].replace("_"," ")
				if taxonname=="unidentified" or taxonname=="":
					continue
				if not taxonname in classes.keys():
					classes.setdefault(taxonname,{})
				classes[taxonname][seqid]=seqrecords[seqid]
	classificationfile.close()	
	return classificationdict,classes

def LoadPrediction(predictionfilename):
	predictiondict={}
	p_id=-1
	p_l=-1
	p_pl=-1
	p_c=-1
	p_p=-1
	p_r=-1
	p_bs=-1
	p_s=-1
	p_co=-1
	p_rank=-1
	p_v=-1
	p_t=-1
	p_cutoff=-1
	p_confidence=-1
	p_numberofrefsequences=-1
	p_branchlength=-1
	p_maxbranchlength=-1
	p_averagebranchlength=-1
	predictionfile=open(predictionfilename)
	headers=next(predictionfile)
	i=0
	for header in headers.rstrip().split("\t"):
		if ("sequenceid" in header.lower()) or ("sequence_id" in header.lower()) or ("sequence id" in header.lower()) :
			p_id=i
		if "given label" in header.lower():
			p_l=i
		if "prediction" in header.lower():
			p_pl=i
		if "full classification" in header.lower():
			p_c=i
		if "probability" in header.lower():
			p_p=i
		if ("referenceid" in header.lower()) or ("reference_id" in header.lower()) or ("reference id" in header.lower()) :
			p_r=i
		if "score" in header.lower():
			p_bs=i	
		if "sim" in header.lower():
			p_s=i
		if "coverage" in header.lower():
			p_co=i	
		if "rank" in header.lower():
			p_rank=i
		if "cut-off" in header.lower():
			p_cutoff=i	
		if "confidence" in header.lower():
			p_confidence=i		
		if "verified label" in header.lower():
			p_v=i
		if "tree filename" in header.lower():
			p_t=i	
		if "number of reference sequences" in header.lower():
			p_numberofrefsequences=i
		if "branch length" in header.lower():
			p_branchlength=i
		if "max branch length" in header.lower():
			p_maxbranchlength=i	
		if "average branch length" in header.lower():
			p_averagebranchlength=i	
		i=i+1	
	for line in predictionfile:
		texts=line.rstrip().split("\t")
		if p_id >=0 and p_id < len(texts):
			seqid=texts[p_id]
			predictiondict.setdefault(seqid,{})
		else:
			continue
		label=""
		if p_l >=0 and p_l < len(texts):
			label=texts[p_l]
		predictiondict[seqid]["givenlabel"]=label
		predlabel=""
		if p_pl >=0 and p_pl < len(texts):
			predlabel=texts[p_pl]
		predictiondict[seqid]["predlabel"]=predlabel
		pred_classification=""
		if p_c >=0 and p_c < len(texts):
			pred_classification=texts[p_c]
		predictiondict[seqid]["classification"]=pred_classification
		proba=1
		if p_p >=0 and p_p <len(texts):
			proba=float(texts[p_p])
		predictiondict[seqid]["proba"]=proba
		refid=""
		if p_r>0 and p_r<len(texts):
			refid=texts[p_r]
		predictiondict[seqid]["refid"]=refid
		score=0
		if p_bs>0 and p_bs<len(texts):
			score=float(texts[p_bs])
		predictiondict[seqid]["score"]=score
		sim=0
		if p_s>0 and p_s<len(texts):
			sim=float(texts[p_s])
		predictiondict[seqid]["sim"]=sim
		cov=""
		if p_co>0 and p_co<len(texts):
			cov=float(texts[p_co])
		predictiondict[seqid]["coverage"]=cov
		rank=""
		if p_rank>0 and p_rank<len(texts):
			rank=texts[p_rank]
		predictiondict[seqid]["rank"]=rank
		cutoff=0
		if p_cutoff>0 and p_cutoff<len(texts):
			cutoff=texts[p_cutoff]
		predictiondict[seqid]["cutoff"]=cutoff
		confidence=0
		if p_confidence>0 and p_confidence<len(texts):
			confidence=texts[p_confidence]
		predictiondict[seqid]["confidence"]=confidence
		verifiedlabel=""
		if p_v>0 and p_v<len(texts):
			verifiedlabel=texts[p_v]
		predictiondict[seqid]["verifiedlabel"]=verifiedlabel
		treefilename=""
		if p_t>0 and p_t<len(texts):
			treefilename=texts[p_t]
		predictiondict[seqid]["treefilename"]=treefilename
		numberofrefsequences=0
		if p_numberofrefsequences>0 and p_numberofrefsequences<len(texts):
			numberofrefsequences=texts[p_numberofrefsequences]
		predictiondict[seqid]["numberofrefsequences"]=numberofrefsequences
		branchlength=0
		if p_branchlength>0 and p_branchlength<len(texts):
			branchlength=texts[p_branchlength]
		predictiondict[seqid]["branchlength"]=branchlength
		maxbranchlength=0
		if p_maxbranchlength>0 and p_maxbranchlength<len(texts):
			maxbranchlength=texts[p_maxbranchlength]
		predictiondict[seqid]["maxbranchlength"]=maxbranchlength
		averagebranchlength=0
		if p_averagebranchlength>0 and p_averagebranchlength<len(texts):
			averagebranchlength=texts[p_averagebranchlength]
		predictiondict[seqid]["averagebranchlength"]=averagebranchlength
		
	return predictiondict

def GetLevel(rank):
	level=-1
	if rank=="species":	
		level=6
	elif rank=="genus":	
		level=5
	elif rank=="family":	
		level=4
	elif rank=="order":	
		level=3
	elif rank=="class":	
		level=2
	elif rank=="phylum":
		level=1
	elif rank=="kingdom":
		level=0	
	return level

def lookup_by_names(tree):
	names = {}
	for clade in tree.find_clades():
		if clade.name:
			if clade.name in names:
				raise ValueError("Duplicate key: %s" % clade.name)
			names[clade.name] = clade
	return names

def verifyBasedOnBranchLengths(seqid,treefilename):
	tree = Phylo.read(treefilename, "newick")
	names = lookup_by_names(tree)
	max_length=0
	length=0
	average=0
	n=0
	for name in names:
		clade=names[name]
		if name==seqid:
			length=clade.branch_length
		else:
			average=average + clade.branch_length
			n=n+1
			if clade.branch_length > max_length:
				max_length=clade.branch_length
	verified=(length<=max_length) and (len(names) >=3)
	if n>0:
		average=round(average/n,2)
	return verified,length,max_length,average

def PrintTree(treefilename,redo):
	svgfilename=GetWorkingBase(fastafilename) + ".tree.svg"
	if (not os.path.exists(svgfilename)) or redo!="":
		tree = Phylo.read(treefilename, "newick")
		Phylo.draw(tree,do_show=False)	
		#Phylo.draw_ascii(tree,do_show=False)	
		pylab.savefig(svgfilename,format='svg', bbox_inches='tight', dpi=300)
		print("A iq-tree and its svg file are saved in " + treefilename + " and " + svgfilename + ".")	
	
def CreateTree(fastafilename,redo):
	alignmentfilename =  GetBase(fastafilename) + ".aligned.fas"
	if (not os.path.exists(alignmentfilename)) or redo!="":
		#make the alignment
		if not os.path.exists(alignmentfilename):
			command="clustalo -i " + fastafilename + " -o " + alignmentfilename
			os.system(command)
	treefilename=GetBase(fastafilename) + ".aligned.fas.treefile"
	if (not os.path.exists(treefilename)) or redo!="":
		#make tree
		command="iqtree -pers 0.2 -n 500 -s " + alignmentfilename
		#command="iqtree -s " + alignmentfilename
		os.system(command)
	#print tree
	PrintTree(treefilename,redo)
	return treefilename	

def CreateFastaFile(seqrecord,taxonname,classeswithsequences,numberofrefsequences,maxseqno,redo):
	if not os.path.exists(outputpath + "/trees"):
		os.system("mkdir " + outputpath + "/trees")
	if sys.version_info[0] < 3:
		taxonname=unicode(taxonname,errors='ignore')
	newfastafilename=outputpath + "/trees/" + seqrecord.id.replace("|","_") + ".fasta"
	if not os.path.exists(newfastafilename) or redo!="":
		if taxonname in classeswithsequences.keys():
			seqrecords=[]
			seqrecords.append(seqrecord)
			sequences=classeswithsequences[taxonname]
			numberofrefsequences=len(sequences.keys())
			if len(sequences) >=2:#only make tree of more than 3 sequences
				if (maxseqno >0) and (len(sequences) > maxseqno):
					#select randomly 100 sequences to compare
					selectedlist=random.sample(range(0, len(sequences)), k=maxseqno)
					for i in selectedlist:
						seqids=list(sequences.keys())
						sequenceid=seqids[i]
						if sequenceid != seqrecord.id:
							seqrecords.append(sequences[sequenceid])
				else:	
					for sequenceid in sequences.keys():
						if sequenceid != seqrecord.id:
							seqrecords.append(sequences[sequenceid])
				if len(seqrecords) >=3:			
					SeqIO.write(seqrecords,newfastafilename,"fasta")
				else:
					newfastafilename=""
			else:
				newfastafilename=""
	if numberofrefsequences==0 and os.path.exists(newfastafilename):
		records=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
		numberofrefsequences=len(records)-1
	return newfastafilename,numberofrefsequences
	
def Verify(seqrecords,predictiondict,refclasses,maxseqno,verifyingrank):
	count=0
	notree_count=0
	total=0
	#create Fasta files
	for seqid in predictiondict.keys():
		predictedname=predictiondict[seqid]["predlabel"]
		rank=predictiondict[seqid]["rank"]
		treefilename=predictiondict[seqid]["treefilename"]
		numberofrefsequences=predictiondict[seqid]["numberofrefsequences"]
		branchlength=predictiondict[seqid]["branchlength"]
		maxbranchlength=predictiondict[seqid]["maxbranchlength"]
		averagebranchlength=predictiondict[seqid]["averagebranchlength"]
		#if sys.version_info[0] < 3:
			#predictedname=unicode(predictedname,'latin1')
		verifiedlabel=""
		if (verifyingrank=="" or (verifyingrank!="" and rank==verifyingrank)):
			#only predict when the tree file name does not exist
			if predictedname in refclasses.keys() and seqid in seqrecords.keys():
				total=total+1
				seqrecord=seqrecords[seqid]
				if treefilename=="" or (redo !=""):
					treefastafilename,numberofrefsequences=CreateFastaFile(seqrecord,predictedname,refclasses,numberofrefsequences,maxseqno,redo)
					if os.path.exists(treefastafilename):
						treefilename= CreateTree(treefastafilename,redo)
				verified=False		
				if 	os.path.exists(treefilename):	
					verified,branchlength,maxbranchlength,averagebrachlength=verifyBasedOnBranchLengths(seqid,treefilename)
				else:
					notree_count=notree_count+1
				if verified==True:
					verifiedlabel=predictedname
					count=count+1
		predictiondict[seqid]["verifiedlabel"]=verifiedlabel
		predictiondict[seqid]["treefilename"]=treefilename		
		predictiondict[seqid]["numberofrefsequences"]=numberofrefsequences
		predictiondict[seqid]["branchlength"]=branchlength
		predictiondict[seqid]["maxbranchlength"]=maxbranchlength
		predictiondict[seqid]["averagebranchlength"]=averagebranchlength
	return count,notree_count,total

def SaveVerification(predictiondict,output,notverifiedoutput,classificationfilename):
	classificationreportfile=open(classificationreportfilename,"w")
	outputfile=open(output,"w")
	notverifiedoutputfile=open(notverifiedoutput,"w")
	notverifiedoutputfile.write("#SequenceID\tGiven label\tPrediction\tFull classification\tProbability\tRank\tCut-off\tConfidence\tReferenceID\tBLAST score\tBLAST sim\tBLAST coverage\tVerified label\tTree filename\tNumber of reference sequences\tBranch length\tMax branch length\tAverage branch length\n")
	outputfile.write("#SequenceID\tGiven label\tPrediction\tFull classification\tProbability\tRank\tCut-off\tConfidence\tReferenceID\tBLAST score\tBLAST sim\tBLAST coverage\tVerified label\tTree filename\tNumber of reference sequences\tBranch length\tMax branch length\tAverage branch length\n")
	classificationreportfile.write("SequenceID\tReferenceID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\trank\tscore\tcutoff\tconfidence\n")
	for seqid in predictiondict.keys():
		prediction=predictiondict[seqid]
		treefilename=prediction["treefilename"]
		predlabel=prediction["predlabel"]
		givenlabel=prediction["givenlabel"]
		classification=prediction["classification"]
		rank=prediction["rank"]
		bestscore=prediction["score"]
		cutoff=prediction["cutoff"]
		confidence=prediction["confidence"]
		proba=prediction["proba"]
		refid=prediction["refid"]
		sim=prediction["sim"]
		coverage=prediction["coverage"]
		verifiedlabel=prediction["verifiedlabel"]
		numberofrefsequences=prediction["numberofrefsequences"]
		branchlength=prediction["branchlength"]
		maxbranchlength=prediction["maxbranchlength"]
		averagebranchlength=prediction["averagebranchlength"]
		if verifiedlabel!="" and verifiedlabel!="unidentified":			
			outputfile.write(seqid + "\t" + givenlabel + "\t"  + predlabel + "\t"+ classification + "\t" + str(proba) + "\t" + rank + "\t" + str(cutoff) + "\t" + str(confidence) + "\t" + refid + "\t" + str(bestscore) + "\t" + str(sim) + "\t" + str(coverage) + "\t" + verifiedlabel + "\t" + treefilename + "\t" + str(numberofrefsequences) + "\t" + str(branchlength) + "\t" + str(maxbranchlength) + "\t" + str(averagebranchlength) + "\n")			
			cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","").replace("_"," ")
			if refid==seqid:
				refid=""
			classificationreportfile.write(seqid + "\t" + refid + "\t" + cleanclassification.replace(";","\t") + "\t" + rank + "\t" + str(bestscore) + "\t" + str(cutoff) + "\t" + str(confidence) + "\n")
		elif verifyingrank!="" and verifyingrank==rank:
			notverifiedoutputfile.write(seqid + "\t" + givenlabel + "\t"  + predlabel + "\t"+ classification + "\t" + str(proba) + "\t" + rank + "\t" + str(cutoff) + "\t" + str(confidence) + "\t" + refid + "\t" + str(bestscore) + "\t" + str(sim) + "\t" + str(coverage) + "\t" + verifiedlabel + "\t" + treefilename + "\t" + str(numberofrefsequences) + "\t" + str(branchlength) + "\t" + str(maxbranchlength) + "\t" + str(averagebranchlength) + "\n")			
	classificationreportfile.close()
	outputfile.close()	
	notverifiedoutputfile.close()	
		

def LoadClassificationForKronaReport(classificationfilename):
	classificationdict={}
	classificationfile= open(classificationfilename)
	next(classificationfile)
	for line in classificationfile:
		if line.startswith("#"):
			continue 
		texts=line.split("\t")
		classification=texts[3]
		classification=classification.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","")
		if ";" in classification:
			texts=classification.split(";")
			classification=texts[0]
			if classification=="":
				classification="unidentified"
			for i in range(1,len(texts)):
				text=texts[i]
				if text=="":
					text="unidentified"
				classification=classification + "\t" + text	
		if classification in classificationdict.keys():
			classificationdict[classification]=classificationdict[classification] + 1
		else:
			classificationdict.setdefault(classification, 1)
	classificationfile.close()	
	return classificationdict

def KronaPieCharts(classification,kronareport,kronahtml):
	kronareportfile=open(kronareport,"w")
	for classname in classification.keys():
		kronareportfile.write(str(classification[classname]) + "\t" + classname + "\n")
	kronareportfile.close()	
	#create kronahtml
	command="ImportText.pl " + kronareport + " -o " + kronahtml
	#print(command)
	os.system(command)
	os.system("firefox " + kronahtml) 
if __name__ == "__main__":
	if prefix=="" or prefix==None:
		prefix=GetBase(predictionfilename)
		if "/" in prefix:
			prefix=prefix[prefix.rindex("/")+1:]	
	outputname=GetWorkingBase(prefix) + ".verified"
	notverifiedoutputname=GetWorkingBase(prefix) + ".notverified"
	classificationreportfilename=GetWorkingBase(prefix) + ".verified.classification"
	if outputname==predictionfilename:
		outputname=outputname+".verified"
	if verifyingrank!="":
		outputname=GetWorkingBase(prefix) + "." + verifyingrank + ".verified"
		notverifiedoutputname=GetWorkingBase(prefix) + "." + verifyingrank + ".notverified"
		classificationreportfilename=GetWorkingBase(prefix) + "." + verifyingrank + ".verified.classification"
#	if outputname==predictionfilename:
#		outputname=outputname+".verified"
	#load prediction
	predictiondict=LoadPrediction(predictionfilename)	
	#load sequences
	seqrecords={}
	if os.path.exists(fastafilename):
		seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	#load reference sequences:	
	refseqrecords={}
	if os.path.exists(referencefastafilename):
		refseqrecords=SeqIO.to_dict(SeqIO.parse(referencefastafilename, "fasta"))
	refclassificationdict,refclasses= LoadClassification(refseqrecords,classificationfilename)
	#verifying...	
	count,notree_count,total=Verify(seqrecords,predictiondict,refclasses,maxseqno,verifyingrank)
	if total >0:
		print("Number of classified sequences: " + str(total))
		print("Number of verified sequences: " + str(count) + "(" + str(round(count*100/total,2)) + " %).")
		print("Number of sequences without alignment: " + str(notree_count) + "(" + str(round(notree_count*100/total,2)) + " %).")
		print("Number of sequences with failt verification: " + str(total - count - notree_count) + "(" + str(round((total - count - notree_count)*100/total,2)) + " %).")
	#print("The results are saved in file  " + outputname)
	SaveVerification(predictiondict,outputname,notverifiedoutputname,classificationfilename)
	print("The results are saved in file  " + outputname + ", " + notverifiedoutputname + " and " + classificationreportfilename + ".")
	#making krona report
	kronareport = GetBase(outputname) + ".krona.report"
	kronahtml=GetBase(kronareport) + ".html"
	classificationdict= LoadClassificationForKronaReport(outputname)
	KronaPieCharts(classificationdict,kronareport,kronahtml)
	print("The krona report and html are saved in files " + kronareport + " and " + kronahtml + ".") 
