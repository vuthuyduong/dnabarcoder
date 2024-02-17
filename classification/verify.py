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
import json
#from sklearn.metrics import precision_recall_fscore_support
#from sklearn.metrics import cohen_kappa_score
#from sklearn.metrics import matthews_corrcoef
#from sklearn.metrics import confusion_matrix
#from sklearn.metrics import accuracy_score
#import json
from Bio import SeqIO
from Bio import Phylo
#import pylab
import matplotlib.pyplot as plt
plt.rc('font',size=4)
import random
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='verify.py',  
							   usage="%(prog)s [options] -i classified file -pred_columnname Prediction -f the fasta file -r referencefastafilename -c classificationfile -o output",
							   description='''Script that verifies the classified sequences of the prediction file to their BLAST best match based on the given cutoffs or phylogenetic trees.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the classified file')
parser.add_argument('-idcolumnname','--idcolumnname', default="ID", help='the column name of sequence id in the classification file.')
parser.add_argument('-f','--fasta', help='the fasta file')
parser.add_argument('-r','--reference', help='the reference fasta file')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', default="", help='the classification file in tab. format.')
parser.add_argument('-maxseqno','--maxseqno', type=int, default=500, help='Maximum number of the sequences of the predicted taxon name from the classification file will be selected for the comparison to find the best match. If it is not given, all the sequences will be selected.')
parser.add_argument('-pred_columnname','--prediction_columnname', default="prediction", help='the colunm name of the prediction.')
#parser.add_argument('-rank','--classificationrank', default="", help='the classification rank')
#parser.add_argument('-fullclassification_columnname','--fullclassification_columnname', default="full classification", help='the colunm name of the predicted full classification.')
#parser.add_argument('-givenlabel_columnname','--givenlabel_columnname', default="given label", help='the colunm name of the given labels, for comparison purpose.')
parser.add_argument('-redo','--redo', default="", help='the classification rank')
parser.add_argument('-prefix','--prefix', help='the prefix of output filenames')
parser.add_argument('-savefig','--savefig', default="no", help='save the figures of the phylogenetic trees or not: yes or no.')
parser.add_argument('-display','--display',default="", help='If display=="yes" then the krona html is displayed.')
parser.add_argument('-cutoff','--globalcutoff', type=float, default=0,help='The global cutoff to assign the sequences to predicted taxa. If the cutoffs file is not given, this value will be taken for sequence assignment.')
parser.add_argument('-confidence','--globalconfidence', type=float,default=0,help='The global confidence to assign the sequences to predicted taxa')
parser.add_argument('-cutoffs','--cutoffs', help='The json file containing the cutoffs to assign the sequences to the predicted taxa.')
parser.add_argument('-minseqno','--minseqno', type=int, default=0, help='the minimum number of sequences for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
parser.add_argument('-mingroupno','--mingroupno', type=int, default=0, help='the minimum number of groups for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
parser.add_argument('-minproba','--minproba', type=float, default=0, help='The minimum probability for verifying the classification results.')
parser.add_argument('-ml','--minalignmentlength', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should probably be set to smaller, 50 for instance.')
parser.add_argument('-alignmentmethod','--alignmentmethod',default="mafft", help='the alignment method: mafft or clustalo.')
parser.add_argument('-saveverifiedonly','--saveverifiedonly',default="yes", help='The option to save only verified sequences (yes) or all (no) in the classification output.')
parser.add_argument('-method','--method', default="cutoff", help='The methods (cutoff,tree) based on the similarity cutoffs or phylogenic trees for the verification of classification.')
parser.add_argument('-seqid','--sequenceid', default="", help='If the sequence id is given, then only classification of this sequence is verified. Otherwise all classifications are verified.')


args=parser.parse_args()
predictionfilename=args.input
fastafilename= args.fasta
referencefastafilename= args.reference
classificationfilename=args.classification
maxseqno=args.maxseqno
mincoverage = args.minalignmentlength
prefix=args.prefix
#verifyingrank=args.classificationrank
minproba=args.minproba
method=args.method
globalcutoff=args.globalcutoff
globalconfidence=args.globalconfidence
cutoffsfilename=args.cutoffs
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

def LoadClassification(seqrecords,classificationfilename,idcolumnname):
	classificationdict={}
	classes={}
	taxonomy={}
	isError=False
	if classificationfilename != "":
		classificationfile= open(classificationfilename, errors='ignore')
		header=next(classificationfile)
		texts=header.rstrip().split("\t")
		i=0
		for text in texts:
			if text.lower()==idcolumnname.lower():
				seqidpos=i
				i=i+1
		if 	seqidpos==-1:
			print("Please specify the sequence id columnname by using -idcolumnname.")
			isError=True
		for line in classificationfile:
			texts=line.split("\t")
			seqid = texts[seqidpos].rstrip()
			classification,taxonname,rank=GetTaxonomicClassification(0,header,texts)
			if classification!="":
				classificationdict.setdefault(seqid,{})
				classificationdict[seqid]["classification"]=classification
				classificationdict[seqid]["taxonname"]=taxonname
				classificationdict[seqid]["rank"]=rank	
			if seqid in seqrecords.keys():
				taxonnames=classification.split(";")
				ranks=["kingdom","phylum","class","order","family","genus","species"]
				level=0
				for taxonname in taxonnames:
					if taxonname=="":
						level=level+1
						continue
					taxonname=taxonname.split("__")[1]
					taxonname=taxonname.replace("_"," ")
					if taxonname=="unidentified" or taxonname=="":
						level=level+1
						continue
					rank=ranks[level]
					if not taxonname in classes.keys():
						classes.setdefault(taxonname,{})
						taxonomy.setdefault(taxonname,{})
					classes[taxonname][seqid]=seqrecords[seqid]
					taxonomy[taxonname]["rank"]=rank
					taxonomy[taxonname]["classification"]=GetRankClassification(level, classification)
					level=level+1
		classificationfile.close()	
	else:
		for seqid in seqrecords.keys():
			description=seqrecords[seqid].description
			species="s__unidentified"
			genus="g__unidentified"
			family="f__unidentified"
			order="o__unidentified"
			bioclass="c__unidentified"
			phylum="p__unidentified"
			kingdom="k__unidentified"
			if " " in description:
				description=description.split(" ")[1]
			texts=description.split("|")
			for text in texts:
				taxa=text.split(";")	
				for taxon in taxa:
					if taxon.startswith("k__"):
						kingdom=taxon
					elif taxon.startswith("p__"):
						phylum=taxon
					elif taxon.startswith("c__"):
						bioclass=taxon
					elif taxon.startswith("o__"):
						order=taxon
					elif taxon.startswith("f__"):
						family=taxon
					elif taxon.startswith("g__"):
						genus=taxon
					elif taxon.startswith("s__") and (" " in taxon.replace("s__","") or "_" in taxon.replace("s__","")):
						species=taxon	
			classification=kingdom + ";" + phylum + ";" + bioclass + ";" + order + ";" + family + ";" + genus + ";" + species
			taxonnames=[kingdom,phylum,bioclass,order,family,genus,species]
			ranks=["kingdom","phylum","class","order","family","genus","species"]
			level=0
			currenttaxonname=""
			for taxonname in taxonnames:
				if taxonname=="":
					level=level+1
					continue
				taxonname=taxonname.split("__")[1]
				taxonname=taxonname.replace("_"," ")		
				if taxonname=="unidentified" or taxonname=="":
					level=level+1
					continue
				rank=ranks[level]
				currenttaxonname=taxonname
				if not (taxonname in classes.keys()):
					classes.setdefault(taxonname,{})
					taxonomy.setdefault(taxonname,{})
				classes[taxonname][seqid]=seqrecords[seqid]
				taxonomy[taxonname]["rank"]=rank
				taxonomy[taxonname]["classification"]=GetRankClassification(level, classification)
				level=level+1	
			if currenttaxonname!="":
				classificationdict.setdefault(seqid,{})
				classificationdict[seqid]["classification"]=classification
				classificationdict[seqid]["taxonname"]=currenttaxonname
				classificationdict[seqid]["rank"]=rank	
	return classificationdict,classes,taxonomy,isError

def LoadPrediction(predictionfilename,idcolumnname,givenseqid):
	isError=False
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
	#p_rank=-1
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
		if header.lower()==idcolumnname.lower():
			p_id=i
		#if "given label" in header.lower():
		#	p_l=i
		if header.lower()==args.prediction_columnname:
			p_pl=i
		#if "full classification" in header.lower():
		#	p_c=i
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
		#if "rank" in header.lower():
		#	p_rank=i
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
		if args.prediction_columnname.lower() == header.lower():
			p_pl=i			   
		i=i+1	
	if p_id==-1:
		print("Please specify the sequence id columnname by using -idcolumnname.")
		isError=True
	if 	p_pl==-1:
		print("Please use columnname -pred_columnname for loading prediction.")
		isError=True
		
	for line in predictionfile:
		texts=line.rstrip().split("\t")
		seqid=""
		if p_id >=0 and p_id < len(texts):
			seqid=texts[p_id]
			if givenseqid!="":
				if givenseqid!=seqid:
					continue
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
		#rank=""
		#if p_rank>0 and p_rank<len(texts):
		#	rank=texts[p_rank]
		#predictiondict[seqid]["rank"]=rank
		cutoff=0
		if p_cutoff>0 and p_cutoff<len(texts):
			cutoff=texts[p_cutoff]
		predictiondict[seqid]["cut-off"]=cutoff
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
		if givenseqid==seqid:
			break
	return predictiondict,isError

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
	figfilename=treefilename + ".png"
	if (not os.path.exists(figfilename)) or redo!="":
		tree = Phylo.read(treefilename, "newick")
		Phylo.draw(tree,do_show=False)	
		#Phylo.draw_ascii(tree,do_show=False)	
		#pylab.savefig(svgfilename,format='svg', bbox_inches='tight', dpi=300)
		plt.tight_layout()
		plt.savefig(figfilename,dpi=500,bbox_inches='tight')
		print("A figure of the tree in png format is saved in  file " + figfilename + ".")	

def CreateAlignment(fastafilename,alignmentmethod):
	alignmentfilename=GetBase(fastafilename) +"." + alignmentmethod + ".aligned.fas"
	#make the alignment
	command=""
	if not os.path.exists(alignmentfilename):
		if alignmentmethod.lower()=="clustalo":
			command="clustalo -i " + fastafilename + " -o " + alignmentfilename
		else:
			command="mafft " + fastafilename + " > " + alignmentfilename
		print(command)
		os.system(command)
	return alignmentfilename

def CreateTree(fastafilename,alignmentmethod,redo):
	alignmentfilename = CreateAlignment(fastafilename,alignmentmethod)
	treefilename= alignmentfilename + ".treefile"
	if (not os.path.exists(treefilename)) or redo!="":
		#make tree
		command="iqtree -pers 0.2 -n 500 -s " + alignmentfilename
		if redo !="":
			command = "iqtree -pers 0.2 -n 500 -s " + alignmentfilename + " -redo"
		#command="iqtree -s " + alignmentfilename
		os.system(command)
		print("A iq-tree in newick format is saved in file " + treefilename + ".")	
	#print tree
	if args.savefig=="yes":
		PrintTree(treefilename,redo)
	return treefilename	

def CreateFastaFileForTrees(seqrecord,taxonname,sequences,maxseqno,redo):
	if not os.path.exists(outputpath + "/verification"):
		os.system("mkdir " + outputpath + "/verification")
	if sys.version_info[0] < 3:
		taxonname=unicode(taxonname,errors='ignore')
	newfastafilename=outputpath + "/verification/" + seqrecord.id.replace("|","_") + ".fasta"
	numberofrefsequences=len(sequences.keys())
	if not os.path.exists(newfastafilename) or redo!="":
		seqrecords=[]
		seqrecords.append(seqrecord)
		numberofrefsequences=len(sequences.keys())
		if len(sequences) >=2:#only make tree of more than 3 sequences
			if (maxseqno >0) and (len(sequences) > maxseqno):
				#select randomly 100 sequences to compare
				selectedlist=random.sample(range(0, len(sequences)), k=maxseqno)
				seqids=list(sequences.keys())
				for i in selectedlist:
					sequenceid=seqids[i]
					if sequenceid != seqrecord.id:
						seqrecords.append(sequences[sequenceid])
			else:	
				for sequenceid in sequences.keys():
					if sequenceid != seqrecord.id:
						seqrecords.append(sequences[sequenceid])
		if len(seqrecords) >=3:			
			SeqIO.write(seqrecords,newfastafilename,"fasta")
			numberofrefsequences=len(seqrecords)
		else:
			newfastafilename=""
	return newfastafilename,numberofrefsequences

def CreateFastaFileForBLAST(seqrecord,taxonname,sequences,maxseqno):
	if not os.path.exists(outputpath + "/verification"):
		os.system("mkdir " + outputpath + "/verification")
	if sys.version_info[0] < 3:
		taxonname=unicode(taxonname,errors='ignore')
	newfastafilename=""
	numberofrefsequences=0
	newfastafilename=outputpath + "/verification/" + taxonname.replace("|","_").replace(" ","_") + ".fasta"
	if not os.path.exists(newfastafilename):
		seqrecords=[]
		numberofrefsequences=len(sequences.keys())
		if (maxseqno >0) and (len(sequences) > maxseqno):
				#select randomly 100 sequences to compare
			selectedlist=random.sample(range(0, len(sequences)), k=maxseqno)
			seqids=list(sequences.keys())
			for i in selectedlist:
				sequenceid=seqids[i]
				seqrecords.append(sequences[sequenceid])
		else:	
			for sequenceid in sequences.keys():
				seqrecords.append(sequences[sequenceid])
	if len(seqrecords) >=1:			
		SeqIO.write(seqrecords,newfastafilename,"fasta")
	else:
		newfastafilename=""			
	return newfastafilename,numberofrefsequences

# def GetLevel(rank):
# 	level=-1
# 	if rank=="species":	
# 		level=6
# 	elif rank=="genus":	
# 		level=5
# 	elif rank=="family":	
# 		level=4
# 	elif rank=="order":	
# 		level=3
# 	elif rank=="class":	
# 		level=2
# 	elif rank=="phylum":
# 		level=1
# 	elif rank=="kingdom":
# 		level=0	
# 	return level

def GetHigherTaxa(rank,classification):
	highertaxa=[]
	taxa=classification.split(";")
	species=""
	genus=""
	family=""
	order=""
	bioclass=""
	phylum=""
	kingdom=""
	for taxon in taxa:
		if taxon.startswith("k__"):
			kingdom=taxon.replace("k__","")
		if taxon.startswith("p__"):
			phylum=taxon.replace("p__","")	
		if taxon.startswith("c__"):
			bioclass=taxon.replace("c__","")	
		if taxon.startswith("o__"):
			order=taxon.replace("o__","")	
		if taxon.startswith("f__"):
			family=taxon.replace("f__","")	
		if taxon.startswith("g__"):
			genus=taxon.replace("g__","")	
		if taxon.startswith("s__"):
			species=taxon.replace("s__","")	
	level=GetLevel(rank)	
	if level >=6 and species!="" and species!="unidentified":
		highertaxa.append(species)
	if level >=5 and genus!="" and genus!="unidentified":
		highertaxa.append(genus)
	if level >=4 and family!="" and family!="unidentified":
		highertaxa.append(family)	
	if level >=3 and order!="" and order!="unidentified":
		highertaxa.append(order)	
	if level >=2 and bioclass!="" and bioclass!="unidentified":
		highertaxa.append(bioclass)	
	if level >=1 and phylum!="" and phylum!="unidentified":
		highertaxa.append(phylum)
	if level >=0 and kingdom!="" and kingdom!="unidentified":
		highertaxa.append(kingdom)
	return highertaxa

def GetCutoffAndConfidence(rank,classification,cutoffs):
	if not rank in cutoffs.keys():
		return [0,0,False]
	#cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","")
	#taxa=cleanclassification.split(";")
	#taxa.append("All") 
	highertaxa=GetHigherTaxa(rank,classification)
	highertaxa.append("All") 
	localcutoff=0
	seqno=0
	groupno=0
	datasets=cutoffs[rank]
	maxconfidence=-1
	bestcutoff=0
	isComputed=False
	for highertaxonname in highertaxa:
		if not highertaxonname in datasets.keys():
			continue
		if "cut-off" in datasets[highertaxonname].keys():
			localcutoff=datasets[highertaxonname]["cut-off"]
			isComputed=True
		confidence=0	
		if "confidence" in datasets[highertaxonname].keys():
			confidence=datasets[highertaxonname]["confidence"]
		seqno=0
		if "sequence number" in datasets[highertaxonname].keys():
			seqno=datasets[highertaxonname]["sequence number"]	
		groupno=0
		if "group number" in datasets[highertaxonname].keys():
			groupno=datasets[highertaxonname]["group number"]	
		if not ((seqno >0 and seqno < args.minseqno) or (groupno >0 and groupno < args.mingroupno)):	
			if maxconfidence < confidence:
				maxconfidence =confidence
				bestcutoff=localcutoff
			if isComputed==True:	
				break
	return [bestcutoff,maxconfidence,isComputed]

def GetRank(taxonname,classification):
	rank=""
	level=-1
	texts=classification.split(";")
	for text in texts:
		if not taxonname in text:
			continue
		if text.startswith("k__"):
			rank="kingdom"
			level=0
		if text.startswith("p__"):
			rank="phylum"
			level=1
		if text.startswith("c__"):
			rank="class"	
			level=2
		if text.startswith("o__"):
			rank="order"	
			level=3
		if text.startswith("f__"):
			rank="family"
			level=4
		if text.startswith("g__"):
			rank="genus"
			level=5
		if text.startswith("s__"):
			rank="species"	
			level=6
	return rank,level

def GetRankClassification(level,classification):
	if classification=="" or level ==-1:
		return "k__unidentified;p__unidentified;c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified"	
	species=classification.split(";")[6].replace("s__","")
	genus=classification.split(";")[5].replace("g__","")
	family=classification.split(";")[4].replace("f__","")
	order=classification.split(";")[3].replace("o__","")
	bioclass=classification.split(";")[2].replace("c__","")
	phylum=classification.split(";")[1].replace("p__","")
	kingdom=classification.split(";")[0].replace("k__","")
	newclassification=""
	if level >=0 and kingdom!="unidentified":
		newclassification="k__" + kingdom +";p__unidentified;c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified"	
	if level >=1 and phylum!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified"
	if level >=2 and bioclass!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__"+bioclass+";o__unidentified;f__unidentified;g__unidentified;s__unidentified"
	if level >=3 and order!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order + ";f__unidentified;g__unidentified;s__unidentified"
	if level >=4 and family!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family +";g__unidentified;s__unidentified"
	if level >=5 and genus!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family + ";g__"+ genus + ";s__unidentified"
	if level >=6 and species!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family + ";g__"+ genus+";s__"+species 	
	return newclassification

def AddCutoffsToTaxonomy(taxonomy,cutoff,confidence,cutoffs):
	for taxonname in taxonomy.keys():
		if cutoffs!={}:
			classification=taxonomy[taxonname]["classification"]
			rank=taxonomy[taxonname]["rank"]
			if taxonname in cutoffs.keys():
				taxonomy[taxonname]["cut-off"]=cutoffs[taxonname]["cut-off"]
				taxonomy[taxonname]["confidence"]=cutoffs[taxonname]["confidence"]
			else:	
				cutoff_confidence=GetCutoffAndConfidence(rank,classification,cutoffs)
				if cutoff_confidence[2]==True:
					taxonomy[taxonname]["cut-off"]=cutoff_confidence[0]
					taxonomy[taxonname]["confidence"]=cutoff_confidence[1]
			if not ("cut-off" in taxonomy[taxonname].keys()) and globalcutoff >=0: #use the globalcutoff
				taxonomy[taxonname]["cut-off"]=globalcutoff
				taxonomy[taxonname]["confidence"]=globalconfidence
		else:
			taxonomy[taxonname]["cut-off"]=cutoff
			taxonomy[taxonname]["confidence"]=confidence

def ComputeBestLocalBLASTScore(testrecord,reffilename,mincoverage):
	#Create fasta file of the test record 
	queryname=testrecord.id + "." + "fasta"
	if "|" in testrecord.id:
		queryname=testrecord.id.split("|")[0] + "." + "fasta"
	SeqIO.write(testrecord, queryname , "fasta")
	#Create fasta file of predictedname
	if reffilename=="":
		return "",0,0,0
	#Blast the test record to fastafilename:
	#makedbcommand = "makeblastdb -in " + refname + " -dbtype \'nucl\' " +  " -out db"
	#os.system(makedbcommand)
	blastoutput = testrecord.id 
	if "|" in testrecord.id:
		blastoutput=blastoutput.split("|")[0]
	#blast
	makedbcommand = "makeblastdb -in " + reffilename + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	
	blastoutput= blastoutput + ".blast.out"
	#blastcommand = "blastn -query " + queryname + " -db  db -outfmt 6 -out " + blastoutput
	blastcommand = "blastn -query " + queryname + " -db  db -task blastn-short -outfmt 6 -out " + blastoutput + " -num_threads " + str(nproc)
	#for long read
	if mincoverage >=300:
		blastcommand = "blastn -query " + queryname + " -db  db -outfmt 6 -out " + blastoutput + " -num_threads " + str(nproc)
	os.system(blastcommand)
	#read output of Blast
	blastoutputfile = open(blastoutput)
	bestlocalscore=0
	bestlocalsim=0
	bestlocalcoverage=0
	bestrefid=""
	for line in blastoutputfile:
		words = line.split("\t")
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		refid=words[1]
		score=sim
		if coverage < mincoverage:
				score=float(score * coverage)/mincoverage
		if score > bestlocalscore:
			bestrefid=refid
			bestlocalscore=score
			bestlocalsim=sim
			bestlocalcoverage=coverage
	os.system("rm " + blastoutput)
	os.system("rm " + queryname)
	return bestrefid,bestlocalscore,bestlocalsim,bestlocalcoverage

def VerifyBasedOnCutoffs(seqrecords,predictiondict,refclasses,maxseqno,taxonomy,redo):
	count=0
	total=0
	#create Fasta files
	for seqid in predictiondict.keys():
		predictedname=predictiondict[seqid]["predlabel"]
		#rank=predictiondict[seqid]["rank"]
		refid=predictiondict[seqid]["refid"]
		coverage=predictiondict[seqid]["coverage"]
		sim=predictiondict[seqid]["sim"]
		score=predictiondict[seqid]["score"]
		proba=predictiondict[seqid]["proba"]
		numberofrefsequences=0
		if score==0:
			score=sim
		verifiedlabel=""
		#if (verifyingrank=="" or (verifyingrank!="" and rank==verifyingrank)) and (proba >=minproba): #verifying only for classification results with probability greater than min proba
		if proba >=minproba:			
		#only predict when the tree file name does not exist
			if predictedname in refclasses.keys() and seqid in seqrecords.keys():
				total=total+1
				seqrecord=seqrecords[seqid]
				if score==0 or (redo !=""):
					sequences=refclasses[predictedname]
					reffilename,numberofrefsequences=CreateFastaFileForBLAST(seqrecord,predictedname,sequences,maxseqno)
					if os.path.exists(reffilename):
						#compute BLAST score
						newrefid,newbestscore,newsim,newcoverage=ComputeBestLocalBLASTScore(seqrecord,reffilename,mincoverage)
						os.system("rm " + reffilename)	
						if newbestscore > score:
							refid=newrefid
							score=newbestscore
							sim=newsim
							coverage=newcoverage
				predicted_cutoff=taxonomy[predictedname]["cut-off"]
				predicted_confidence=taxonomy[predictedname]["confidence"]
				predicted_classification=taxonomy[predictedname]["classification"]
				if 	score >= predicted_cutoff:
					verifiedlabel=predictedname
					count=count+1
					predictiondict[seqid]["verifiedlabel"]=verifiedlabel
				predictiondict[seqid]["numberofrefsequences"]=numberofrefsequences
				predictiondict[seqid]["classification"]=predicted_classification
				predictiondict[seqid]["refid"]=refid
				predictiondict[seqid]["sim"]=sim
				predictiondict[seqid]["score"]=score
				predictiondict[seqid]["coverage"]=coverage
				predictiondict[seqid]["cut-off"]=predicted_cutoff
				predictiondict[seqid]["confidence"]=predicted_confidence
	return count,total
	
def VerifyBasedOnTrees(seqrecords,predictiondict,refclasses,maxseqno,alignmentmethod,redo):
	count=0
	notree_count=0
	total=0
	#create Fasta files
	for seqid in predictiondict.keys():
		predictedname=predictiondict[seqid]["predlabel"]
		#rank=predictiondict[seqid]["rank"]
		treefilename=predictiondict[seqid]["treefilename"]
		numberofrefsequences=predictiondict[seqid]["numberofrefsequences"]
		branchlength=predictiondict[seqid]["branchlength"]
		maxbranchlength=predictiondict[seqid]["maxbranchlength"]
		averagebranchlength=predictiondict[seqid]["averagebranchlength"]
		proba=predictiondict[seqid]["proba"]
		#if sys.version_info[0] < 3:
			#predictedname=unicode(predictedname,'latin1')
		verifiedlabel=""
		#if (verifyingrank=="" or (verifyingrank!="" and rank==verifyingrank)) and (proba >=minproba):#verifying only for classification results with probability greater than min proba
		if proba >=minproba:#verifying only for classification results with probability greater than min proba
			#only predict when the tree file name does not exist
			if (predictedname in refclasses.keys()) and (seqid in seqrecords.keys()):
				total=total+1
				seqrecord=seqrecords[seqid]
				verified = False
				if treefilename=="" or (redo !=""):
					sequences=refclasses[predictedname]
					if not (seqid in sequences.keys()):
						treefastafilename,numberofrefsequences=CreateFastaFileForTrees(seqrecord,predictedname,sequences,maxseqno,redo)
						if os.path.exists(treefastafilename):
							treefilename= CreateTree(treefastafilename,alignmentmethod,redo)
					else:
						verified=True
						print("The sequence " + seqid + " is in the reference file. No verification needed.")
				if verified==False:
					if 	os.path.exists(treefilename):
						verified,branchlength,maxbranchlength,averagebranchlength=verifyBasedOnBranchLengths(seqid,treefilename)
				if verified==True:
					verifiedlabel=predictedname
					count=count+1
				else:
					notree_count = notree_count + 1
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
	notverifiedoutputfile.write("ID\tGiven label\tPrediction\tFull classification\tProbability\tCut-off\tConfidence\tReferenceID\tBLAST score\tBLAST sim\tBLAST coverage\tVerified label\tTree filename\tNumber of reference sequences\tBranch length\tMax branch length\tAverage branch length\n")
	outputfile.write("ID\tGiven label\tPrediction\tFull classification\tProbability\tCut-off\tConfidence\tReferenceID\tBLAST score\tBLAST sim\tBLAST coverage\tVerified label\tTree filename\tNumber of reference sequences\tBranch length\tMax branch length\tAverage branch length\n")
	classificationreportfile.write("SequenceID\tReferenceID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tscore\tcutoff\tconfidence\n")
	for seqid in predictiondict.keys():
		prediction=predictiondict[seqid]
		treefilename=prediction["treefilename"]
		predlabel=prediction["predlabel"]
		givenlabel=prediction["givenlabel"]
		classification=prediction["classification"]
		#rank=prediction["rank"]
		bestscore=prediction["score"]
		cutoff=prediction["cut-off"]
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
		if args.saveverifiedonly!="yes":
			outputfile.write(seqid + "\t" + givenlabel + "\t"  + verifiedlabel + "\t"+ classification + "\t" + str(proba) + "\t" + str(cutoff) + "\t" + str(confidence) + "\t" + refid + "\t" + str(bestscore) + "\t" + str(sim) + "\t" + str(coverage) + "\t" + verifiedlabel + "\t" + treefilename + "\t" + str(numberofrefsequences) + "\t" + str(branchlength) + "\t" + str(maxbranchlength) + "\t" + str(averagebranchlength) + "\n")			
			#cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","").replace("_"," ")
		else:
			if verifiedlabel!="" and verifiedlabel!="unidentified":			
				outputfile.write(seqid + "\t" + givenlabel + "\t"  + predlabel + "\t"+ classification + "\t" + str(proba) + "\t" + str(cutoff) + "\t" + str(confidence) + "\t" + refid + "\t" + str(bestscore) + "\t" + str(sim) + "\t" + str(coverage) + "\t" + verifiedlabel + "\t" + treefilename + "\t" + str(numberofrefsequences) + "\t" + str(branchlength) + "\t" + str(maxbranchlength) + "\t" + str(averagebranchlength) + "\n")			
				#cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","").replace("_"," ")
		#save un verified sequences
		#if (verifiedlabel=="" or verifiedlabel=="unidentified") and ((verifyingrank!="" and verifyingrank==rank) or verifyingrank==""):
		if (verifiedlabel=="" or verifiedlabel=="unidentified"):
			notverifiedoutputfile.write(seqid + "\t" + givenlabel + "\t"  + predlabel + "\t"+ classification + "\t" + str(proba)  + "\t" + str(cutoff) + "\t" + str(confidence) + "\t" + refid + "\t" + str(bestscore) + "\t" + str(sim) + "\t" + str(coverage) + "\t" + verifiedlabel + "\t" + treefilename + "\t" + str(numberofrefsequences) + "\t" + str(branchlength) + "\t" + str(maxbranchlength) + "\t" + str(averagebranchlength) + "\n")			
		classificationreportfile.close()
	outputfile.close()	
	notverifiedoutputfile.close()	
		
def LoadClassificationForKronaReport(classificationfilename):
	classificationdict={}
	classificationfile= open(classificationfilename)
	next(classificationfile)
	for line in classificationfile:
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
	if args.display=="yes":
		os.system("firefox " + kronahtml) 
if __name__ == "__main__":
	if prefix=="" or prefix==None:
		prefix=GetBase(predictionfilename)
		if "/" in prefix:
			prefix=prefix[prefix.rindex("/")+1:]	
	outputname=GetWorkingBase(prefix) + ".verified"
	notverifiedoutputname=GetWorkingBase(prefix) + ".unverified"
	classificationreportfilename=GetWorkingBase(prefix) + ".verified.classification"
	if outputname==predictionfilename:
		outputname=outputname+".verified"
# 	if verifyingrank!="":
# 		outputname=GetWorkingBase(prefix) + "." + verifyingrank + ".verified"
# 		notverifiedoutputname=GetWorkingBase(prefix) + "." + verifyingrank + ".unverified"
# 		classificationreportfilename=GetWorkingBase(prefix) + "." + verifyingrank + ".verified.classification"
	#load prediction
	predictiondict,isError=LoadPrediction(predictionfilename,args.idcolumnname,args.sequenceid)	
	if isError==True:
		sys.exit()
	#load sequences
	seqrecords={}
	if os.path.exists(fastafilename):
		seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	#load reference sequences:	
	refseqrecords={}
	if os.path.exists(referencefastafilename):
		refseqrecords=SeqIO.to_dict(SeqIO.parse(referencefastafilename, "fasta"))
	refclassificationdict,refclasses,taxonomy,isError= LoadClassification(refseqrecords,classificationfilename,args.idcolumnname)
	if isError==True or refclasses=={} or refclassificationdict=={}:
		print("Please check the classification of reference sequences.")
		sys.exit()
	#verifying...	
	if args.method=="tree":
		count,notree_count,total=VerifyBasedOnTrees(seqrecords,predictiondict,refclasses,maxseqno,args.alignmentmethod,args.redo)
		if total >0:
			print("Number of classified sequences: " + str(total))
			print("Number of verified sequences: " + str(count) + "(" + str(round(count*100/total,2)) + " %).")
			print("Number of sequences without alignment: " + str(notree_count) + "(" + str(round(notree_count*100/total,2)) + " %).")
			print("Number of sequences with failt verification: " + str(total - count - notree_count) + "(" + str(round((total - count - notree_count)*100/total,2)) + " %).")
			if notree_count==0:
				print("Please check if the trees have been created, or mafft/clustalo and/or iqtree have been installed.")
	else:
		cutoffs={}
		if cutoffsfilename!="" and cutoffsfilename!=None:
			with open(cutoffsfilename) as cutoffsfile:
				cutoffs = json.load(cutoffsfile)
		#add cutoffs to taxa for sequence identification		
		AddCutoffsToTaxonomy(taxonomy,globalcutoff,globalconfidence,cutoffs)		
		count,total=VerifyBasedOnCutoffs(seqrecords,predictiondict,refclasses,maxseqno,taxonomy,args.redo)
		if total >0:
			print("Number of classified sequences: " + str(total))
			print("Number of verified sequences: " + str(count) + "(" + str(round(count*100/total,2)) + " %).")	
	#print("The results are saved in file  " + outputname)
	SaveVerification(predictiondict,outputname,notverifiedoutputname,classificationfilename)
	print("The results are saved in file  " + outputname + ", " + notverifiedoutputname + " and " + classificationreportfilename + ".")
	#making krona report
	kronareport = GetBase(outputname) + ".krona.report"
	kronahtml=GetBase(kronareport) + ".html"
	classificationdict= LoadClassificationForKronaReport(outputname)
	KronaPieCharts(classificationdict,kronareport,kronahtml)
	print("The krona report and html are saved in files " + kronareport + " and " + kronahtml + ".") 
