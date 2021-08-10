#!/usr/bin/env python
# FILE: classify.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
import sys
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
from Bio import SeqIO
import json
import multiprocessing

nproc=multiprocessing.cpu_count()
#from keras.utils import np_utils

parser=argparse.ArgumentParser(prog='classify.py',  
							   usage="%(prog)s [options] -i fastafile -r referencefastafile -t optthreshold -c classificationfilename -p classificationposition -mp minproba -mc mincoverage ",
							   description='''Script that classifies the sequences of the fasta files using BLAST with a cut-off value or the cut-off values given in the cutoffs file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file to be classified.')
parser.add_argument('-r','--reference', required=True, help='the reference fasta file.')
parser.add_argument('-mc','--mincoverage', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, mc should probably be set to 100.')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', default="", help='the classification file in tab. format.')#optinal
parser.add_argument('-cutoffs','--cutoffs', help='The json file containing the cutoffs to assign the sequences to the predicted taxa.')
parser.add_argument('-cutoff','--cutoff', type=float, default=0,help='The cutoff to assign the sequences to predicted taxa. If the cutoffs file is not given, this value will be taken for sequence assignment.')
parser.add_argument('-confidence','--confidence', type=float,default=0.7,help='The confidence of the cutoff to assign the sequences to predicted taxa')
parser.add_argument('-rank','--classificationrank', default="", help='the classification rank')
parser.add_argument('-prefix','--prefix', help='the prefix of output filenames')

args=parser.parse_args()
testdataset= args.input
traindataset = args.reference
mincoverage = args.mincoverage
classificationfilename=args.classification
classificationrank=args.classificationrank
cutoff=args.cutoff
globalconfidence=args.confidence
cutoffsfilename=args.cutoffs
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

def GetRankClassification(level,classification):
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
		classification="k__" + kingdom +";p__"+phylum +";c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified"
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
	species=""
	genus=""
	family=""
	order="" 
	bioclass=""
	phylum=""
	kingdom=""
	if p_s< len(texts):
		species=texts[p_s].rstrip()
	if p_g< len(texts):
		genus=texts[p_g].rstrip()
	if p_f< len(texts):
		family=texts[p_f].rstrip()
	if p_o< len(texts):
		order=texts[p_o].rstrip()
	if p_c< len(texts):
		bioclass=texts[p_c].rstrip()
	if p_p< len(texts):
		phylum=texts[p_p].rstrip()
	if p_k< len(texts):
		kingdom=texts[p_k].rstrip()
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
		if seqid in seqrecords.keys():
			classificationdict.setdefault(seqid,{})
			classification,taxonname,rank=GetTaxonomicClassification(0,header,texts)
			classificationdict[seqid]["classification"]=classification
			classificationdict[seqid]["taxonname"]=taxonname
			classificationdict[seqid]["rank"]=rank
			if taxonname in classes.keys():
				classes[taxonname].append(seqrecords[seqid])
			else:
				classes.setdefault(taxonname,[seqrecords[seqid]])		
	classificationfile.close()	
	return classificationdict,classes

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
	os.system(makedbcommand)
	#for short read
	blastcommand = "blastn -query " + indexed_query + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	#for long read
	if mincoverage >=300:
		blastcommand = "blastn -query " + indexed_query + " -db  db -outfmt 6 -out out.txt -num_threads " + str(nproc)
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
		if score > bestscorelist[i]:
			bestscorelist[i]= score
			bestrefidlist[i]=refid
			bestsimlist[i]=sim
			bestcoveragelist[i]=coverage
	os.system("rm " + indexed_query)		
	os.system("rm out.txt")
	return bestrefidlist,bestscorelist,bestsimlist,bestcoveragelist

def GetCutoffAndConfidence(rank,classification,cutoffs):
	if not rank in cutoffs.keys():
		return [0,0]
	cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","")
	taxa=cleanclassification.split(";")
	taxa.append("All") 
	localcutoff=0
	confidence=0
	seqno=0
	groupno=0
	datasets=cutoffs[rank]
	maxconfidence=0
	bestcutoff=0
	for highertaxonname in taxa:
		if not highertaxonname in datasets.keys():
			continue
		if "cut-off" in datasets[highertaxonname].keys():
			localcutoff=datasets[highertaxonname]["cut-off"]
		if "confidence" in datasets[highertaxonname].keys():
			confidence=datasets[highertaxonname]["confidence"]
		if "sequence number" in datasets[highertaxonname].keys():
			seqno=datasets[highertaxonname]["sequence number"]	
		if "group number" in datasets[highertaxonname].keys():
			groupno=datasets[highertaxonname]["group number"]	
		if not ((seqno >0 and seqno < args.minseqno) or (groupno >0 and groupno < args.mingroupno)):	
			if maxconfidence < confidence:
				maxconfidence =confidence
				bestcutoff=localcutoff
	return [bestcutoff,maxconfidence]

def GetCutoffs(classification,cutoffs):
	species="unindentified"
	genus="unindentified"
	family="unindentified"
	order="unindentified"
	bioclass="unindentified"
	phylum="unindentified"
	kingdom="unindentified"
	texts=classification.split(";")
	for text in texts:
		if text.startswith("s__"):
			species=text.replace("s__","")
		if text.startswith("g__"):
			genus=text.replace("g__","")	
		if text.startswith("f__"):
			family=text.replace("f__","")	
		if text.startswith("o__"):
			order=text.replace("o__","")	
		if text.startswith("c__"):
			bioclass=text.replace("c__","")	
		if text.startswith("p__"):
			phylum=text.replace("p__","")	
		if text.startswith("k__"):
			kingdom=text.replace("k__","")	
	taxacutoffs={}		
	taxacutoffs["species"]=GetCutoffAndConfidence("species",classification,cutoffs)
	taxacutoffs["genus"]=GetCutoffAndConfidence("genus",classification,cutoffs)
	taxacutoffs["family"]=GetCutoffAndConfidence("family",classification,cutoffs)
	taxacutoffs["order"]=GetCutoffAndConfidence("order",classification,cutoffs)
	taxacutoffs["class"]=GetCutoffAndConfidence("class",classification,cutoffs)
	taxacutoffs["phylum"]=GetCutoffAndConfidence("phylum",classification,cutoffs)
	taxacutoffs["kingdom"]=GetCutoffAndConfidence("kingdom",classification,cutoffs)
	return taxacutoffs,kingdom,phylum,bioclass,order,family,genus,species

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

def GetAssignment(refid,classificationdict,bestscore,cutoffs,cutoff,globalconfidence,classificationrank):
	confidence=globalconfidence
	rank=classificationrank
	localcutoff=cutoff
	taxonname=""
	level=-1
	classification=""
	if cutoffs!={} and refid!="" and (refid in classificationdict.keys()):
		refclassification=classificationdict[refid]['classification']
		taxacutoffs,kingdom,phylum,bioclass,order,family,genus,species=GetCutoffs(refclassification,cutoffs)
		if bestscore >=taxacutoffs["species"][0] and taxacutoffs["species"][0] >0 and species!="unidentified" and (classificationrank=="species" or classificationrank== ""):
			rank="species"
			localcutoff=taxacutoffs["species"][0]
			confidence=taxacutoffs["species"][1]
			taxonname=species		
			level=6
		elif bestscore >=taxacutoffs["genus"][0] and taxacutoffs["genus"][0] >0 and genus!="unidentified" and (classificationrank=="genus" or classificationrank== ""):
			rank="genus"
			localcutoff=taxacutoffs["genus"][0]
			confidence=taxacutoffs["genus"][1]
			taxonname=genus
			level=5
		elif bestscore >=taxacutoffs["family"][0] and taxacutoffs["family"][0] >0 and family!="unidentified" and (classificationrank=="family" or classificationrank== ""):
			rank="family"
			localcutoff=taxacutoffs["family"][0]
			confidence=taxacutoffs["family"][1]
			taxonname=family	
			level=4
		elif bestscore >=taxacutoffs["order"][0] and taxacutoffs["order"][0] >0 and order!="unidentified" and (classificationrank=="order" or classificationrank== ""):
			rank="order"
			localcutoff=taxacutoffs["order"][0]
			confidence=taxacutoffs["order"][1]
			taxonname=order
			level=3
		elif bestscore >=taxacutoffs["class"][0] and taxacutoffs["class"][0] >0 and bioclass!="unidentified" and (classificationrank=="class" or classificationrank== ""):
			rank="class"
			localcutoff=taxacutoffs["class"][0]
			confidence=taxacutoffs["class"][1]
			taxonname=bioclass
			level=2
		elif bestscore >=taxacutoffs["phylum"][0] and taxacutoffs["phylum"][0] >0 and phylum!="unidentified" and (classificationrank=="phylum" or classificationrank== ""):
			rank="phylum"
			localcutoff=taxacutoffs["phylum"][0]
			confidence=taxacutoffs["phylum"][1]
			taxonname=phylum
			level=1
		elif bestscore >=taxacutoffs["kingdom"][0] and taxacutoffs["kingdom"][0] >0 and kingdom!="unidentified" and (classificationrank=="kingdom" or classificationrank== ""):
			rank="kingdom"
			localcutoff=taxacutoffs["kingdom"][0]
			confidence=taxacutoffs["kingdom"][1]
			taxonname=kingdom
			level=0	
		level=GetLevel(rank)
		classification=GetRankClassification(level,refclassification)
	elif refid!="" and (refid in classificationdict.keys()):
		refclassification=classificationdict[refid]['classification']
		if bestscore >=cutoff:
			if classificationrank=="":
				taxonname=classificationdict[refid]['taxonname']
				rank=classificationdict[refid]['rank']
				level=GetLevel(rank)	
			else:	
				level=GetLevel(classificationrank)	
				rank=classificationrank
				taxonname=(refclassification.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","")).split(";")[level]
		classification=GetRankClassification(level,refclassification)		
	return classification,taxonname,rank,level,localcutoff,confidence
	
def SavePrediction(classificationdict,refseqIDs,testclassification,testseqIDs,cutoffs,cutoff,globalconfidence,bestscorelist,bestsimlist,bestcoveragelist,bestrefidlist,outputname,classificationreportfilename):
	classificationlevel=GetLevel(classificationrank)
	output=open(outputname,"w")
	classificationreportfile=open(classificationreportfilename,"w")
	output.write("SequenceID\tGiven label\tPrediction\tFull classification\tRank\tCut-off\tConfidence\tReferenceID\tBLAST score\tBLAST sim\tBLAST coverage\n")
	classificationreportfile.write("SequenceID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\trank\tscore\tcutoff\tconfidence\n")
	i=0
	count=0
	given_labels=[]
	assigned_labels=[]
	unclassifiedseqids=[]
	for seqid in testseqIDs:
		refid=bestrefidlist[i]
		bestscore=bestscorelist[i]
		classification,taxonname,rank,level,localcutoff,confidence=GetAssignment(refid,classificationdict,bestscore,cutoffs,cutoff,globalconfidence,classificationrank)	
		if taxonname!="":
			if sys.version_info[0] < 3:
				taxonname=unicode(taxonname,'latin1')
		if level>=classificationlevel and level!=-1:
			cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","")
			classificationreportfile.write(seqid + "\t" + cleanclassification.replace(";","\t") + "\t" + rank + "\t" + str(bestscore) + "\t" + str(localcutoff) + "\t" + str(confidence) + "\n")
			count=count+1
		else:
			unclassifiedseqids.append(seqid)
		testlabel=""
		if seqid in testclassificationdict.keys():
			#testlabel=testclassificationdict[seqid]['taxonname']
			testclassification=testclassificationdict[seqid]['classification']
			testclassification=testclassification.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","")
			testlabel=testclassification.split(";")[level]
		testlabel=testlabel.replace("_"," ")
		taxonname=taxonname.replace("_"," ")	
		if testlabel!="" and testlabel!="unidentified":		
			given_labels.append(testlabel)
			assigned_labels.append(taxonname)	
		output.write(seqid + "\t" + testlabel + "\t"  + taxonname + "\t" + classification + "\t" + rank + "\t" + str(localcutoff) + "\t" + str(confidence) + "\t"  + bestrefidlist[i] + "\t" +  str(bestscorelist[i]) + "\t" + str(bestsimlist[i]) + "\t" + str(bestcoveragelist[i]) +"\n")
		i=i+1
	output.close()
	classificationreportfile.close()
	return given_labels,assigned_labels,unclassifiedseqids,count 

def CalculateMetrics(test_labels,pred_labels,labels): 
        accuracy=accuracy_score(test_labels,pred_labels)
        precision,recall,fscore,support=precision_recall_fscore_support(test_labels,pred_labels,average='macro')
        precisionvector,recallvector,fscorevector,support=precision_recall_fscore_support(test_labels,pred_labels,labels=labels)
        mcc=matthews_corrcoef(test_labels,pred_labels)
        #cohenkappascore=cohen_kappa_score(test_labels,pred_labels,labels=labels)
        confusionmatrix=confusion_matrix(test_labels,pred_labels,labels=labels)
        return round(accuracy,4),round(precision,4),round(recall,4),round(fscore,4),precisionvector,recallvector,fscorevector,round(mcc,4),confusionmatrix

def LoadClassificationForKronaReport(classificationfilename):
	classificationdict={}
	classificationfile= open(classificationfilename)
	next(classificationfile)
	for line in classificationfile:
		if line.startswith("#"):
			continue 
		texts=line.split("\t")
		classname=texts[2]
		if classname=="":
			classname="unidentified"
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
	
def LoadTaxa(classificationfilename):
	taxa=[]
	classificationfile=open(classificationfilename)
	for line in classificationfile:
		line=line.rstrip()
		texts=line.split("\t")
		for text in texts:
			if text!="":
				taxa.append(text)
	classificationfile.close()			
	return taxa	
	
def CalculateClassificationMetrics(givenlabels,predlabels,reftaxa,reportname):
	tmp_givenlabels=[]
	tmp_predlabels=[]
	numberofspeciescomplexes=0
	i=0
	for label in givenlabels:
		genus=""
		predgenus=""
		predlabel=predlabels[i]
		if " " in label:
			genus=label.split(" ")[0]
			label=label.split(" ")[1]
		if " " in predlabels[i]:
			predgenus=predlabel.split(" ")[0]	
			predlabel=predlabel.split(" ")[1]
		if predgenus==genus and genus!="" and label!=predlabel:
			numberofspeciescomplexes=numberofspeciescomplexes+1
		tmp_givenlabels.append(label)		
		tmp_predlabels.append(predlabel)
		i=i+1		
	accuracy,precision,recall,fscore,precisionvector,recallvector,fscorevector,mcc,confusionmatrix=CalculateMetrics(tmp_givenlabels,tmp_predlabels,tmp_predlabels)
	filteredgivenlabels=[]
	filteredpredlabels=[]
	filterednumberofspeciescomplexes=0
	i=0
	for label in givenlabels:
		if label in reftaxa:
			genus=""
			predgenus=""
			predlabel=predlabels[i]
			if " " in label:
				genus=label.split(" ")[0]
				label=label.split(" ")[1]
			if " " in predlabels[i]:
				predgenus=predlabel.split(" ")[0]	
				predlabel=predlabel.split(" ")[1]
			if predgenus==genus and genus!="" and label!=predlabel:
				filterednumberofspeciescomplexes=filterednumberofspeciescomplexes+1	
			filteredgivenlabels.append(label)
			filteredpredlabels.append(predlabel)
		i=i+1			
	filteredaccuracy,filteredprecision,filteredrecall,filteredfscore,filteredprecisionvector,filteredrecallvector,filterdfscorevector,filteredmcc,filteredconfusionmatrix=CalculateMetrics(filteredgivenlabels,filteredpredlabels,filteredpredlabels)
	if numberofspeciescomplexes>0:
		print("Number of species complexes: " + str(numberofspeciescomplexes) + "(" +  str(round(numberofspeciescomplexes*100/len(given_labels),2)) + "%).")
	if filterednumberofspeciescomplexes>0:
		print("Number of species complexes that are present in the reference dataset: " + str(filterednumberofspeciescomplexes) + "(" +  str(round(filterednumberofspeciescomplexes*100/len(filteredgivenlabels),2)) + "%).")	
	print("Number of sequences with a given label: " + str(len(given_labels)) + "\t" + str(len(filteredgivenlabels)))
	print("The mcc of assigning the sequences with a given label: " + str(mcc) + "\t" + str(filteredmcc))
	print("The accuracy of assigning the sequence with a given label: " + str(accuracy) + "\t" + str(filteredaccuracy))
	print("The precision of assigning the sequence with a given label: " + str(precision) + "\t" + str(filteredprecision))
	print("The fscore of assigning the sequence with a given label: " + str(fscore) + "\t" + str(filteredfscore))
		
	report=open(reportname,"w")
	report.write("Number of sequences with a given label\tMcc\tAccuracy\tPrecision\tFscore\tNumber of sequences with a given label present in the references\tMcc\tAccuracy\tPrecision\tFscore\n")
	report.write(str(len(givenlabels))+ "\t" + str(mcc) + "\t" + str(accuracy) + "\t" + str(precision) + "\t" + str(fscore) + "\t")
	report.write(str(len(filteredgivenlabels))+ "\t" + str(filteredmcc) + "\t" + str(filteredaccuracy) + "\t" + str(filteredprecision) + "\t" + str(filteredfscore) + "\n")
	
	report.write("Taxonname\tAccuracy\tPrecision\tFscore\n")
	i=0
	existinglabels=[]
	for label in predlabels:
		if label in existinglabels:
			i=i+1
			continue
		existinglabels.append(label)
		report.write(label + "\t" + str(recallvector[i]) + "\t" + str(precisionvector[i]) + "\t" + str(fscorevector[i]) + "\n")
		i=i+1
	report.close()
	print("Accuracy, precision, and fscore of the taxa are given in file " + reportname +  ".")		
##############################################################################
# MAIN
##############################################################################
path=sys.argv[0]
path=path[:-(len(path)-path.rindex("/")-1)]

#load ref seq records
refseqrecords = SeqIO.to_dict(SeqIO.parse(traindataset, "fasta"))

#load test seq records
testseqrecords = SeqIO.to_dict(SeqIO.parse(testdataset, "fasta"))

#Load classes, classification:
classificationdict,classes= LoadClassification(refseqrecords,classificationfilename)
testclassificationdict,testclasses= LoadClassification(testseqrecords,classificationfilename)
		
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
reportfilename=GetWorkingBase(prefix) + "." + basename + "_BLAST.classified"
classificationreportfilename=GetBase(reportfilename) + ".classification"
unclassifiedfastafilename=GetBase(reportfilename) + ".unclassified.fasta"
cutoffs={}
if cutoffsfilename!="" and cutoffsfilename!=None:
	with open(cutoffsfilename) as cutoffsfile:
		cutoffs = json.load(cutoffsfile)			
#save prediction
given_labels,assigned_labels,unclassifiedseqids,count=SavePrediction(classificationdict,refseqrecords.keys(),testclassificationdict,testseqrecords.keys(),cutoffs,cutoff,globalconfidence,bestscorelist,bestsimlist,bestcoveragelist,bestmatchlist,reportfilename,classificationreportfilename)
print("Number of classified sequences: " + str(count))
print("The results are saved in file  " + reportfilename + " and " + classificationreportfilename + ".")
if len(unclassifiedseqids) > 0:
	unclassifiedseqrecords=[]
	for seqid in unclassifiedseqids:
		unclassifiedseqrecords.append(testseqrecords[seqid])
	#write to fasta file	
	SeqIO.write(unclassifiedseqrecords, unclassifiedfastafilename, "fasta")	
	print("The unclassified sequences are saved in the file " +   unclassifiedfastafilename + ".")
#compute metrices
if len(given_labels) >0:
	reportname=GetBase(reportfilename) + ".report"
	reftaxa=LoadTaxa(classificationfilename)
	CalculateClassificationMetrics(given_labels,assigned_labels,reftaxa,reportname)
#generate krona report	
kronareport = GetBase(reportfilename) + ".krona.report"
kronahtml=GetBase(kronareport) + ".html"
classificationdict= LoadClassificationForKronaReport(reportfilename)
KronaPieCharts(classificationdict,kronareport,kronahtml)
print("The krona report and html are saved in files " + kronareport + " and " + kronahtml + ".") 	
	

