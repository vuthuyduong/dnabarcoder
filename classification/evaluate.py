#!/usr/bin/env python
# FILE: evaluate.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
#from sklearn.model_selection import StratifiedKFold
#from sklearn.preprocessing import StandardScaler
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from Bio import SeqIO

parser=argparse.ArgumentParser(prog='evaluate.py',  
							   usage="%(prog)s [options] -i predictionfile -qc queryclassification -rc referenceclassification",
							   description='''Script that computes metrices for the classification results.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the assigment/classified file')
parser.add_argument('-o','--out', help='The metrics output .')
parser.add_argument('-qc','--queryclassification', default="", help='the given classification file of the query sequences if exists, in tab. format to compute classification metrices.')
parser.add_argument('-rc','--refclassification', default="", help='the classification file os the reference sequences in the training dataset if exists, in tab. format to compute classification metrices for sequences having labels in the training dataset.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')

parser.add_argument('-label_columnname','--label_columnname',default="given label", help='the column name of the given labels in the classification file.')
parser.add_argument('-pred_columnname','--pred_columnname', default="prediction", help='the colunm name of the prediction.')
parser.add_argument('-proba_columnname','--proba_columnname', default="probability", help='the colunm name of the probability or confidence measure of the prediction.')
parser.add_argument('-minproba','--minproba', type=float, default=0, help='the colunm name of the probability or confidence measure of the prediction.')

parser.add_argument('-fullclassificationcolumnname','--fullclassificationcolumnname',default="full classification", help='the column name of the predicted full classifications in the classification file.')
parser.add_argument('-rankcolumnname','--rankcolumnname',default="rank", help='the column name of the ranks in the classification file.')
parser.add_argument('-rank','--rank',default="", help='the taxonomic of the prediction, in order to load the taxonomic name at the given rank in the case that the query labels are not given.')

args=parser.parse_args()
predictionfilename=args.input
queryclassificationfilename=args.queryclassification
refclassificationfilename=args.refclassification
outputfilename=args.out


def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]
	
# def GetWorkingBase(filename):
# 	basename=os.path.basename(filename)
# 	if "." in basename:
# 		basename=basename[:-(len(basename)-basename.rindex("."))] 
# 	path=outputpath + "/" + filename
# 	return path

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
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family + ";g__"+ genus + "s__unidentified"
	if level <1 and species!="unidentified":
		taxonname=species
		rank="species"
		classification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family + ";g__"+ genus+";s__"+species 	
	return classification,taxonname,rank

def LoadClassification(classificationfilename,idcolumnname):
	classificationdict={}
	if not os.path.exists(classificationfilename):
		return classificationdict
	classificationfile=open(classificationfilename, errors='ignore')
	header=next(classificationfile)
	i=0
	isError=False
	p_seqid=-1
	for text in header.split("\t"):
		if text.rstrip().lower()==idcolumnname.lower():
			p_seqid=i
		i=i+1	
	if 	p_seqid==-1:
		print("Please specify the sequence id columnname by using -idcolumnname.")	
		isError=True
	for line in classificationfile:
		texts=line.split("\t")
		seqid = texts[p_seqid].rstrip()
		classificationdict.setdefault(seqid,"")
		classification,taxonname,rank=GetTaxonomicClassification(0,header,texts)
		classificationdict[seqid]=classification
	classificationfile.close()	
	return classificationdict,isError

def LoadClassificationFromDescription(fastafilename):
	classificationdict={}
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
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
		classificationdict[seqid]=classification
	return classificationdict

def LoadTaxa(fastafilename,classificationfilename):
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	taxa={}
	classificationfile=open(classificationfilename)
	for line in classificationfile:
		line=line.rstrip()
		texts=line.split("\t")
		seqid=texts[0]
		if seqid in seqrecords.keys():
			for text in texts:
				if text!="":
					if not (text in taxa.keys()):
						taxa.setdefault(text,1)
					else:
						taxa[text]=taxa[text] + 1
	classificationfile.close()	
	return taxa

def LoadTaxaFromDescription(fastafilename):
	alltaxa={}
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	for seqid in seqrecords.keys():
		description=seqrecords[seqid].description
		if " " in description:
			description=description.split(" ")[1]
		texts=description.split("|")
		for text in texts:
			taxa=text.split(";")	
			for taxon in taxa:
				if "__" in taxon:
					taxon=taxon.split("__")[1]
					taxon=taxon.replace("_"," ")
				if taxon!="" and taxon!="unidentified":
					if  not (taxon in alltaxa.keys()):
						alltaxa.setdefault(taxon,1)
					else:	
						alltaxa[taxon]=alltaxa[taxon] + 1
	return alltaxa

def is_fasta(filename):
	try:
		with open(filename, "r") as handle:
			fasta = SeqIO.parse(handle, "fasta")
			return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file
	except ValueError:
		return False
	    
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

def LoadPrediction(predictionfilename,queryclassificationdict,outputname):
	if outputname!="":
		outputfile=open(outputname,"w")
	given_labels=[]
	pred_labels=[]
	probas=[]
	p_seqid=-1	
	p_rank=-1
	p_givenlabel=-1
	p_predlabel=-1
	p_proba=-1
	p_fullclassification=-1
	rank=args.rank
	predictionfile= open(predictionfilename)
	header=next(predictionfile)
	texts=header.split("\t")
	i=0
	for text in texts:
		text=text.rstrip()
		if args.label_columnname.lower() == text.lower():
			p_givenlabel=i
		if args.pred_columnname.lower() == text.lower():
			p_predlabel=i	
		if args.proba_columnname.lower() == text.lower():
			p_proba=i	
		if args.rankcolumnname.lower() == text.lower():
			p_rank=i
		if args.fullclassificationcolumnname.lower() == text.lower():
			p_fullclassification=i	
		if text.lower()==args.idcolumnname.lower():
			p_seqid=i	
		i=i+1	
# 	if p_seqid==-1 or p_givenlabel==-1 or p_predlabel==-1:
# 		print("Please enter columnnames given in the input file for sequence id using -idcolumnname, or given labels using -givenlabelcolumnname, and/or predicted labels using -predictedlabelcolumnname.")
# 		return [],[]
	if outputname!="":
		if p_givenlabel==-1:
			header=header.replace("\n","\t") + "Given label\n"
		outputfile.write(header)
	for line in predictionfile:
		texts=line.split("\t")
		seqid=texts[p_seqid]
		classname="unidentified"
		if p_predlabel >-1:
			classname=texts[p_predlabel].rstrip()
		if classname=="":
			classname="unidentified"
		givenlabel="unidentified"
		if p_givenlabel >-1:
			givenlabel=texts[p_givenlabel].rstrip()
		if givenlabel=="":
			givenlabel="unidentified"
		proba=0
		if p_proba >-1:
			proba=float(texts[p_proba].rstrip())
		if rank=="" and p_rank >-1:
			rank=texts[p_rank]
		level=GetLevel(rank)
		if seqid in queryclassificationdict.keys():
			queryclassification=queryclassificationdict[seqid]
			if level >=0:
				queryclassification=queryclassification.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","")
				givenlabel=queryclassification.split(";")[level]
		if givenlabel!="" and givenlabel!="unidentified":
			givenlabel=givenlabel.replace("_"," ")
			classname=classname.replace("_"," ")
			if outputname!="":
				texts[p_predlabel]=classname
				if p_givenlabel >=0:
					texts[p_givenlabel]=givenlabel
				newline=""
				for text in texts:
					newline=newline + text + "\t"
				if p_givenlabel==-1:
					newline=newline + givenlabel +"\t"
				newline=newline[:-1]
				if "\n" not in newline:
					newline=newline + "\n"
				outputfile.write(newline)
		if proba >= args.minproba:		
			given_labels.append(givenlabel)
			pred_labels.append(classname)
			probas.append(proba)
		classification=""
		if p_fullclassification > -1:
			classification=texts[p_fullclassification]
			classification=classification.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","")

	if outputname!="":
		outputfile.close()
	return given_labels,pred_labels,probas

def CalculateMetrics(test_labels,pred_labels,labels): 
	if len(test_labels)==0:
		return 0,0,0,0,[],[],[],0,[]
	accuracy=accuracy_score(test_labels,pred_labels)
	precision,recall,fscore,support=precision_recall_fscore_support(test_labels,pred_labels,average='macro')
	precisionvector,recallvector,fscorevector,support=precision_recall_fscore_support(test_labels,pred_labels,labels=labels)
	mcc=matthews_corrcoef(test_labels,pred_labels)
	#cohenkappascore=cohen_kappa_score(test_labels,pred_labels,labels=labels)
	#print("Cohenkappascore "+ str(cohenkappascore))
	confusionmatrix=confusion_matrix(test_labels,pred_labels,labels=labels)
	return round(accuracy,4),round(precision,4),round(recall,4),round(fscore,4),precisionvector,recallvector,fscorevector,round(mcc,4),confusionmatrix
	
def CalculateClassificationMetrics(givenlabels,predlabels,reftaxa,reportname,outputfilename):
	tmp_givenlabels=[]
	tmp_predlabels=[]
	i=0
	for label in givenlabels:
		predlabel=predlabels[i]
		if predlabel!="" and predlabel !="unidentified":
			tmp_givenlabels.append(label)		
			tmp_predlabels.append(predlabel)
# 		else:
# 			print(label + "\t" + predlabel)
		i=i+1		
	accuracy,precision,recall,fscore,precisionvector,recallvector,fscorevector,mcc,confusionmatrix=CalculateMetrics(givenlabels,predlabels,predlabels)
	tmpaccuracy,tmpprecision,tmprecall,tmpfscore,tmpprecisionvector,tmprecallvector,tmpfscorevector,tmpmcc,tmpconfusionmatrix=CalculateMetrics(tmp_givenlabels,tmp_predlabels,tmp_predlabels)
	filteredgivenlabels=[]
	filteredpredlabels=[]
	i=0
	for label in givenlabels:
		if label in reftaxa.keys():
			predlabel=predlabels[i]
			filteredgivenlabels.append(label)
			filteredpredlabels.append(predlabel)
		i=i+1			
	filteredaccuracy,filteredprecision,filteredrecall,filteredfscore,filteredprecisionvector,filteredrecallvector,filterdfscorevector,filteredmcc,filteredconfusionmatrix=CalculateMetrics(filteredgivenlabels,filteredpredlabels,filteredpredlabels)
	if os.path.exists(outputfilename):	
		outputfile=open(outputfilename,"a")
		outputfile.write(args.input + "\t" + str(len(givenlabels))+ "\t" + str(mcc) + "\t" + str(accuracy) + "\t" + str(recall) + "\t" + str(precision) + "\t" + str(fscore) + "\t")
		outputfile.write(str(len(tmp_givenlabels))+ "\t" + str(tmpmcc) + "\t" + str(tmpaccuracy) + "\t" + str(tmprecall) + "\t" + str(tmpprecision) + "\t" + str(tmpfscore) + "\t")
		outputfile.write(str(len(filteredgivenlabels))+ "\t" + str(filteredmcc) + "\t" + str(filteredaccuracy) + "\t" + str(filteredrecall) + "\t" + str(filteredprecision) + "\t" + str(filteredfscore) + "\n")
		outputfile.close()
	else:
		outputfile=open(outputfilename,"w")
		outputfile.write("Dataset\tNumber of sequences with a given label\tMcc\tAccuracy\tRecall\tPrecision\tFscore\tNumber of sequences with an identified predicted label\tMcc\tAccuracy\tRecall\tPrecision\tFscore\tNumber of sequences with a given label present in the training dataset\tMcc\tAccuracy\tRecall\tPrecision\tFscore\n")
		outputfile.write(args.input + "\t" + str(len(givenlabels))+ "\t" + str(mcc) + "\t" + str(accuracy) + "\t" + str(recall) + "\t" + str(precision) + "\t" + str(fscore) + "\t")
		outputfile.write(str(len(tmp_givenlabels))+ "\t" + str(tmpmcc) + "\t" + str(tmpaccuracy) + "\t" + str(tmprecall) + "\t" + str(tmpprecision) + "\t" + str(tmpfscore) + "\t")
		outputfile.write(str(len(filteredgivenlabels))+ "\t" + str(filteredmcc) + "\t" + str(filteredaccuracy) + "\t" + str(filteredrecall) + "\t" + str(filteredprecision) + "\t" + str(filteredfscore) + "\n")
		outputfile.close()
		
	report=open(reportname,"w")
	report.write("Taxonname\tAccuracy\tPrecision\tFscore\tNumber of reference sequences in the training dataset\n")
	i=0
	existinglabels=[]
	for label in pred_labels:
		if label in existinglabels:
			i=i+1
			continue
		existinglabels.append(label)
		count=0
		if label in reftaxa.keys():
			count=reftaxa[label]
		report.write(label + "\t" + str(recallvector[i]) + "\t" + str(precisionvector[i]) + "\t" + str(fscorevector[i]) + "\t" + str(count) + "\n")
		i=i+1
	report.close()
	

if __name__ == "__main__":
	queryclassificationdict={}
	if queryclassificationfilename!="":
		if is_fasta(queryclassificationfilename):
			queryclassificationdict=LoadClassificationFromDescription(queryclassificationfilename)
		else:	
			queryclassificationdict,isError=LoadClassification(queryclassificationfilename,args.idcolumnname)
			if isError==True:
				sys.exit()
	
	outputname=GetBase(predictionfilename) + ".labeled"
	reftaxa={}
	if refclassificationfilename!="":
		if is_fasta(refclassificationfilename):
			reftaxa=LoadTaxaFromDescription(refclassificationfilename)	
		else:	
			reftaxa=LoadTaxa(refclassificationfilename)	
	given_labels,pred_labels,probas= LoadPrediction(predictionfilename,queryclassificationdict,outputname)	
	if len(given_labels) >0:
		reportname=GetBase(outputname) + ".report"
		CalculateClassificationMetrics(given_labels,pred_labels,reftaxa,reportname,outputfilename)
		print("Prediction number: " + str(len(pred_labels)) + ".")
		print("Accuracy, precision, and fscore of the " + predictionfilename + " are given in file " + outputfilename +  ".")
		print("The assigned sequences with given labels are saved in file " + outputname + ".") 
		print("Accuracy, precision, and fscore of the taxa are given in file " + reportname +  ".")
	