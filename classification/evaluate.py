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

parser=argparse.ArgumentParser(prog='evaluate.py',  
							   usage="%(prog)s [options] -i predictionfile -c queryclassification -r referenceclassification",
							   description='''Script that visualizes using Krona and computes metrices for the classification results.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the assigment/classified file')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--queryclassification', required=True, help='the given classification file of the query if exists, in tab. format to compute classification metrices.')
parser.add_argument('-r','--refclassification', required=True, help='the given classification file of the query if exists, in tab. format to compute classification metrices.')
parser.add_argument('-seqidpos','--sequenceidposition', type=int,default=0, help='the position of sequence id in the classification file.')
parser.add_argument('-givenlabelpos','--givenlabelposition', type=int,default=-1, help='the position of given labels in the prediction file.')
parser.add_argument('-predpos','--predictionposition', type=int,default=-1, help='the position of predicted labels in the prediction file.')
parser.add_argument('-rankpos','--rankposition', type=int,default=-1, help='the position of ranks in the prediction file.')


args=parser.parse_args()
predictionfilename=args.input
queryclassificationfilename=args.queryclassification
refclassificationfilename=args.refclassification


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

def LoadClassification(classificationfilename):
	classificationdict={}
	if not os.path.exists(classificationfilename):
		return classificationdict
	classificationfile=open(classificationfilename)
	header=next(classificationfile)
	i=0
	p_seqid=0
	for text in header.split("\t"):
		text=text.rstrip().replace(" ","")
		if ("seqid" in text) or ("sequenceid" in text):
			p_seqid=i
		i=i+1	
	for line in classificationfile:
		texts=line.split("\t")
		seqid = texts[p_seqid].replace(">","").rstrip()
		classificationdict.setdefault(seqid,"")
		classification,taxonname,rank=GetTaxonomicClassification(0,header,texts)
		classificationdict[seqid]=classification
	classificationfile.close()	
	return classificationdict

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

def LoadPrediction(predictionfilename,queryclassificationdict,outputname,reftaxa):
	if outputname!="":
		outputfile=open(outputname,"w")
	given_labels=[]
	pred_labels=[]
	p_seqid=args.sequenceidposition
	p_rank=args.rankposition
	p_givenlabel=args.givenlabelposition
	p_predlabel=args.predictionposition
	predictionfile= open(predictionfilename, "r")
	header=next(predictionfile)
	if p_givenlabel==-1 or p_predlabel==-1 or p_rank==-1:
		texts=header.split("\t")
		i=0
		for text in texts:
			text=text.rstrip()
			if "givenlabel" in text.lower().replace(" ",""):
				p_givenlabel=i
			if "prediction" in text.lower().replace(" ",""):
				p_predlabel=i	
			if "rank" in text.lower().replace(" ",""):
				p_rank=i		
			i=i+1	
	if outputname!="":
		outputfile.write(header)
	for line in predictionfile:
		texts=line.split("\t")
		seqid=texts[p_seqid]
		classname=texts[p_predlabel]
		if classname=="":
			classname="unidentified"
		if classname=="unidentified":
			continue
		givenlabel=texts[p_givenlabel]
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
				texts[p_givenlabel]=givenlabel
				newline=""
				for text in texts:
					newline=newline + text + "\t"
				newline=newline[:-1]
				outputfile.write(newline)
			given_labels.append(givenlabel)
			pred_labels.append(classname)	
		classification=texts[3]
		classification=classification.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","")

	if outputname!="":
		outputfile.close()
	return given_labels,pred_labels

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
	for label in pred_labels:
		if label in existinglabels:
			i=i+1
			continue
		existinglabels.append(label)
		report.write(label + "\t" + str(recallvector[i]) + "\t" + str(precisionvector[i]) + "\t" + str(fscorevector[i]) + "\n")
		i=i+1
	report.close()
	print("Accuracy, precision, and fscore of the taxa are given in file " + reportname +  ".")

if __name__ == "__main__":
	queryclassificationdict=LoadClassification(queryclassificationfilename)
	outputname=GetBase(predictionfilename) + ".labeled"
	reftaxa=LoadTaxa(refclassificationfilename)	
	given_labels,pred_labels= LoadPrediction(predictionfilename,queryclassificationdict,outputname,reftaxa)
	if outputname!="":
		print("The assigned sequences with given labels are saved in file " + outputname + ".") 
	if len(given_labels) >0:
		reportname=GetBase(outputname) + ".report"
		CalculateClassificationMetrics(given_labels,pred_labels,reftaxa,reportname)
	