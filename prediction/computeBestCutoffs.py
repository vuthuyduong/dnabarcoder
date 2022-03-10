#!/usr/bin/env python
# FILE: computeLocalCutoffs.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
from sklearn.metrics import precision_recall_fscore_support
from sklearn.metrics import cohen_kappa_score
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
import json
from Bio import SeqIO
import random
import multiprocessing
parser=argparse.ArgumentParser(prog='computeBestCutoffs.py',  
							   usage="%(prog)s [options] -i cutoffs -c classificationfile -o output",
							   description='''Script that computes best cutoffs of the taxa given in the cutoffs and classification file for sequence identification at different taxonomic levels.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input',required=True, help='the cutoffs file')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-cutoffs','--cutoffs', help='The json file containing the cutoffs to assign the sequences to the predicted taxa.')
#parser.add_argument('-minSeqNo','--minseqno', type=int, default=50, help='the minimum number of sequences for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
#arser.add_argument('-minGroupNo','--mingroupno', type=int, default=10, help='the minimum number of groups for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
parser.add_argument('-prefix','--prefix', help='the prefix of output filenames')
parser.add_argument('-savebestcutoffsascutoffs','--savebestcutoffsascutoffs', default="yes", help='the prefix of output filenames')

args=parser.parse_args()
classificationfilename=args.classification
cutoffsfilename=args.input
outputpath=args.out
prefix=args.prefix
#minGroupNo=args.mingroupno
#minSeqNo=args.minseqno

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

def GetRankTaxonomicClassification(level,classification):
	species=classification.split(";")[6].replace("s__","")
	genus=classification.split(";")[5].replace("g__","")
	family=classification.split(";")[4].replace("f__","")
	order=classification.split(";")[3].replace("o__","")
	bioclass=classification.split(";")[2].replace("c__","")
	phylum=classification.split(";")[1].replace("p__","")
	kingdom=classification.split(";")[0].replace("k__","")
	newclassification=""
	rank=""
	if level >=0 and kingdom!="unidentified":
		newclassification="k__" + kingdom 
		rank="kingdom"
	if level >=1 and phylum!="unidentified":
		newclassification="k__" + kingdom +";p__"+ phylum
		rank="phylum"
	if level >=2 and bioclass!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__"+bioclass
		rank="class"
	if level >=3 and order!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order 
		rank="order"
	if level >=4 and family!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family 
		rank="family"
	if level >=5 and genus!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family + ";g__"+ genus 
		rank="genus"
	if level >=6 and species!="unidentified":
		newclassification="k__" + kingdom +";p__"+phylum +";c__"+bioclass +";o__"+ order+";f__"+family + ";g__"+ genus+";s__"+species 		
		rank="species"
	return newclassification,rank

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
	if classificationfilename == "":
		return {}
	classificationfile= open(classificationfilename)
	header=next(classificationfile)
	for line in classificationfile:
		texts=line.split("\t")
		classification,taxonname,rank=GetTaxonomicClassification(0,header,texts)
		taxonnames=classification.split(";")
		level=0
		for taxonname in taxonnames:
			taxonname=taxonname.split("__")[1]
			if taxonname=="unidentified" or taxonname=="":
				level=level+1
				continue
			newclassification,rank=GetRankTaxonomicClassification(level,classification)	
			if taxonname in classificationdict.keys():
				if len(newclassification.replace("unidentified","")) > len(classificationdict[taxonname]):
					classificationdict[taxonname]["classification"]=newclassification
					classificationdict[taxonname]["rank"]=rank
			else:
				item={}
				item["classification"]=newclassification
				item["rank"]=rank
				classificationdict.setdefault(taxonname,item)		
			level=level+1	
	classificationfile.close()	
	return classificationdict

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

def GetHigherTaxa(rank,classification):
	highertaxa=[]
	if classification=="":
		return highertaxa
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
	highertaxa=GetHigherTaxa(rank,classification)
	highertaxa.append("All") 
	#seqno=0
	#groupno=0
	datasets=cutoffs[rank]
	maxconfidence=0
	bestcutoff=0
	bestminalignmentlength=0
	besttaxon=""
	for highertaxonname in highertaxa:
		if not highertaxonname in datasets.keys():
			continue
		cutoff=0
		if "cut-off" in datasets[highertaxonname].keys():
			cutoff=datasets[highertaxonname]["cut-off"]
		confidence=0
		if "confidence" in datasets[highertaxonname].keys():
			confidence=datasets[highertaxonname]["confidence"]
		minalignmentlength=0	
		if "min alignment length" in datasets[highertaxonname].keys():
			minalignmentlength=datasets[highertaxonname]["min alignment length"]		
#		if "sequence number" in datasets[highertaxonname].keys():
#			seqno=datasets[highertaxonname]["sequence number"]	
#		if "group number" in datasets[highertaxonname].keys():
#			groupno=datasets[highertaxonname]["group number"]	
		#if not ((seqno >0 and seqno < args.minseqno) or (groupno >0 and groupno < args.mingroupno)):	
		#if not (groupno < minGroupNo or seqno < minSeqNo):	
		if maxconfidence < confidence:
			maxconfidence =confidence
			bestcutoff=cutoff
			bestminalignmentlength=minalignmentlength
			besttaxon=highertaxonname	
	return bestcutoff,maxconfidence,besttaxon,bestminalignmentlength

def SaveBestCutoffsAsCutoffs(cutoffs,classificationdict,jsonoutputname,txtoutputname,problematicoutputname,problematicoutputname1):
	count=0
	total=0
	count1=0
	total1=0
	for rank in cutoffs.keys():
		datasets=cutoffs[rank]
		for datasetname in datasets.keys():
			dataset=datasets[datasetname]
			cutoff=dataset["cut-off"]
			confidence=dataset["confidence"]
			minalignmentlength=0
			if "min alignment length" in dataset.keys():
				minalignmentlength=dataset["min alignment length"]
			bestcutoff=cutoff
			maxconfidence=confidence
			bestminalignmentlength=minalignmentlength
			besttaxon=datasetname
			dataset["original cut-off"]=cutoff
			dataset["original confidence"]=confidence
			dataset["dataset with max confidence"]=datasetname
			dataset["original min alignment length"]=minalignmentlength
			if datasetname !="All":
				total=total +1	
			if datasetname in classificationdict.keys():
				classification=classificationdict[datasetname]["classification"]
				bestcutoff,maxconfidence,besttaxon,bestminalignmentlength=GetCutoffAndConfidence(rank,classification,cutoffs)					
			dataset["cut-off"]=bestcutoff
			dataset["confidence"]=maxconfidence
			dataset["dataset with max confidence"]=besttaxon
			dataset["min alignment length"]=bestminalignmentlength
	#save in json			
	with open(jsonoutputname,"w") as json_file:
		if sys.version_info[0] >= 3:
			json.dump(cutoffs,json_file,indent=2)
		else:
			json.dump(cutoffs,json_file,encoding='latin1',indent=2)
	#save as tab. format
	textfile=open(txtoutputname,"w")
	problematicfile=open(problematicoutputname,"w")
	problematicfile1=open(problematicoutputname1,"w")
	textfile.write("Rank\tDataset\toriginal cut-off\toriginal confidence\toriginal min alignment length\tDataset with max. confidence\tcut-off\tconfidence\tmin alignment length\tsequence number\tgroup number\tfasta filename\tclassification filename\n")
	problematicfile.write("Rank\tDataset\toriginal cut-off\toriginal confidence\toriginal min alignment length\tDataset with max. confidence\tcut-off\tconfidence\tmin alignment length\tsequence number\tgroup number\tfasta filename\tclassification filename\n")
	problematicfile1.write("Rank\tDataset\toriginal cut-off\toriginal confidence\toriginal min alignment length\tDataset with max. confidence\tcut-off\tconfidence\tmin alignment length\tsequence number\tgroup number\tfasta filename\tclassification filename\n")
	globaltotal=0
	globalcount=0
	for rank in cutoffs.keys():
		datasets=cutoffs[rank]
		#globalcutoff=0
		globalconfidence=0
		if "All" in datasets.keys():
			dataset=datasets["All"]
			#globalcutoff= dataset["cut-off"]
			globalconfidence=dataset["confidence"]
		for datasetname in datasets.keys():
			dataset=datasets[datasetname]
			cutoff=dataset["cut-off"]
			confidence=dataset["confidence"]
			minalignmentlength=dataset["min alignment length"]
			originalcutoff=dataset["original cut-off"]
			originalconfidence=dataset["original confidence"]
			originalminalignmentlength=dataset["original min alignment length"]
			highertaxon=dataset["dataset with max confidence"]
			SeqNo=dataset["sequence number"]
			GroupNo=dataset["group number"]
			fastafilename=dataset["fasta filename"]
			classificationfilename=dataset["classification filename"]
			#classificationposition=dataset["classification position"]
			textfile.write(rank+"\t" + datasetname + "\t"+str(originalcutoff)+"\t"+str(originalconfidence) + "\t" + str(originalminalignmentlength)+ "\t" + highertaxon +"\t" +str(cutoff)+"\t"+str(confidence)+ "\t" + str(minalignmentlength)+"\t"+str(SeqNo)+"\t"+str(GroupNo)+"\t"+fastafilename+"\t"+classificationfilename+"\n")
			if datasetname=="All":
				continue
			datasetname_rank=""
			if datasetname in classificationdict.keys():
				datasetname_rank=classificationdict[datasetname]["rank"]
			datasetname_level=GetLevel(datasetname_rank)
			level=GetLevel(rank)
			if datasetname_level==level - 1:
				total1=total1+1
			if confidence > originalconfidence and cutoff != originalcutoff:
				count=count+1
				problematicfile.write(rank+"\t" + datasetname + "\t"+str(originalcutoff)+"\t"+str(originalconfidence) + "\t" + str(originalminalignmentlength)+ "\t" + highertaxon +"\t" +str(cutoff)+"\t"+str(confidence)+ "\t" + str(minalignmentlength)+"\t"+str(SeqNo)+"\t"+str(GroupNo)+"\t"+fastafilename+"\t"+classificationfilename+"\n")
				if datasetname_level==level - 1:
					count1=count1+1
					problematicfile1.write(rank+"\t" + datasetname + "\t"+str(originalcutoff)+"\t"+str(originalconfidence) + "\t" + str(originalminalignmentlength)+ "\t" + highertaxon +"\t" +str(cutoff)+"\t"+str(confidence)+ "\t" + str(minalignmentlength)+"\t"+str(SeqNo)+"\t"+str(GroupNo)+"\t"+fastafilename+"\t"+classificationfilename+"\n")		
			if originalconfidence >= globalconfidence:
				globalcount=globalcount+1
			globaltotal=globaltotal+1				
	textfile.close()		
	problematicfile.close()
	problematicfile1.close()
	return count,total,count1,total1,globaltotal,globalcount

def SaveBestCutoffs(cutoffs,classificationdict,jsonoutputname,txtoutputname,problematicoutputname,problematicoutputname1):
	count=0
	total=0
	count1=0
	total1=0
	for rank in cutoffs.keys():
		datasets=cutoffs[rank]
		for datasetname in datasets.keys():
			dataset=datasets[datasetname]
			originalcutoff=dataset["cut-off"]
			originalconfidence=dataset["confidence"]
			originalminalignmentlength=0
			if "min alignment length" in dataset.keys():
				originalminalignmentlength=dataset["min alignment length"]
			bestcutoff=originalcutoff
			maxconfidence=originalconfidence
			bestminalignmentlength=originalminalignmentlength
			besttaxon=datasetname
			if datasetname !="All":
				total=total +1	
			if datasetname in classificationdict.keys():
				classification=classificationdict[datasetname]["classification"]
				bestcutoff,maxconfidence,besttaxon,bestminalignmentlength=GetCutoffAndConfidence(rank,classification,cutoffs)					
			dataset["best cut-off"]=bestcutoff
			dataset["max confidence"]=maxconfidence
			dataset["dataset with max confidence"]=besttaxon
			dataset["best min alignment length"]=bestminalignmentlength
	#save in json			
	with open(jsonoutputname,"w") as json_file:
		if sys.version_info[0] >= 3:
			json.dump(cutoffs,json_file,indent=2)
		else:
			json.dump(cutoffs,json_file,encoding='latin1',indent=2)
	#save as tab. format
	textfile=open(txtoutputname,"w")
	problematicfile=open(problematicoutputname,"w")
	problematicfile1=open(problematicoutputname1,"w")
	textfile.write("Rank\tDataset\tcut-off\tconfidence\tmin alignment length\tDataset with max. confidence\tbest cut-off\tmax confidence\tbest min alignment length\tsequence number\tgroup number\tfasta filename\tclassification filename\n")
	problematicfile.write("Rank\tDataset\tcut-off\tconfidence\tmin alignment length\tDataset with max. confidence\tbest cut-off\tmax confidence\tbest min alignment length\tsequence number\tgroup number\tfasta filename\tclassification filename\n")
	problematicfile1.write("Rank\tDataset\tcut-off\tconfidence\tmin alignment length\tDataset with max. confidence\tcut-off\tconfidence\tmin alignment length\tsequence number\tgroup number\tfasta filename\tclassification filename\n")
	globaltotal=0
	globalcount=0
	for rank in cutoffs.keys():
		datasets=cutoffs[rank]
		#globalcutoff=0
		globalconfidence=0
		if "All" in datasets.keys():
			dataset=datasets["All"]
			#globalcutoff= dataset["cut-off"]
			globalconfidence=dataset["confidence"]
		for datasetname in datasets.keys():
			dataset=datasets[datasetname]
			bestcutoff=dataset["best cut-off"]
			maxconfidence=dataset["max confidence"]
			bestminalignmentlength=dataset["best min alignment length"]
			originalcutoff=dataset["cut-off"]
			originalconfidence=dataset["confidence"]
			originalminalignmentlength=dataset["min alignment length"]
			highertaxon=dataset["dataset with max confidence"]
			SeqNo=dataset["sequence number"]
			GroupNo=dataset["group number"]
			fastafilename=dataset["fasta filename"]
			classificationfilename=dataset["classification filename"]
			#classificationposition=dataset["classification position"]
			textfile.write(rank+"\t" + datasetname + "\t"+str(originalcutoff)+"\t"+str(originalconfidence) + "\t" + str(originalminalignmentlength)+ "\t" + highertaxon +"\t" +str(bestcutoff)+"\t"+str(maxconfidence)+ "\t" + str(bestminalignmentlength)+"\t"+str(SeqNo)+"\t"+str(GroupNo)+"\t"+fastafilename+"\t"+classificationfilename+"\n")
			if datasetname=="All":
				continue
			datasetname_rank=""
			if datasetname in classificationdict.keys():
				datasetname_rank=classificationdict[datasetname]["rank"]
			datasetname_level=GetLevel(datasetname_rank)
			level=GetLevel(rank)
			if datasetname_level==level - 1:
				total1=total1+1
			if maxconfidence > originalconfidence and bestcutoff != originalcutoff:
				count=count+1
				problematicfile.write(rank+"\t" + datasetname + "\t"+str(originalcutoff)+"\t"+str(originalconfidence) + "\t" + str(originalminalignmentlength)+ "\t" + highertaxon +"\t" +str(bestcutoff)+"\t"+str(maxconfidence)+ "\t" + str(bestminalignmentlength)+"\t"+str(SeqNo)+"\t"+str(GroupNo)+"\t"+fastafilename+"\t"+classificationfilename+"\n")
				if datasetname_level==level - 1:
					count1=count1+1
					problematicfile1.write(rank+"\t" + datasetname + "\t"+str(originalcutoff)+"\t"+str(originalconfidence) + "\t" + str(originalminalignmentlength)+ "\t" + highertaxon +"\t" +str(bestcutoff)+"\t"+str(maxconfidence)+ "\t" + str(bestminalignmentlength)+"\t"+str(SeqNo)+"\t"+str(GroupNo)+"\t"+fastafilename+"\t"+classificationfilename+"\n")		
			if originalconfidence >= globalconfidence:
				globalcount=globalcount+1
			globaltotal=globaltotal+1				
	textfile.close()		
	problematicfile.close()
	problematicfile1.close()
	return count,total,count1,total1,globaltotal,globalcount
			
if __name__ == "__main__":
	if prefix=="" or prefix==None:
		prefix=GetBase(os.path.basename(cutoffsfilename))
	jsonoutputname=GetWorkingBase(prefix) + ".best.json"		
	txtoutputname=GetWorkingBase(prefix) + ".best.txt"	
	problematicoutputname=GetWorkingBase(prefix) + ".cutoffs.lowerconfidence.txt"
	problematicoutputname1=GetWorkingBase(prefix) + ".cutoffs.lowerconfidence.immediatehigherlevel.txt"	
	cutoffs={}
	if cutoffsfilename!="" and cutoffsfilename!=None:
		with open(cutoffsfilename) as cutoffsfile:
			cutoffs = json.load(cutoffsfile)
	classificationdict= LoadClassification(classificationfilename)
	if args.savebestcutoffsascutoffs=="yes":
		count,total,count1,total1,globaltotal,globalcount=SaveBestCutoffsAsCutoffs(cutoffs,classificationdict,jsonoutputname,txtoutputname,problematicoutputname,problematicoutputname1)
	else:
		count,total,count1,total1,globaltotal,globalcount=SaveBestCutoffs(cutoffs,classificationdict,jsonoutputname,txtoutputname,problematicoutputname,problematicoutputname1)
	if globaltotal >0:
		print("The best similarity cut-offs are saved in json and text format files " + jsonoutputname + " and " + txtoutputname + ".")
	print("The number of taxa having higher or equal prediction confidence than the global confidence: " + str(globalcount) + "/" + str(globaltotal) + " (" + str(round(globalcount*100/globaltotal,2)) +"%).")
	print("The number of taxa having the best cutoff: " + str(total-count) + "/" + str(total) + " (" + str(round((total-count)*100/total,2)) +"%).")
	print("The number of immedidate taxa having the best cutoff: " + str(total1-count1) + "/" + str(total1) + " (" + str(round((total1-count1)*100/total1,2)) +"%).")
	print("The local cut-offs for each taxon name at different taxonomic level are saved in file " + jsonoutputname + " and " + txtoutputname + ".")
	print("The taxa having lower prediction confidence and different cut-off than their higher taxa are saved in file " + problematicoutputname + ".")
	print("The immediate taxa having lower prediction confidence and different cut-off than their higher taxa are saved in file " + problematicoutputname1 + ".")
		
	