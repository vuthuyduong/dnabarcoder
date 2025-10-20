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
#from sklearn.metrics import precision_recall_fscore_support
#from sklearn.metrics import cohen_kappa_score
#from sklearn.metrics import matthews_corrcoef
#from sklearn.metrics import confusion_matrix
#from sklearn.metrics import accuracy_score
import json
from Bio import SeqIO

def main():
	parser=argparse.ArgumentParser(prog='classify.py',  
								   usage="%(prog)s [options] -i bestmatch/classified file -r referencefastafilename -c classificationfile -ml minalignment -cutoffs cutoffsfile -o output",
								   description='''Script that assigns the classified sequences of the prediction file to their BLAST best match based on the given cutoffs.''',
								   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
	   )

	parser.add_argument('-i','--input', required=True, help='the classified file.')
	parser.add_argument('-f','--fasta', default="", help='The fasta file of the sequences for saving unidentified sequences. Optional.')
	parser.add_argument('-c','--classification', default="", help='the classification file in tab. format or the fasta file where classifications are given in the headers.')
	#parser.add_argument('-r','--reference', default="", help='the reference fasta file, in case the classification of the sequences is given in the sequence headers.')
	parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
	parser.add_argument('-fmt','--inputformat', default="tab delimited", help='the format of the classified file. The inputfmt can have two values "tab delimited" and "blast". The value "tab delimited" is given as default, and the "blast" fmt is the format of the BLAST output with outfmt=6.')
	parser.add_argument('-cutoff','--globalcutoff', type=float, default=-1,help='The global cutoff to assign the sequences to predicted taxa. If the cutoffs file is not given, this value will be taken for sequence assignment.')
	parser.add_argument('-confidence','--globalconfidence', type=float,default=-1,help='The global confidence to assign the sequences to predicted taxa')
	parser.add_argument('-rank','--classificationrank', default="", help='the classification rank')
	parser.add_argument('-prefix','--prefix', help='the prefix of output filenames')
	parser.add_argument('-cutoffs','--cutoffs', help='The json file containing the local cutoffs to assign the sequences to the predicted taxa.')
	parser.add_argument('-minseqno','--minseqno', type=int, default=0, help='the minimum number of sequences for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
	parser.add_argument('-mingroupno','--mingroupno', type=int, default=0, help='the minimum number of groups for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
	parser.add_argument('-ml','--minalignmentlength', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should probably be set to smaller, 50 for instance.')
	parser.add_argument('-saveclassifiedonly','--saveclassifiedonly',default=False, help='The option to save all (False) or only classified sequences (True) in the classification output.')
	parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')
	parser.add_argument('-display','--display',default="", help='If display=="yes" then the krona html is displayed.')

	args=parser.parse_args()
	predictionfilename=args.input
	globalcutoff=args.globalcutoff
	globalconfidence=args.globalconfidence
	cutoffsfilename=args.cutoffs
	classificationfilename=args.classification
	classificationrank=args.classificationrank
	fastafilename= args.fasta
	#referencefastafilename= args.reference
	mincoverage = args.minalignmentlength
	prefix=args.prefix
	outputpath=args.out
	classify(predictionfilename,args.idcolumnname,classificationfilename,globalcutoff,globalconfidence,cutoffsfilename,classificationrank,args.inputformat,mincoverage,args.minseqno,args.mingroupno,fastafilename,args.saveclassifiedonly,outputpath,prefix,args.display)

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]
	
def GetWorkingBase(filename,outputpath):
	basename=os.path.basename(filename)
	if "." in basename:
		basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + filename
	return path

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

def LoadClassification(classificationfilename,idcolumnname):
	classificationdict={}
	taxonomy={}
	classificationfile= open(classificationfilename,errors='ignore')
	header=next(classificationfile)
	seqidpos=-1
	isError=False
	texts=header.rstrip().split("\t")
	i=0
	for text in texts:
		if text.lower()==idcolumnname.lower():
			seqidpos=i
		i=i+1
	if 	seqidpos==-1:
		print("Please specify the sequence id columnname by using -idcolumnname or check the format of the classification file.")
		isError=True
		return classificationdict, taxonomy, isError
	for line in classificationfile:
		texts=line.split("\t")
		seqid = texts[seqidpos].rstrip()
		classification,taxonname,rank=GetTaxonomicClassification(0,header,texts)
		if classification!="":
			classificationdict.setdefault(seqid,{})
			classificationdict[seqid]["classification"]=classification
			classificationdict[seqid]["taxonname"]=taxonname
			classificationdict[seqid]["rank"]=rank
		taxonnames=classification.split(";")
		ranks=["kingdom","phylum","class","order","family","genus","species"]
		level=0
		for taxonname in taxonnames:				
			if taxonname=="":
				level=level+1
				continue
			taxonname=taxonname.split("__")[1].replace("_"," ")
			if taxonname=="unidentified" or taxonname=="":
				level=level+1
				continue
			rank=ranks[level]
			if not (taxonname in taxonomy.keys()):
				taxonomy.setdefault(taxonname,{})	
			rank=ranks[level]
			taxonomy[taxonname]["rank"]=rank
			taxonomy[taxonname]["classification"]=GetRankClassification(level, classification)	
			level=level+1	
	classificationfile.close()	
	return classificationdict,taxonomy,isError

def LoadClassificationFromDescription(classificationfilename):
	seqrecords={}
	if os.path.exists(classificationfilename):
		seqrecords=SeqIO.to_dict(SeqIO.parse(classificationfilename, "fasta"))
	classificationdict={}
	taxonomy={}
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
		rank=""
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
			if not (taxonname in taxonomy.keys()):
				taxonomy.setdefault(taxonname,{})
			taxonomy[taxonname]["rank"]=rank
			taxonomy[taxonname]["classification"]=GetRankClassification(level, classification)
			level=level+1		
		if currenttaxonname!="":
			classificationdict.setdefault(seqid,{})
			classificationdict[seqid]["classification"]=classification
			classificationdict[seqid]["taxonname"]=currenttaxonname
			classificationdict[seqid]["rank"]=rank	
	return classificationdict,taxonomy

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

def GetCutoffAndConfidence(rank,classification,cutoffs,minseqno,mingroupno):
	if not rank in cutoffs.keys():
		return [0,0,False]
	if classification=="k__Eukaryota;p__Streptophyta;c__Magnoliopsida;o__Fagales;f__Fagaceae;g__Castanea;s__unidentified": 
		print(classification)
		print(rank)
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
		if not ((seqno >0 and seqno < minseqno) or (groupno >0 and groupno < mingroupno)):	
			if maxconfidence < confidence:
				maxconfidence =confidence
				bestcutoff=localcutoff
			if isComputed==True:	
				break
	return [bestcutoff,maxconfidence,isComputed]

def GetCutoffs(classification,taxonomy):
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
			species=species.replace("_"," ")
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
	try:
		taxacutoffs["species"]=[taxonomy[species]["cut-off"],taxonomy[species]["confidence"],True]
	except KeyError:
		taxacutoffs["species"]=[0,0,False]	
	try:
		taxacutoffs["genus"]=[taxonomy[genus]["cut-off"],taxonomy[genus]["confidence"],True]
	except KeyError:
		taxacutoffs["genus"]=[0,0,False]		
	try:
		taxacutoffs["family"]=[taxonomy[family]["cut-off"],taxonomy[family]["confidence"],True]
	except KeyError:
		taxacutoffs["family"]=[0,0,False]
	try:
		taxacutoffs["order"]=[taxonomy[order]["cut-off"],taxonomy[order]["confidence"],True]
	except KeyError:
		taxacutoffs["order"]=[0,0,False]	
	try:
		taxacutoffs["class"]=[taxonomy[bioclass]["cut-off"],taxonomy[bioclass]["confidence"],True]
	except KeyError:
		taxacutoffs["class"]=[0,0,False]		
	try:	
		taxacutoffs["phylum"]=[taxonomy[phylum]["cut-off"],taxonomy[phylum]["confidence"],True]
	except KeyError:
		taxacutoffs["phylum"]=[0,0,False]	
	try:	
		taxacutoffs["kingdom"]=[taxonomy[kingdom]["cut-off"],taxonomy[kingdom]["confidence"],True]
	except KeyError:
		taxacutoffs["kingdom"]=[0,0,False]		
	return taxacutoffs,kingdom,phylum,bioclass,order,family,genus,species

def GetAssignment(refid,classificationdict,bestscore,taxonomy,classificationrank):
	localcutoff=taxonomy["unidentified"]["cut-off"]
	confidence=taxonomy["unidentified"]["confidence"]
	rank=classificationrank
	taxonname=""
	level=-1
	classification=""
	try:		
		refclassification=classificationdict[refid]['classification']
		taxacutoffs,kingdom,phylum,bioclass,order,family,genus,species=GetCutoffs(refclassification,taxonomy)
		if bestscore >=taxacutoffs["species"][0] and taxacutoffs["species"][2]==True and species!="unidentified" and (classificationrank=="species" or classificationrank== ""):
			rank="species"
			localcutoff=taxacutoffs["species"][0]
			confidence=taxacutoffs["species"][1]
			taxonname=species		
			level=6
		elif bestscore >=taxacutoffs["genus"][0] and taxacutoffs["genus"][2]==True and genus!="unidentified" and (classificationrank=="genus" or classificationrank== ""):
			rank="genus"
			localcutoff=taxacutoffs["genus"][0]
			confidence=taxacutoffs["genus"][1]
			taxonname=genus
			level=5
		elif bestscore >=taxacutoffs["family"][0] and taxacutoffs["family"][2]==True and family!="unidentified" and (classificationrank=="family" or classificationrank== ""):
			rank="family"
			localcutoff=taxacutoffs["family"][0]
			confidence=taxacutoffs["family"][1]
			taxonname=family	
			level=4
		elif bestscore >=taxacutoffs["order"][0] and taxacutoffs["order"][2]==True and order!="unidentified" and (classificationrank=="order" or classificationrank== ""):
			rank="order"
			localcutoff=taxacutoffs["order"][0]
			confidence=taxacutoffs["order"][1]
			taxonname=order
			level=3
		elif bestscore >=taxacutoffs["class"][0] and taxacutoffs["class"][2]==True and bioclass!="unidentified" and (classificationrank=="class" or classificationrank== ""):
			rank="class"
			localcutoff=taxacutoffs["class"][0]
			confidence=taxacutoffs["class"][1]
			taxonname=bioclass
			level=2
		elif bestscore >=taxacutoffs["phylum"][0] and taxacutoffs["phylum"][2]==True and phylum!="unidentified" and (classificationrank=="phylum" or classificationrank== ""):
			rank="phylum"
			localcutoff=taxacutoffs["phylum"][0]
			confidence=taxacutoffs["phylum"][1]
			taxonname=phylum
			level=1
		elif bestscore >=taxacutoffs["kingdom"][0] and taxacutoffs["kingdom"][2]==True and kingdom!="unidentified" and (classificationrank=="kingdom" or classificationrank== ""):
			rank="kingdom"
			localcutoff=taxacutoffs["kingdom"][0]
			confidence=taxacutoffs["kingdom"][1]
			taxonname=kingdom
			level=0	
		level=GetLevel(rank)
		classification=GetRankClassification(level,refclassification)
	except KeyError:
		classification=GetRankClassification(-1,classification)
	return classification,taxonname,rank,level,localcutoff,confidence

def Assign(refclassificationdict,taxonomy,bestmatchdict,outputname,classificationreportfilename,classificationrank,saveclassifiedonly):
	#classificationlevel=GetLevel(classificationrank)
	output=open(outputname,"w")
	classificationreportfile=open(classificationreportfilename,"w")
	output.write("ID\tGiven label\tPrediction\tFull classification\tRank\tCut-off\tConfidence\tReferenceID\tBLAST score\tBLAST sim\tBLAST coverage\n")
	classificationreportfile.write("ID\tReferenceID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\trank\tscore\tcutoff\tconfidence\n")
	given_labels=[]
	assigned_labels=[]
	unclassifiedseqids=[]
	count=0
	for seqid in bestmatchdict.keys():
		rank=""
		level=-1
		predictedname= ""
		if "predlabel" in bestmatchdict[seqid].keys():
			predictedname=bestmatchdict[seqid]["predlabel"]
		giventaxonname=""
		if "label" in bestmatchdict[seqid].keys():
			giventaxonname=bestmatchdict[seqid]["label"]
		classification=""	
		if "predclassification" in bestmatchdict[seqid].keys():
			classification=bestmatchdict[seqid]["predclassification"]
		if classification=="":
			classification="k__unidentified;p__unidentified;c__unidentified;o__unidentified;f__unidentified;g__unidentified;s__unidentified"		
		refid=bestmatchdict[seqid]["refid"]
		bestscore=bestmatchdict[seqid]["score"]
		sim=bestmatchdict[seqid]["sim"]
		coverage=bestmatchdict[seqid]["alignmentlength"]
		confidence=-1
		cutoff=-1
		if refid!="":
			classification,predictedname,rank,level,cutoff,confidence=GetAssignment(refid,refclassificationdict,bestscore,taxonomy,classificationrank)		
		cutoff_str=str(cutoff)	
		if cutoff==-1:
			cutoff_str="N/A"	
		confidence_str=str(confidence)	
		if confidence==-1:
			confidence_str="N/A"
		if sys.version_info[0] < 3:
			predictedname=unicode(predictedname,'latin1')
		try:
			giventaxadict=refclassificationdict[seqid]
		except KeyError:
			pass
		else:
			giventaxa=giventaxadict['classification']
			giventaxa=giventaxa.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","")
			giventaxonname=giventaxa.split(";")[level]		
		if giventaxonname!="" and giventaxonname!="unidentified":	
			giventaxonname=giventaxonname.replace("_"," ")
			predictedname=predictedname.replace("_"," ")
			given_labels.append(giventaxonname)
			assigned_labels.append(predictedname)	
		if predictedname!="" and predictedname!="unidentified":
			count=count+1
		else:
			unclassifiedseqids.append(seqid)	
			rank=""
		cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","").replace("o__","").replace("f__","").replace("g__","").replace("s__","").replace("_"," ")
		#save all classification"
		if saveclassifiedonly==False:
			#save all including unidentified sequences in the classification file
			output.write(seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t"+ classification + "\t" + rank + "\t" + cutoff_str + "\t" + confidence_str + "\t" + refid + "\t" + str(bestscore) + "\t" + str(sim) + "\t" + str(coverage) + "\n")			
			classificationreportfile.write(seqid + "\t" + refid + "\t" + cleanclassification.replace(";","\t") + "\t" + rank + "\t" + str(bestscore) + "\t" + cutoff_str + "\t" + confidence_str + "\n")
		else:
			#save only the classified sequences in the classification file
			if predictedname!="" and predictedname!="unidentified":
				output.write(seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t"+ classification + "\t" + rank + "\t" + cutoff_str + "\t" + confidence_str + "\t" + refid + "\t" + str(bestscore) + "\t" + str(sim) + "\t" + str(coverage) + "\n")			
				classificationreportfile.write(seqid + "\t" + refid + "\t" + cleanclassification.replace(";","\t") + "\t" + rank + "\t" + str(bestscore) + "\t" + cutoff_str + "\t" + confidence_str + "\n")
	output.close()
	classificationreportfile.close()
	return count,given_labels,assigned_labels,unclassifiedseqids
	
def LoadPrediction(predictionfilename,mincoverage,idcolumnname):
	bestmatchdict={}
	p_id=-1
	p_l=-1
	p_pl=-1
	p_c=-1
	p_p=-1
	p_r=-1
	p_bs=-1
	p_s=-1
	p_co=-1
	predictionfile=open(predictionfilename)
	headers=next(predictionfile)
	i=0
	for header in headers.rstrip().split("\t"):
		if header.lower()==idcolumnname.lower() or header.lower()=="id":
			p_id=i
		if "given label" in header.lower():
			p_l=i
		if "prediction" in header.lower():
			p_pl=i
		if "full classification" in header.lower():
			p_c=i
		if "probability" in header.lower():
			p_p=i
		if ("referenceid" in header.lower()) or ("reference_id" in header.lower()) or ("reference id" in header.lower()) or ("refid" in header.lower()) :
			p_r=i
		if "score" in header.lower():
			p_bs=i	
		if "sim" in header.lower():
			p_s=i
		if "coverage" in header.lower():
			p_co=i	
		i=i+1	
	for line in predictionfile:
		texts=line.rstrip().split("\t")
		seqid=""
		if p_id >=0 and p_id < len(texts):
			seqid=texts[p_id]
		if seqid=="":
			continue
		label=""
		if p_l >=0 and p_l < len(texts):
			label=texts[p_l]
		predlabel=""
		if p_pl >=0 and p_pl < len(texts):
			predlabel=texts[p_l]
		predclassification=""
		if p_c >=0 and p_c < len(texts):
			predclassification=texts[p_c]
		proba=1
		if p_p >=0 and p_p <len(texts):
			proba=float(texts[p_p])
		refid=""
		if p_r>0 and p_r<len(texts):
			refid=texts[p_r]
		sim=0
		if p_s>0 and p_s<len(texts):
			sim=float(texts[p_s])
		alignmentlength=""
		if p_co>0 and p_co<len(texts):
			alignmentlength=float(texts[p_co])
		score=0
		if p_bs>0 and p_bs<len(texts):
			score=float(texts[p_bs])
		else:
			score=sim
			if alignmentlength < mincoverage:
				score=float(score * alignmentlength)/mincoverage
		bestmatchdict.setdefault(seqid,{})
		bestmatchdict[seqid]["refid"]=refid
		bestmatchdict[seqid]["label"]=label
		bestmatchdict[seqid]["predlabel"]=predlabel
		bestmatchdict[seqid]["predclassification"]=predclassification
		bestmatchdict[seqid]["proba"]=proba
		bestmatchdict[seqid]["score"]=score
		bestmatchdict[seqid]["sim"]=sim
		bestmatchdict[seqid]["alignmentlength"]=alignmentlength
	return bestmatchdict

def LoadBlastOutput(blastoutput,mincoverage):
	bestmatchdict={}
	#read blast output
	blastoutputfile = open(blastoutput)
	refid = ""
	score=0
	seqid=""
	for line in blastoutputfile:
		words = line.split("\t")
		seqid=words[0]
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		refid=words[1]
		score=sim
		if coverage < mincoverage:
			score=float(score * coverage)/mincoverage
		if not (seqid in bestmatchdict.keys()):
			bestmatchdict.setdefault(seqid,{})
			bestmatchdict[seqid]["refid"]=""
			bestmatchdict[seqid]["score"]=0
			bestmatchdict[seqid]["sim"]=0
			bestmatchdict[seqid]["alignmentlength"]=0
		#if (score > bestmatchdict[seqid]["score"]):
		if (score > bestmatchdict[seqid]["score"]) or (score == bestmatchdict[seqid]["score"] and coverage > bestmatchdict[seqid]["alignmentlength"]):
			bestmatchdict[seqid]["refid"]=refid
			bestmatchdict[seqid]["score"]=score
			bestmatchdict[seqid]["sim"]=sim
			bestmatchdict[seqid]["alignmentlength"]=coverage
	return bestmatchdict
	
def GetClassificationpos(pred_labels,classificationfilename):
	classificationpos=0
	rank=""
	classificationfile=open(classificationfilename)
	for line in classificationfile:
		if line.startswith("#"):
			continue
		texts=line.rstrip().split("\t")	 
		i=0
		for text in texts:
			if text!="":
				if text in pred_labels:
					classificationpos=i
					break
			i=i+1	
		if classificationpos >0:
			break
	classificationfile.close()	
	return classificationpos,rank

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

def AddCutoffsToTaxonomy(taxonomy,globalcutoff,globalconfidence,cutoffs,minseqno,mingroupno):
	for taxonname in taxonomy.keys():
		if cutoffs!={}:
			classification=taxonomy[taxonname]["classification"]
			rank=taxonomy[taxonname]["rank"]
			if taxonname in cutoffs.keys():
				taxonomy[taxonname]["cut-off"]=cutoffs[taxonname]["cut-off"]
				taxonomy[taxonname]["confidence"]=cutoffs[taxonname]["confidence"]
			else:
				cutoff_confidence=GetCutoffAndConfidence(rank,classification,cutoffs,minseqno,mingroupno)
				if cutoff_confidence[2]==True:
					taxonomy[taxonname]["cut-off"]=cutoff_confidence[0]
					taxonomy[taxonname]["confidence"]=cutoff_confidence[1]
			if not ("cut-off" in taxonomy[taxonname].keys()) and globalcutoff >=0: #use the globalcutoff
				taxonomy[taxonname]["cut-off"]=globalcutoff
				taxonomy[taxonname]["confidence"]=globalconfidence
		else:
			taxonomy[taxonname]["cut-off"]=globalcutoff
			taxonomy[taxonname]["confidence"]=globalconfidence
	taxonomy.setdefault("unidentified",{})
	taxonomy["unidentified"]["cut-off"]=globalcutoff
	taxonomy["unidentified"]["confidence"]=globalconfidence		

def KronaPieCharts(classification,kronareport,kronahtml,display):
	kronareportfile=open(kronareport,"w")
	for classname in classification.keys():
		kronareportfile.write(str(classification[classname]) + "\t" + classname + "\n")
	kronareportfile.close()	
	#create kronahtml
	if os.path.exists("ImportText.pl"):
		command="ImportText.pl " + kronareport + " -o " + kronahtml
		#print(command)
		os.system(command)
	else:
		command="ktImportText " + kronareport + " -o " + kronahtml
		#print(command)
		os.system(command)
	if display=="yes":
		os.system("firefox " + kronahtml) 
		
def is_fasta(filename):
	try:
		with open(filename, "r") as handle:
			fasta = SeqIO.parse(handle, "fasta")
			return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file
	except ValueError:
		return False		
		
def classify(predictionfilename,idcolumnname,classificationfilename,globalcutoff,globalconfidence,cutoffsfilename,classificationrank,inputformat,mincoverage,minseqno,mingroupno,fastafilename,saveclassifiedonly,outputpath,prefix,display):
	if not os.path.exists(outputpath):
		os.system("mkdir " + outputpath)
	if prefix=="" or prefix==None:
		prefix=GetBase(predictionfilename)
		if "/" in prefix:
			prefix=prefix[prefix.rindex("/")+1:]	
		if globalcutoff >0 and cutoffsfilename=="":
			prefix =prefix + "." + str(globalcutoff)
	outputname=GetWorkingBase(prefix,outputpath) + ".classified"
	if classificationrank!="":
		outputname=GetWorkingBase(prefix,outputpath) + "." + classificationrank + ".classified"
	if outputname==predictionfilename:
		outputname=outputname+".classified"
	classificationreportfilename=GetBase(outputname) + ".classification"
	unclassifiedfastafilename=GetBase(outputname)  + ".unclassified.fasta"	
	#load prediction
	bestmatchdict={}
	if inputformat=="blast":
		bestmatchdict=LoadBlastOutput(predictionfilename,mincoverage)	
	else:
		bestmatchdict=LoadPrediction(predictionfilename,mincoverage,idcolumnname)	
	refclassificationdict={}
	#load classification for the sequences
	if not is_fasta(classificationfilename):
		refclassificationdict,taxonomy,isError = LoadClassification(classificationfilename,idcolumnname)
		if isError==True:
			sys.exit()
	else:
		#load reference sequences, in case the classification of the sequences is given in sequence headers
		
		refclassificationdict,taxonomy = LoadClassificationFromDescription(classificationfilename)

	cutoffs={}
	if cutoffsfilename!="" and cutoffsfilename!=None:
		with open(cutoffsfilename) as cutoffsfile:
			cutoffs = json.load(cutoffsfile)
	#add cutoffs to taxa for sequence identification
	AddCutoffsToTaxonomy(taxonomy,globalcutoff,globalconfidence,cutoffs,minseqno,mingroupno)
	count,given_labels,assigned_labels,unclassifiedseqids=Assign(refclassificationdict,taxonomy,bestmatchdict,outputname,classificationreportfilename,classificationrank,saveclassifiedonly)
	print("Number of classified sequences: " + str(count))
	#print("The results are saved in file  " + outputname)
	print("The results are saved in file  " + outputname + " and " + classificationreportfilename + ".")
	unclassifiedseqrecords=[]
	if len(unclassifiedseqids) > 0:
		#load sequences if the fasta file of the sequences is given, to save unidentified sequences
		seqrecords={}
		if os.path.exists(fastafilename):
			seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))	
		for seqid in unclassifiedseqids:
			if seqid in seqrecords.keys():
				unclassifiedseqrecords.append(seqrecords[seqid])
		#write to fasta file	
		if len(unclassifiedseqrecords)>0:
			SeqIO.write(unclassifiedseqrecords, unclassifiedfastafilename, "fasta")	
			print("The unclassified sequences are saved in the file " +   unclassifiedfastafilename + ".")
	#making krona report
	if count > 0:
		kronareport = GetBase(outputname) + ".krona.report"
		kronahtml=GetBase(kronareport) + ".html"
		classificationdict= LoadClassificationForKronaReport(outputname)
		KronaPieCharts(classificationdict,kronareport,kronahtml,display)
		print("The krona report and html are saved in files " + kronareport + " and " + kronahtml + ".") 
		
	
if __name__ == "__main__":
	main()
