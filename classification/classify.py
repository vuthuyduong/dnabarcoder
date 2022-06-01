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
import random
import multiprocessing
parser=argparse.ArgumentParser(prog='classify.py',  
							   usage="%(prog)s [options] -i bestmatch/classified file -f the fasta file -r referencefastafilename -c classificationfile -mp minproba -mc mincoverage -cutoffs cutoffsfile -o output",
							   description='''Script that assigns the classified sequences of the prediction file to their BLAST best match based on the given cutoffs.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the classified file.')
parser.add_argument('-fmt','--inputformat', default="tab delimited", help='the format of the classified file. The inputfmt can have two values "tab delimited" and "blast". The value "tab delimited" is given as default, and the "blast" fmt is the format of the BLAST output with outfmt=6.')
parser.add_argument('-f','--fasta', default="", help='the fasta file')
parser.add_argument('-r','--reference', default="", help='the reference fasta file')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', default="", help='the classification file in tab. format.')
parser.add_argument('-mp','--minproba', type=float, default=0, help='The minimum probability for verifying the classification results.')
parser.add_argument('-ml','--minalignmentlength', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should probably be set to smaller, 50 for instance.')
parser.add_argument('-m','--maxseqno', type=int, default=0, help='Maximum number of the sequences of the predicted taxon name from the classification file will be selected for the comparison to find the best match. If it is not given, all the sequences will be selected.')
parser.add_argument('-cutoff','--cutoff', type=float, default=0,help='The cutoff to assign the sequences to predicted taxa. If the cutoffs file is not given, this value will be taken for sequence assignment.')
parser.add_argument('-confidence','--confidence', type=float,default=0,help='The confidence of the cutoff to assign the sequences to predicted taxa')
parser.add_argument('-rank','--classificationrank', default="", help='the classification rank')
parser.add_argument('-prefix','--prefix', help='the prefix of output filenames')
parser.add_argument('-cutoffs','--cutoffs', help='The json file containing the cutoffs to assign the sequences to the predicted taxa.')
parser.add_argument('-minseqno','--minseqno', type=int, default=0, help='the minimum number of sequences for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
parser.add_argument('-mingroupno','--mingroupno', type=int, default=0, help='the minimum number of groups for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
parser.add_argument('-save','--save',default="", help='The option to save all (default) or only classified sequences (-save classified) in the classification output.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')

args=parser.parse_args()
predictionfilename=args.input
fastafilename= args.fasta
referencefastafilename= args.reference
minprobaforBlast=args.minproba
mincoverage = args.minalignmentlength
cutoff=args.cutoff
classificationfilename=args.classification
classificationrank=args.classificationrank
maxseqno=args.maxseqno
prefix=args.prefix
outputpath=args.out
cutoffsfilename=args.cutoffs
globalconfidence=args.confidence

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

def LoadClassification(seqrecords,classificationfilename,idcolumnname):
	classificationdict={}
	classes={}
	classificationfile= open(classificationfilename)
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
			for taxonname in taxonnames:
				if taxonname=="":
					continue
				taxonname=taxonname.split("__")[1].replace("_"," ")
				if taxonname=="unidentified" or taxonname=="":
					continue
				if taxonname in classes.keys():
					classes[taxonname].append(seqrecords[seqid])
				else:
					classes.setdefault(taxonname,[seqrecords[seqid]])		
	classificationfile.close()	
	return classificationdict,classes,isError

def LoadClassificationFromDescription(seqrecords):
	classificationdict={}
	classes={}
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
		taxonname=""
		for t in taxonnames:
			t=t.split("__")[1]
			t=t.replace("_"," ")		
			if t=="unidentified" or t=="":
				level=level+1
				continue
			rank=ranks[level]
			taxonname=t
			if t in classes.keys():
				classes[t].append(seqrecords[seqid])
			else:
				classes.setdefault(t,[seqrecords[seqid]])
			level=level+1		
		if taxonname!="":
			classificationdict.setdefault(seqid,{})
			classificationdict[seqid]["classification"]=classification
			classificationdict[seqid]["taxonname"]=taxonname
			classificationdict[seqid]["rank"]=rank	
	return classificationdict,classes

#def LoadTaxa(classificationfilename):
#	taxa=[]
#	classificationfile=open(classificationfilename)
#	for line in classificationfile:
#		line=line.rstrip()
#		texts=line.split("\t")
#		for text in texts:
#			if text!="":
#				taxa.append(text)
#	classificationfile.close()			
#	return taxa

def CreateFastaFile(taxonname,classeswithsequences,maxseqno):
	if sys.version_info[0] < 3:
		taxonname=unicode(taxonname,errors='ignore')
	newfastafilename=""
	sequences=[]
	if taxonname in classeswithsequences.keys():
		sequences=classeswithsequences[taxonname]
		if len(sequences) >0:
			newfastafilename=taxonname.replace(" ","_").replace(".","_") + ".fasta"
			fastafile=open(newfastafilename,"w")
			if (maxseqno >0) and (len(sequences) > maxseqno):
				#select randomly 100 sequences to compare
				selectedlist=random.sample(range(0, len(sequences)), k=maxseqno)
				for i in selectedlist:
					sequence=sequences[i]
					seqid=sequence.id
					seq=str(sequence.seq)
					fastafile.write(">" + seqid + "\n")
					fastafile.write(seq + "\n")
			else:	
				for sequence in sequences:
					seqid=sequence.id
					seq=str(sequence.seq)
					fastafile.write(">" + seqid + "\n")
					fastafile.write(seq + "\n")
			fastafile.close()
	return newfastafilename,len(sequences)

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

def GetCutoffAndConfidence(rank,classification,cutoffs):
	if not rank in cutoffs.keys():
		return [0,0]
	#cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","")
	#taxa=cleanclassification.split(";")
	#taxa.append("All") 
	highertaxa=GetHigherTaxa(rank,classification)
	highertaxa.append("All") 
	localcutoff=0
	seqno=0
	groupno=0
	datasets=cutoffs[rank]
	maxconfidence=0
	bestcutoff=0
	for highertaxonname in highertaxa:
		if not highertaxonname in datasets.keys():
			continue
		if "cut-off" in datasets[highertaxonname].keys():
			localcutoff=datasets[highertaxonname]["cut-off"]
		else:
			continue
		confidence=1	
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
				break
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
	taxacutoffs["species"]=GetCutoffAndConfidence("species",classification,cutoffs)
	taxacutoffs["genus"]=GetCutoffAndConfidence("genus",classification,cutoffs)
	taxacutoffs["family"]=GetCutoffAndConfidence("family",classification,cutoffs)
	taxacutoffs["order"]=GetCutoffAndConfidence("order",classification,cutoffs)
	taxacutoffs["class"]=GetCutoffAndConfidence("class",classification,cutoffs)
	taxacutoffs["phylum"]=GetCutoffAndConfidence("phylum",classification,cutoffs)
	taxacutoffs["kingdom"]=GetCutoffAndConfidence("kingdom",classification,cutoffs)
	return taxacutoffs,kingdom,phylum,bioclass,order,family,genus,species

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
		elif bestscore >=taxacutoffs["genus"][0] and taxacutoffs["genus"][0]>0 and genus!="unidentified" and (classificationrank=="genus" or classificationrank== ""):
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
		elif bestscore >=taxacutoffs["phylum"][0] and taxacutoffs["phylum"][0]>0 and phylum!="unidentified" and (classificationrank=="phylum" or classificationrank== ""):
			rank="phylum"
			localcutoff=taxacutoffs["phylum"][0]
			confidence=taxacutoffs["phylum"][1]
			taxonname=phylum
			level=1
		elif bestscore >=taxacutoffs["kingdom"][0] and taxacutoffs["kingdom"][0]>0 and kingdom!="unidentified" and (classificationrank=="kingdom" or classificationrank== ""):
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
	else:
		classification=GetRankClassification(-1,classification)
	return classification,taxonname,rank,level,localcutoff,confidence

def Assign(classeswithsequences,refclassificationdict,queryclassificationdict,predictedclassificationdict,minprobaforBlast,mincoverage,cutoffs,cutoff,globalconfidence,seqdict,bestmatchdict,maxseqno,outputname,classificationreportfilename):
	#classificationlevel=GetLevel(classificationrank)
	output=open(outputname,"w")
	classificationreportfile=open(classificationreportfilename,"w")
	output.write("ID\tGiven label\tPrediction\tFull classification\tProbability\tRank\tCut-off\tConfidence\tReferenceID\tBLAST score\tBLAST sim\tBLAST coverage\n")
	classificationreportfile.write("ID\tReferenceID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\trank\tscore\tcutoff\tconfidence\n")
	given_labels=[]
	assigned_labels=[]
	unclassifiedseqids=[]
	count=0
	for seqid in bestmatchdict.keys():
		proba =1
		rank=""
		if "proba" in bestmatchdict[seqid].keys():
			proba=bestmatchdict[seqid]["proba"]
		predictedname= ""
		if "predlabel" in bestmatchdict[seqid].keys():
			predictedname=bestmatchdict[seqid]["predlabel"]
		giventaxonname=""
		if "label" in bestmatchdict[seqid].keys():
			giventaxonname=bestmatchdict[seqid]["label"]
		classification=""
		if "predclassification" in bestmatchdict[seqid].keys():
			classification=bestmatchdict[seqid]["predclassification"]
		refid=bestmatchdict[seqid]["refid"]
		bestscore=bestmatchdict[seqid]["score"]
		sim=bestmatchdict[seqid]["sim"]
		coverage=bestmatchdict[seqid]["alignmentlength"]
		if  proba >= minprobaforBlast:
			if refid=="":
				if seqid in seqdict.keys():
					reffilename,numberofsequences=CreateFastaFile(predictedname,classeswithsequences,maxseqno)
					if reffilename!="":
						newrefid,newbestscore,newsim,newcoverage=ComputeBestLocalBLASTScore(seqdict[seqid],reffilename,mincoverage)
						os.system("rm " + reffilename)	
						if newbestscore > bestscore:
							refid=newrefid
							bestscore=newbestscore
							sim=newsim
							coverage=newcoverage
			if refid!="":
				classification,predictedname,rank,level,cutoff,confidence=GetAssignment(refid,refclassificationdict,bestscore,cutoffs,cutoff,globalconfidence,classificationrank)			
			else:
				classification,predictedname,rank,level,cutoff,confidence=GetAssignment(seqid,predictedclassificationdict,bestscore,cutoffs,cutoff,globalconfidence,classificationrank)			
			if sys.version_info[0] < 3:
				predictedname=unicode(predictedname,'latin1')
		if seqid in refclassificationdict.keys():
			giventaxa=refclassificationdict[seqid]['classification']
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
		cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","").replace("_"," ")
		#save all classification"
		if args.save=="":
			#save all including unidentified sequences in the classification file
			output.write(seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t"+ classification + "\t" + str(proba) + "\t" + rank + "\t" + str(cutoff) + "\t" + str(confidence) + "\t" + refid + "\t" + str(bestscore) + "\t" + str(sim) + "\t" + str(coverage) + "\n")			
			classificationreportfile.write(seqid + "\t" + refid + "\t" + cleanclassification.replace(";","\t") + "\t" + rank + "\t" + str(bestscore) + "\t" + str(cutoff) + "\t" + str(confidence) + "\n")
		elif args.save=="classified":
			#save only the classified sequences in the classification file
			if predictedname!="" and predictedname!="unidentified":
				output.write(seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t"+ classification + "\t" + str(proba) + "\t" + rank + "\t" + str(cutoff) + "\t" + str(confidence) + "\t" + refid + "\t" + str(bestscore) + "\t" + str(sim) + "\t" + str(coverage) + "\n")			
				classificationreportfile.write(seqid + "\t" + refid + "\t" + cleanclassification.replace(";","\t") + "\t" + rank + "\t" + str(bestscore) + "\t" + str(cutoff) + "\t" + str(confidence) + "\n")
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
		if header.lower()==idcolumnname.lower():
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

def mergeDict(dict1,dict2):
	dict1.update(dict2)
	for key, value in dict1.items():
		if key in dict1 and key in dict2:
			dict1[key] = value + dict2[key]
	return dict1

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
		if cutoff >0:
			prefix =prefix + "." + str(cutoff)
	outputname=GetWorkingBase(prefix) + ".classified"
	if classificationrank!="":
		outputname=GetWorkingBase(prefix) + "." + classificationrank + ".classified"
	if outputname==predictionfilename:
		outputname=outputname+".classified"
	classificationreportfilename=GetBase(outputname) + ".classification"
	unclassifiedfastafilename=GetBase(outputname)  + ".unclassified.fasta"	
	#load prediction
	bestmatchdict={}
	if args.inputformat=="blast":
		bestmatchdict=LoadBlastOutput(predictionfilename,mincoverage)	
	else:
		bestmatchdict=LoadPrediction(predictionfilename,mincoverage,args.idcolumnname)	
	#load sequences
	seqrecords={}
	if os.path.exists(fastafilename):
		seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	#load reference sequences:	
	refseqrecords={}
	if os.path.exists(referencefastafilename):
		refseqrecords=SeqIO.to_dict(SeqIO.parse(referencefastafilename, "fasta"))
	refclassificationdict={}
	refclasses	={}
	queryclassificationdict={}
	queryclasses={}
	#load classification for the sequences
	if classificationfilename!="":
		refclassificationdict,refclasses,isError = LoadClassification(refseqrecords,classificationfilename,args.idcolumnname)
		if isError==True:
			sys.exit()
		queryclassificationdict,queryclasses,isError = LoadClassification(seqrecords,classificationfilename,args.idcolumnname)	
		if isError==True:
			sys.exit()
	else:
		refclassificationdict,refclasses = LoadClassificationFromDescription(refseqrecords)		
		queryclassificationdict,queryclasses = LoadClassificationFromDescription(seqrecords)	
	predictedclassificationdict,predictedclasses,isError= LoadClassification(seqrecords,predictionfilename,args.idcolumnname)	
	if isError==True:
		sys.exit()
	cutoffs={}
	if cutoffsfilename!="" and cutoffsfilename!=None:
		with open(cutoffsfilename) as cutoffsfile:
			cutoffs = json.load(cutoffsfile)		
	count,given_labels,assigned_labels,unclassifiedseqids=Assign(refclasses,refclassificationdict,queryclassificationdict,predictedclassificationdict,minprobaforBlast,mincoverage,cutoffs,cutoff,globalconfidence,seqrecords,bestmatchdict,maxseqno,outputname,classificationreportfilename)
	print("Number of classified sequences: " + str(count))
	#print("The results are saved in file  " + outputname)
	print("The results are saved in file  " + outputname + " and " + classificationreportfilename + ".")
	unclassifiedseqrecords=[]
	if len(unclassifiedseqids) > 0:
		for seqid in unclassifiedseqids:
			if seqid in seqrecords.keys():
				unclassifiedseqrecords.append(seqrecords[seqid])
		#write to fasta file	
		if len(unclassifiedseqrecords)>0:
			SeqIO.write(unclassifiedseqrecords, unclassifiedfastafilename, "fasta")	
			print("The unclassified sequences are saved in the file " +   unclassifiedfastafilename + ".")
	
#	if len(given_labels) >0:
#		reportname=GetBase(outputname) + ".report"
#		reftaxa=LoadTaxa(classificationfilename)
#		CalculateClassificationMetrics(given_labels,pred_labels,reftaxa,reportname)
	#Compute classification metrices
	#making krona report
	if count > 0:
		kronareport = GetBase(outputname) + ".krona.report"
		kronahtml=GetBase(kronareport) + ".html"
		classificationdict= LoadClassificationForKronaReport(outputname)
		KronaPieCharts(classificationdict,kronareport,kronahtml)
		print("The krona report and html are saved in files " + kronareport + " and " + kronahtml + ".") 
