#!/usr/bin/env python
# FILE: assign.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
#from sklearn.model_selection import StratifiedKFold
#from sklearn.preprocessing import StandardScaler
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
parser=argparse.ArgumentParser(prog='assign.py',  
							   usage="%(prog)s [options] -i fastafile -p prediction -r referencefastafilename -c classificationfile -mp minproba -mc mincoverage -cutoffs cutoffsfile -o output",
							   description='''Script that assigns the classified sequences of the prediction file to their BLAST best match based on the given cutoffs.''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', help='the classified file')
parser.add_argument('-f','--fasta', help='the fasta file')
parser.add_argument('-r','--reference', required=True, help='the reference fasta file')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-c','--classification', required=True, help='the classification file in tab. format.')
parser.add_argument('-mp','--minproba', type=float, default=0, help='The minimum probability for verifying the classification results.')
parser.add_argument('-mc','--mincoverage', type=int, default=400, help='Minimum coverage required for the identitiy of the BLAST comparison.')
parser.add_argument('-m','--maxseqno', type=int, default=0, help='Maximum number of the sequences of the predicted taxon name from the classification file will be selected for the comparison to find the best match. If it is not given, all the sequences will be selected.')
parser.add_argument('-cutoff','--cutoff', type=float, default=0,help='The cutoff to assign the sequences to predicted taxa. If the cutoffs file is not given, this value will be taken for sequence assignment.')
parser.add_argument('-confidence','--confidence', type=float,default=0,help='The confidence of the cutoff to assign the sequences to predicted taxa')
parser.add_argument('-rank','--classificationrank', default="", help='the classification rank')
parser.add_argument('-prefix','--prefix', help='the prefix of output filenames')
parser.add_argument('-cutoffs','--cutoffs', help='The json file containing the cutoffs to assign the sequences to the predicted taxa.')
parser.add_argument('-minseqno','--minseqno', type=int, default=50, help='the minimum number of sequences for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')
parser.add_argument('-mingroupno','--mingroupno', type=int, default=10, help='the minimum number of groups for using the predicted cut-offs to assign sequences. Only needed when the cutoffs file is given.')

#parser.add_argument('-minconfidence','--minconfidence', type=float,default=0,help='The minimum confidence to assign the sequences to predicted taxa. If the cutoffs file is not given, this value will be taken for sequence assignment.')
#parser.add_argument('-minspeciescutoff','--minspeciescutoff', type=float, default=0, help='the minimum species cut-off to assign sequences. Only needed when the cutoffs file is given.')
#parser.add_argument('-mingenuscutoff','--mingenuscutoff', type=float, default=0, help='the minimum genus cut-off to assign sequences. Only needed when the cutoffs file is given.')
#parser.add_argument('-minfamilycutoff','--minfamilycutoff', type=float, default=0, help='the minimum family cut-off to assign sequences. Only needed when the cutoffs file is given.')
#parser.add_argument('-minordercutoff','--minordercutoff', type=float, default=0, help='the minimum order cut-off to assign sequences. Only needed when the cutoffs file is given.')
#parser.add_argument('-minclasscutoff','--minclasscutoff', type=float, default=0, help='the minimum class cut-off to assign sequences. Only needed when the cutoffs file is given.')
#parser.add_argument('-minphylumcutoff','--minphylumcutoff', type=float, default=0, help='the minimum phylum cut-off to assign sequences. Only needed when the cutoffs file is given.')
#parser.add_argument('-minkingdomcutoff','--minkingdomcutoff', type=float, default=0, help='the minimum kingdom cut-off to assign sequences. Only needed when the cutoffs file is given.')

args=parser.parse_args()
predictionfilename=args.input
fastafilename= args.fasta
referencefastafilename= args.reference
minprobaforBlast=args.minproba
mincoverage = args.mincoverage
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
			taxonnames=classification.split(";")
			for taxonname in taxonnames:
				taxonname=taxonname.split("__")[1]
				if taxonname=="unidentified" or taxonname=="":
					continue
				if taxonname in classes.keys():
					classes[taxonname].append(seqrecords[seqid])
				else:
					classes.setdefault(taxonname,[seqrecords[seqid]])		
	classificationfile.close()	
	return classificationdict,classes

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

#def GetAverageCutoffAndConfidence(rank,cutoffs):
#	if not rank in cutoffs.keys():
#		return [0,0]
#	datasets=cutoffs[rank]
#	averagecutoff=0
#	averageconfidence=0
#	n=0
#	for t in datasets.keys():
#		if t=="All":
#			continue
#		c=1
#		if "confidence" in datasets[t].keys():
#			c=datasets[t]["confidence"]
#		sn=1
#		if "sequence number" in datasets[t].keys():
#			sn=datasets[t]["sequence number"]	
#		if "cut-off" in datasets[t].keys():
#			averagecutoff=averagecutoff + datasets[t]["cut-off"]*sn
#			averageconfidence= averageconfidence + c*sn
#			n=n+sn
#	if n >0:
#		averagecutoff=float(averagecutoff/n)		
#		averageconfidence=float(averageconfidence/n)	
#	if (classificationrank==rank and classificationrank!="") or (classificationrank==""):	
#		print("At " + rank + " level: average cutoff: " + str(averagecutoff) + "; average confidence: " + str(averageconfidence) + ".")
#	return [averagecutoff,averageconfidence]
#
#def GetAverageCutoffs(cutoffs):
#	averagecutoffs={}
#	averagecutoffs["species"]=GetAverageCutoffAndConfidence("species",cutoffs)
#	averagecutoffs["genus"]=GetAverageCutoffAndConfidence("genus",cutoffs)
#	averagecutoffs["family"]=GetAverageCutoffAndConfidence("family",cutoffs)
#	averagecutoffs["order"]=GetAverageCutoffAndConfidence("order",cutoffs)
#	averagecutoffs["class"]=GetAverageCutoffAndConfidence("class",cutoffs)
#	averagecutoffs["phylum"]=GetAverageCutoffAndConfidence("phylum",cutoffs)
#	averagecutoffs["kingdom"]=GetAverageCutoffAndConfidence("kingdom",cutoffs)
#	averagecutoffs["minseqno"]=args.minseqno
#	averagecutoffs["mingroupno"]=args.mingroupno
#	averagecutoffs["minconfidence"]=args.confidence
#	return averagecutoffs

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
	#level=GetLevel(rank)							
	#p=level-1				
#	while p>=0:
#		highertaxonname=taxa[p]
#		if highertaxonname in datasets.keys():
#			if "cut-off" in datasets[highertaxonname].keys():
#				localcutoff=datasets[highertaxonname]["cut-off"]
#			if "confidence" in datasets[highertaxonname].keys():
#				confidence=datasets[highertaxonname]["confidence"]
#			if "sequence number" in datasets[highertaxonname].keys():
#				seqno=datasets[highertaxonname]["sequence number"]	
#			if "group number" in datasets[highertaxonname].keys():
#				groupno=datasets[highertaxonname]["group number"]	
##			if not (confidence < averagecutoffs["minconfidence"] or (seqno >0 and seqno < averagecutoffs["minseqno"]) or (groupno >0 and groupno < averagecutoffs["mingroupno"])):	
##				break
#		p=p-1		
#	if ((localcutoff==0) or confidence < averagecutoffs["minconfidence"] or (seqno >0 and seqno < averagecutoffs["minseqno"]) or (groupno >0 and groupno < averagecutoffs["mingroupno"])) and ("All" in datasets.keys()):
#		highertaxonname="All"
#		if "cut-off" in datasets[highertaxonname].keys():
#			localcutoff=datasets[highertaxonname]["cut-off"]
#		if "confidence" in datasets[highertaxonname].keys():
#			confidence=datasets[highertaxonname]["confidence"]
#		if "sequence number" in datasets[highertaxonname].keys():
#			seqno=datasets[highertaxonname]["sequence number"]	
#		if "group number" in datasets[highertaxonname].keys():
#			groupno=datasets[highertaxonname]["group number"]
#	#if even the prediction for All is not provided then the average cutoff will be used
#	if localcutoff==0:
#		localcutoff=averagecutoffs[rank][0]
#		confidence=averagecutoffs[rank][1]
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
#	taxacutoffs["species"]=GetCutoffAndConfidence("species",classification,cutoffs,averagecutoffs)
#	taxacutoffs["genus"]=GetCutoffAndConfidence("genus",classification,cutoffs,averagecutoffs)
#	taxacutoffs["family"]=GetCutoffAndConfidence("family",classification,cutoffs,averagecutoffs)
#	taxacutoffs["order"]=GetCutoffAndConfidence("order",classification,cutoffs,averagecutoffs)
#	taxacutoffs["class"]=GetCutoffAndConfidence("class",classification,cutoffs,averagecutoffs)
#	taxacutoffs["phylum"]=GetCutoffAndConfidence("phylum",classification,cutoffs,averagecutoffs)
#	taxacutoffs["kingdom"]=GetCutoffAndConfidence("kingdom",classification,cutoffs,averagecutoffs)
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
	return classification,taxonname,rank,level,localcutoff,confidence

#def GetMinCutoffs():
#	mincutoffs={}
#	mincutoffs["species"]={"cut-off":args.minspeciescutoff,"confidence":args.minconfidence,"sequence number":args.minseqno,"group number": args.mingroupno}
#	mincutoffs["genus"]={"cut-off":args.mingenuscutoff,"confidence":args.minconfidence,"sequence number":args.minseqno,"group number": args.mingroupno}
#	mincutoffs["family"]={"cut-off":args.minfamilycutoff,"confidence":args.minconfidence,"sequence number":args.minseqno,"group number": args.mingroupno}
#	mincutoffs["order"]={"cut-off":args.minordercutoff,"confidence":args.minconfidence,"sequence number":args.minseqno,"group number": args.mingroupno}
#	mincutoffs["class"]={"cut-off":args.minclasscutoff,"confidence":args.minconfidence,"sequence number":args.minseqno,"group number": args.mingroupno}
#	mincutoffs["phylum"]={"cut-off":args.minphylumcutoff,"confidence":args.minconfidence,"sequence number":args.minseqno,"group number": args.mingroupno}
#	mincutoffs["kingdom"]={"cut-off":args.minkingdomcutoff,"confidence":args.minconfidence,"sequence number":args.minseqno,"group number": args.mingroupno}
#	return mincutoffs

def Assign(classeswithsequences,classificationdict,testclassificationdict,minprobaforBlast,mincoverage,cutoffs,cutoff,globalconfidence,seqids,seqdict,labels,pred_labels,pred_classifications,probas,refids,bestscores,sims,coverages,maxseqno,outputname,classificationreportfilename):
	classificationlevel=GetLevel(classificationrank)
	output=open(outputname,"w")
	classificationreportfile=open(classificationreportfilename,"w")
	output.write("#SequenceID\tGiven label\tPrediction\tFull classification\tProbability\tRank\tCut-off\tConfidence\tReferenceID\tBLAST score\tBLAST sim\tBLAST coverage\n")
	classificationreportfile.write("SequenceID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\trank\tscore\tcutoff\tconfidence\n")
	i=0
	given_labels=[]
	assigned_labels=[]
	unclassifiedseqids=[]
	count=0
	averagecutoffs={}
#	if cutoffs!={}:
#		averagecutoffs=GetAverageCutoffs(cutoffs)
	for seqid in seqids:
		proba =1
		if i<len(probas):
			proba=probas[i]
		predictedname=pred_labels[i]
		rank=""
		giventaxonname=labels[i]
		classification=pred_classifications[i]
		bestscore=0
		sim=0
		coverage=0
		refid=""
		if  i<len(probas) and proba >= minprobaforBlast:
			if len(refids)==0:
				if seqid in seqdict.keys():
					reffilename,numberofsequences=CreateFastaFile(predictedname,classeswithsequences,maxseqno)
					refid,bestscore,sim,coverage=ComputeBestLocalBLASTScore(seqdict[seqid],reffilename,mincoverage)
					os.system("rm " + reffilename)	
			else:
				refid=refids[i]
				bestscore=bestscores[i]
				sim=sims[i]
				coverage=coverages[i]
			classification,predictedname,rank,level,cutoff,confidence=GetAssignment(refid,classificationdict,bestscore,cutoffs,cutoff,globalconfidence,classificationrank)			
			if sys.version_info[0] < 3:
				predictedname=unicode(predictedname,'latin1')
			if level>=classificationlevel and level !=-1:
				cleanclassification=classification.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","")
				classificationreportfile.write(seqid + "\t" + cleanclassification.replace(";","\t") + "\t" + rank + "\t" + str(bestscore) + "\t" + str(cutoff) + "\t" + str(confidence) + "\n")
			else:
				unclassifiedseqids.append(seqid)	
		if seqid in classificationdict.keys():
			giventaxa=classificationdict[seqid]['classification']
			giventaxa=giventaxa.replace("k__","").replace("p__","").replace("c__","")	.replace("o__","").replace("f__","").replace("g__","").replace("s__","")
			giventaxonname=giventaxa.split(";")[level]		
		if giventaxonname!="" and giventaxonname!="unidentified":	
			giventaxonname=giventaxonname.replace("_"," ")
			predictedname=predictedname.replace("_"," ")
			given_labels.append(giventaxonname)
			assigned_labels.append(predictedname)	
		if predictedname!="" and predictedname!="unidentified":			
			output.write(seqid + "\t" + giventaxonname + "\t"  + predictedname + "\t"+ classification + "\t" + str(proba) + "\t" + rank + "\t" + str(cutoff) + "\t" + str(confidence) + "\t" + refid + "\t" + str(bestscore) + "\t" + str(sim) + "\t" + str(coverage) + "\n")			
			count=count+1
		i=i+1
	output.close()
	classificationreportfile.close()
	return count,given_labels,assigned_labels,unclassifiedseqids
	
def LoadPrediction(predictionfilename):
	seqids=[]
	p_id=-1
	labels=[]
	p_l=-1
	pred_labels=[]
	p_pl=-1
	pred_classifications=[]
	p_c=-1
	probas=[]
	p_p=-1
	refids=[]
	p_r=-1
	bestscores=[]
	p_bs=-1
	sims=[]
	p_s=-1
	coverages=[]
	p_co=-1
	predictionfile=open(predictionfilename)
	headers=next(predictionfile)
	i=0
	for header in headers.rstrip().split("\t"):
		if "SequenceID" in header:
			p_id=i
		if "Given label" in header:
			p_l=i
		if "Prediction" in header:
			p_pl=i
		if "Full classification" in header:
			p_c=i
		if "Probability" in header:
			p_p=i
		if "ReferenceID" in header:
			p_r=i
		if "BLAST score" in header:
			p_bs=i	
		if "BLAST sim" in header:
			p_s=i
		if "BLAST coverage" in header:
			p_co=i	
		i=i+1	
	for line in predictionfile:
		texts=line.rstrip().split("\t")
		seqids.append(texts[p_id])
		labels.append(texts[p_l])
		pred_labels.append(texts[p_pl])
		pred_classifications.append(texts[p_c])
		proba=1
		if p_p >=0 and p_p <len(texts):
			proba=float(texts[p_p])
		probas.append(proba)
		if p_r>0:
			refids.append(texts[p_r])
			bestscores.append(float(texts[p_bs]))
			sims.append(float(texts[p_s]))
			coverages.append(float(texts[p_co]))
	return seqids,labels,pred_labels,pred_classifications,probas,refids,bestscores,sims,coverages
	
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

def CalculateMetrics(test_labels,pred_labels,labels): 
        accuracy=accuracy_score(test_labels,pred_labels)
        precision,recall,fscore,support=precision_recall_fscore_support(test_labels,pred_labels,average='macro')
        precisionvector,recallvector,fscorevector,support=precision_recall_fscore_support(test_labels,pred_labels,labels=labels)
        mcc=matthews_corrcoef(test_labels,pred_labels)
        #cohenkappascore=cohen_kappa_score(test_labels,pred_labels,labels=labels)
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
	outputname=GetWorkingBase(prefix) + ".assigned"
	if classificationrank!="":
		outputname=GetWorkingBase(prefix) + "." + classificationrank + ".assigned"
	if outputname==predictionfilename:
		outputname=outputname+".assigned"
	classificationreportfilename=GetBase(outputname) + ".classification"
	unclassifiedfastafilename=GetBase(outputname)  + ".unclassified.fasta"	
	#load prediction
	seqids,labels,pred_labels,pred_classifications,probas,refids,bestscores,sims,coverages=LoadPrediction(predictionfilename)	
	#load reference sequences
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	refseqrecords=SeqIO.to_dict(SeqIO.parse(referencefastafilename, "fasta"))
	classificationdict,classes= LoadClassification(refseqrecords,classificationfilename)
	testclassificationdict,testclasses= LoadClassification(seqrecords,classificationfilename)	
	cutoffs={}
	if cutoffsfilename!="" and cutoffsfilename!=None:
		with open(cutoffsfilename) as cutoffsfile:
			cutoffs = json.load(cutoffsfile)
	count,given_labels,assigned_labels,unclassifiedseqids=Assign(classes,classificationdict,testclassificationdict,minprobaforBlast,mincoverage,cutoffs,cutoff,globalconfidence,seqids,seqrecords,labels,pred_labels,pred_classifications,probas,refids,bestscores,sims,coverages,maxseqno,outputname,classificationreportfilename)
	print("Number of assigned sequences: " + str(count))
	print("The results are saved in file  " + outputname + " and " + classificationreportfilename + ".")
	unclassifiedseqrecords=[]
	if len(unclassifiedseqids) > 0:
		for seqid in unclassifiedseqids:
			unclassifiedseqrecords.append(seqrecords[seqid])
		#write to fasta file	
		SeqIO.write(unclassifiedseqrecords, unclassifiedfastafilename, "fasta")	
		print("The unclassified sequences are saved in the file " +   unclassifiedfastafilename + ".")
	
	if len(given_labels) >0:
		reportname=GetBase(outputname) + ".report"
		reftaxa=LoadTaxa(classificationfilename)
		CalculateClassificationMetrics(given_labels,pred_labels,reftaxa,reportname)
	#Compute classification metrices
	#making krona report
	kronareport = GetBase(outputname) + ".krona.report"
	kronahtml=GetBase(kronareport) + ".html"
	classificationdict= LoadClassificationForKronaReport(outputname)
	KronaPieCharts(classificationdict,kronareport,kronahtml)
	print("The krona report and html are saved in files " + kronareport + " and " + kronahtml + ".") 
