#!/usr/bin/env python
# FILE: assignment2funguildformat.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2021
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
from Bio import SeqIO

parser=argparse.ArgumentParser(prog='assignment2funguildformat.py',  
							   usage="%(prog)s [options] -i assignment/classificationfile",
							   description='''Script that creates a funguild format file from an assignment/classification file''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--assignmentfilename', required=True, help='the classification filename 1')
parser.add_argument('-o','--out',default="dnabarcoder", help='The output folder.')

args=parser.parse_args()
assignmentfilename= args.assignmentfilename
outputpath=args.out
if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)	

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + basename
	return path

def GetClassification(fullclassification):
	kingdom=""
	phylum=""
	bioclass=""
	order=""
	family=""
	genus=""
	species=""
	texts=fullclassification.split(";")
	for text in texts:
		if text.startswith("k__"):
			kingdom=text.replace("k__","")
		if text.startswith("p__"):
			phylum=text.replace("p__","")	
		if text.startswith("c__"):
			bioclass=text.replace("c__","")	
		if text.startswith("o__"):
			order=text.replace("o__","")	
		if text.startswith("f__"):
			family=text.replace("f__","")	
		if text.startswith("g__"):
			genus=text.replace("g__","")	
		if text.startswith("s__"):
			species=text.replace("s__","")	
	return kingdom,phylum,bioclass,order,family,genus,species

def Level2Number(level):
	number=0
	if level.lower()=="species":
		number=7
	elif level.lower()=="genus":
		number=6
	elif level.lower()=="family":
		number=5	
	elif level.lower()=="order":
		number=4
	elif level.lower()=="class":
		number=3
	elif level.lower()=="phylum":
		number=2
	elif level.lower()=="kingdom":
		number=1	
	return number
#######MAIN#######################
#outputclassificationclean= GetWorkingBase(assignmentfilename) + ".classification"
#outputfilename= GetWorkingBase(assignmentfilename) + ".funguild"
outputclassificationclean= assignmentfilename + ".clean"
outputfilename=assignmentfilename + ".funguild"
outputfile=open(outputfilename,"w")
outputclassificationcleanfile=open(outputclassificationclean,"w")
inputfile=open(assignmentfilename)
header=inputfile.readline()
words=header.split("\t")
p_s=-1
p_g=-1
p_f=-1
p_o=-1
p_c=-1
p_p=-1
p_k=-1
p_score=-1
p_level=-1
p_fullclassification=-1
i=0
for word in words:
	word=word.rstrip()
	if word.lower()=="kingdom":
		p_k=i
	if word.lower()=="phylum":
		p_p=i	
	if word.lower()=="class":
		p_c=i	
	if word.lower()=="order":
		p_o=i	
	if word.lower()=="family":
		p_f=i	
	if word.lower()=="genus":
		p_g=i	
	if word.lower()=="species":
		p_s=i
	if "score" in word.lower():
		p_score=i	
	if word.lower()=="level":
		p_level=i
	if word.lower()=="rank":
		p_level=i	
	if word.lower()=="full classification":
		p_fullclassification=i		
	i=i+1
outputfile.write("OTU_ID\ttaxonomy\n")
outputclassificationcleanfile.write("OTU_ID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tscore\n")
for line in inputfile:
	kingdom=""
	phylum=""
	bioclass=""
	order=""
	family=""
	genus=""
	species=""
	words=line.split("\t")
	seqid=words[0].split("|")[0]
	score=1
	number=8
	if p_score>-1:
		score=float(words[p_score])
		if score >1:
			score=score/100
	if p_level >-1:
		number=Level2Number(words[p_level])	
	fullclassification=""
	if p_fullclassification > -1:
		fullclassification=words[p_fullclassification]
		kingdom,phylum,bioclass,order,family,genus,species=GetClassification(fullclassification)
	else:
		if p_k< len(words) and p_k>-1:
			kingdom=words[p_k].rstrip()
		if kingdom=="":
			kingdom="Fungi"		   
		if p_p < len(words) and p_p>-1 and number >=2:
			phylum=words[p_p].rstrip()
		if phylum=="" or  "unculture" in phylum:
			phylum="unidentified"	
		if p_c<len(words) and p_c>-1 and number >=3:
			bioclass=words[p_c].rstrip()
		if bioclass=="" or  "unculture" in bioclass:	
			bioclass="unidentified"
		if p_o<len(words) and p_o>-1 and number >=4:
			order=words[p_o].rstrip()
		if order=="" or  "unculture" in order:
			order="unidentified"
		if p_f<len(words) and p_f>-1 and number >=5:
			family=words[p_f].rstrip()
		if family=="" or  "unculture" in family:
			family="unidentified"
		if p_g<len(words) and p_g > -1 and number >=6:
			genus=words[p_g].rstrip()
		if genus=="" or ("unculture" in genus):
			genus="unidentified"
		if p_s<len(words) and p_s >-1 and number >=7:
			species=(words[p_s]).rstrip()
		if species=="" or ("unculture" in species) or "_sp_" in species:
			species="unidentified"
		if " " in species:
			species=species.replace(" ","_")
		fullclassification= "k__" + kingdom + ";" + "p__" + phylum + ";" + "c__" + bioclass + ";" + "o__" + order + ";" + "f__" + family + ";" + "g__" + genus + ";" + "s__" + species + ";"
	outputclassificationcleanfile.write(seqid +"\t" + kingdom + "\t" + phylum + "\t" + bioclass +"\t" + order + "\t" + family + "\t" + genus + "\t" + species.replace("_"," ") + "\t" + str(score) + "\n") 
	outputfile.write(seqid + "\t" + str(round(score*100,2)) + "%|" + fullclassification + "\n")
inputfile.close()
outputfile.close()
outputclassificationcleanfile.close()
print("The outputs are saved in file: " + outputfilename + " and " + outputclassificationclean)
