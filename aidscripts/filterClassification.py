#!/usr/bin/env python
# FILE: overview.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020

import sys
import os, argparse
from Bio import SeqIO

parser=argparse.ArgumentParser(prog='filterTaxonomicClassification.py',  
							   usage="%(prog)s [options] -i fastafile/subclassificationfile -c classificationfilename -o outputname",
							   description='''The script that filter the classification for the sequence ids given in the input file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta input file.')
parser.add_argument('-c','--classification', required=True, help='The taxonomic classification file.')
parser.add_argument('-seqidpos','--sequenceidposition', type=int, default=0, help='the position of the sequence ids in the classification file.')
parser.add_argument('-selectedpositions','--selectedpositions', default="", help='the selected positions for loading data from the classification file.')

parser.add_argument('-o','--out', default="", help='The output filename.')

args=parser.parse_args()
inputfilename= args.input
classificationfilename= args.classification
outputfilename=args.out
seqidpos= args.sequenceidposition
selectedpositions=args.selectedpositions


def GetBase(filename):
	if not ("." in filename):
		return filename
	return filename[:-(len(filename)-filename.rindex("."))]

#if not os.path.exists(outputpath):
#	os.system("mkdir " + outputpath)
	
#def GetWorkingBase(filename):
#	basename=os.path.basename(filename)
#	basename=basename[:-(len(basename)-basename.rindex("."))] 
#	path=outputpath + "/" + basename
#	return path

######MAIN
#load sequence ids
seqids =[]
inputfile=open(inputfilename)
for line in inputfile:
	if line.startswith(">"):
		seqid=line.rstrip().split(" ")[0]
		seqid=seqid.replace(">","")
		seqids.append(seqid)
	elif "\t" in line:
		seqid=line.split("\t")[seqidpos]
		seqids.append(seqid.rstrip())	
positionlist=[]
if selectedpositions !="":
	positions=selectedpositions.split(",")
	for position in positions:
		if position=="":
			continue
		positionlist.append(int(position))
		
if outputfilename==None or outputfilename=="":
	outputfilename=GetBase(inputfilename) + ".classification"		
taxafilename=GetBase(outputfilename) + ".taxonomic.classification"	
classificationfile=None
if sys.version_info[0] ==3 :
	classificationfile=open(classificationfilename, encoding='latin-1')
else:
	classificationfile=open(classificationfilename)
header=next(classificationfile)
if len(positionlist)>0:
	texts=header.split("\t")
	header=""
	for pos in positionlist:
		header=header + texts[pos].rstrip() + "\t" 
	header=header[:-1] + "\n"	
p_s=0
p_g=0
p_f=0
p_o=0
p_c=0
p_p=0
p_k=0
words=header.split("\t")
i=0
for word in words:
   if word.rstrip().lower()=="species":
	   p_s=i
   if word.rstrip().lower()=="genus":
	   p_g=i
   if word.rstrip().lower()=="family":
	   p_f=i	
   if word.rstrip().lower()=="order":
	   p_o=i
   if word.rstrip().lower()=="class":
	   p_c=i	
   if word.rstrip().lower()=="phylum":
	   p_p=i
   if word.rstrip().lower()=="kingdom":
	   p_k=i	        
   i=i+1
outputfile=open(outputfilename,"w")	            	   
taxafile=open(taxafilename,"w")
if (p_s+p_g+p_g+p_f+p_o+p_c+p_p+p_k)>0:
	taxafile.write("#Sequence ID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n")
outputfile.write(header)
for line in classificationfile:
	words=line.split("\t")	   
	seqid=words[seqidpos].rstrip()			
	if not seqid in seqids:
		continue
	words=line.split("\t")
	species=""
	if p_s>0 and p_s<len(words):
		species=words[p_s].rstrip()
	genus=""
	if p_g>0 and p_g<len(words):
		genus=words[p_g].rstrip()	
	family=""
	if p_f>0 and p_f<len(words):
		family=words[p_f].rstrip()	
	order=""
	if p_o>0 and p_o<len(words):
		order=words[p_o].rstrip()	
	bioclass=""
	if p_c>0 and p_c<len(words):
		bioclass=words[p_c].rstrip()	
	phylum=""
	if p_p>0 and p_p<len(words):
		phylum=words[p_p].rstrip()
	if p_k>0 and p_k<len(words):
		kingdom=words[p_k].rstrip()	
	if (p_s+p_g+p_g+p_f+p_o+p_c+p_p+p_k)>0:	
		taxclassification=seqid + "\t" + kingdom + "\t" + phylum + "\t" + bioclass + "\t" + order + "\t" + family + "\t" + genus + "\t" + species  + "\n"
		taxafile.write(taxclassification)
	classification=line	
	if len(positionlist)>0:
		texts=line.split("\t")
		classification=""
		for pos in positionlist:
			classification=classification + texts[pos].rstrip() + "\t" 
		classification=classification[:-1] + "\n"	
	outputfile.write(classification)
classificationfile.close()
taxafile.close()
outputfile.close()
if (p_s+p_g+p_g+p_f+p_o+p_c+p_p+p_k)>0:	
	print("The classifications are saved in  files " + outputfilename + " and " + taxafilename + ".")
else:
	os.system("rm " + taxafilename)
	print("The classifications are saved in  files " + outputfilename + ".")
	
