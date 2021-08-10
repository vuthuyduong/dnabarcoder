#!/usr/bin/env python
# FILE: overview.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020

import sys
import os, argparse
from Bio import SeqIO

parser=argparse.ArgumentParser(prog='overview.py',  
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -out outputname",
							   description='''The script that summarizes the number of taxonomic groups of the sequences of the given fasta file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', default="", help='the fasta input file.')
parser.add_argument('-c','--classification', required=True, help='The taxonomic classification file.')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename= args.classification
outputpath=args.out
if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)

def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + basename
	return path

def ReportAtLevel(level,seqids):
	if level <0:
		return 0,0,0
	classificationfile = open(classificationfilename)
	next(classificationfile)
	taxa=[]
	seqNumbers=[]
	seqNo=0
	count=0
	for line in classificationfile:
		count=count+1
		words=line.split("\t")
		seqid=words[0].replace(".1","").replace(">","")
		if len(seqids) > 0 and (seqid not in seqids):
			continue
		if level >=  len(words):
			continue
		name=line.split("\t")[level].rstrip()
		if name=="" or name=="unidentified" or ("uncultured" in name) or ("_sp_" in name):
			continue
		if name in taxa:
			i=taxa.index(name)
			seqNumbers[i]=seqNumbers[i]+1
			seqNo=seqNo+1
		elif name !="":
			taxa.append(name)
			seqNumbers.append(1)	
			seqNo=seqNo+1
	numberoftaxa=len(taxa)		
	classificationfile.close()
	return numberoftaxa,seqNo,count

######MAIN
outputfilename=""
seqids=[]
if fastafilename != "":
	seqrecords = SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	seqids=seqrecords.keys()
	outputfilename=GetWorkingBase(fastafilename) + ".overview"	
else:
	outputfilename=GetWorkingBase(classificationfilename) + ".overview"	
classificationfile=open(classificationfilename)
header=next(classificationfile)
classificationfile.close()
p_s=-1
p_g=-1
p_f=-1
p_o=-1
p_c=-1
p_p=-1
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
	i=i+1           
outputfile=open(outputfilename,"w")
speciesnumber,speciesseqnumber,count=ReportAtLevel(p_s,seqids)
genusnumber,genusseqnumber,count=ReportAtLevel(p_g,seqids)
familynumber,familyseqnumber,count=ReportAtLevel(p_f,seqids)
ordernumber,orderseqnumber,count=ReportAtLevel(p_o,seqids)
classnumber,classseqnumber,count=ReportAtLevel(p_c,seqids)
phylumnumber,phylumseqnumber,count=ReportAtLevel(p_p,seqids)
if len(seqids) >0:
	count=len(seqids)
print("Number of sequences: " + str(count))
print("Taxonomic level\tNumber of taxa\tNumber of sequences")
print("Species" + "\t" + str(speciesnumber) + "\t" + str(speciesseqnumber))
print("Genus" + "\t" + str(genusnumber) + "\t" + str(genusseqnumber))
print("Family" + "\t" + str(familynumber) + "\t" + str(familyseqnumber))
print("Order" + "\t" + str(ordernumber) + "\t" + str(orderseqnumber))
print("Class" + "\t" + str(classnumber) + "\t" + str(classseqnumber))
print("Phylum" + "\t" + str(phylumnumber) + "\t" + str(phylumseqnumber))	
outputfile.write("Number of sequences: " + str(count) + "\n")
outputfile.write("Taxonomic level\tNumber of taxa\tNumber of sequences\n")
outputfile.write("Species" + "\t" + str(speciesnumber) + "\t" + str(speciesseqnumber) + "\n")
outputfile.write("Genus" + "\t" + str(genusnumber) + "\t" + str(genusseqnumber) + "\n")
outputfile.write("Family" + "\t" + str(familynumber) + "\t" + str(familyseqnumber) + "\n")
outputfile.write("Order" + "\t" + str(ordernumber) + "\t" + str(orderseqnumber) + "\n")
outputfile.write("Class" + "\t" + str(classnumber) + "\t" + str(classseqnumber) + "\n")
outputfile.write("Phylum" + "\t" + str(phylumnumber) + "\t" + str(phylumseqnumber) + "\n")

outputfile.close()
print("The overview is saved in  file " + outputfilename + ".")
