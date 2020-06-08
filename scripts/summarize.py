#!/usr/bin/env python
# FILE: summarize.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019

import sys
import os
from Bio import SeqIO

fastafilename = sys.argv[1]
classificationfilename = sys.argv[2] # the classification file name
outputfilename = sys.argv[3] # the output containing taxonomic information given in the sequence headers

def ReportAtLevel(level,seqids):
	classificationfile = open(classificationfilename)
	taxa=[]
	seqNumbers=[]
	seqNo=0
	for line in classificationfile:
		words=line.split("\t")
		seqid=words[0].replace(".1","").replace(">","")
		if seqid not in seqids:
			continue
		if level >=  len(words):
			continue
		name=line.split("\t")[level].rstrip()
		if name=="":
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
	return numberoftaxa,seqNo

######MAIN
seqrecords = list(SeqIO.parse(fastafilename, "fasta"))
seqids=[]
for rec in seqrecords:
	seqids.append(rec.id)
classificationfile=open(classificationfilename)
line=classificationfile.readline()
classificationfile.close()
p_s=1
p_g=2
p_f=3
p_o=4
p_c=5
p_p=6
if line.startswith("#"):
   words=line.split("\t")
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
speciesnumber,speciesseqnumber=ReportAtLevel(p_s,seqids)
genusnumber,genusseqnumber=ReportAtLevel(p_g,seqids)
familynumber,familyseqnumber=ReportAtLevel(p_f,seqids)
ordernumber,orderseqnumber=ReportAtLevel(p_o,seqids)
classnumber,classseqnumber=ReportAtLevel(p_c,seqids)
phylumnumber,phylumseqnumber=ReportAtLevel(p_p,seqids)
outputfile.write("Number of sequences: " + str(len(seqrecords)) + "\n")
outputfile.write("Taxonomic level\tNumber of taxa\tNumber of sequences\n")
outputfile.write("Species" + "\t" + str(speciesnumber) + "\t" + str(speciesseqnumber) + "\n")
outputfile.write("Genus" + "\t" + str(genusnumber) + "\t" + str(genusseqnumber) + "\n")
outputfile.write("Family" + "\t" + str(familynumber) + "\t" + str(familyseqnumber) + "\n")
outputfile.write("Order" + "\t" + str(ordernumber) + "\t" + str(orderseqnumber) + "\n")
outputfile.write("Class" + "\t" + str(classnumber) + "\t" + str(classseqnumber) + "\n")
outputfile.write("Phylum" + "\t" + str(phylumnumber) + "\t" + str(phylumseqnumber) + "\n")

outputfile.close()
