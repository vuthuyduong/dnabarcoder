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
parser.add_argument('-seqidpos','--sequenceidposition', type=int,default=0, help='the position of sequence id in the classification file.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename= args.classification
outputpath=args.out
seqidpos=args.sequenceidposition

if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)

def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + basename
	return path

def ReportAtLevel(level,seqids,higherlevel,seqidpos):
	if level <0:
		return 0,0,0
	classificationfile = open(classificationfilename)
	next(classificationfile)
	taxa=[]
	seqNo=0
	count=0
	highertaxa={}
	for line in classificationfile:
		count=count+1
		words=line.split("\t")
		seqid=words[seqidpos].replace(".1","").replace(">","")
		if len(seqids) > 0 and (seqid not in seqids):
			continue
		if level >=  len(words):
			continue
		texts=line.split("\t")
		taxon=""
		highertaxon=""
		if level >-1:
			taxon=texts[level].rstrip()
		if higherlevel >-1:	
			highertaxon=texts[higherlevel].rstrip()
		if taxon=="" or taxon=="unidentified" or ("uncultured" in taxon) or ("_sp_" in taxon):
			continue
		seqNo=seqNo+1
		if not (taxon in taxa):
			taxa.append(taxon)
		if highertaxon!="":
			if not (highertaxon in highertaxa.keys()):
				highertaxa.setdefault(highertaxon,[])
			if not (taxon in highertaxa[highertaxon]):
				highertaxa[highertaxon].append(taxon)
	numberoftaxa=len(taxa)		
	classificationfile.close()
	return numberoftaxa,seqNo,count,highertaxa

def SaveOverview(rank,taxa,outputname):
	outfile=open(outputname,"w")
	outfile.write("Taxon\t"+ rank + " number\n")
	for taxon in taxa.keys():
		outfile.write(taxon + "\t" + str(len(taxa[taxon])) + "\n")
	outfile.close()

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
p_k=-1
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
seqnumber,seqnumber,count,species=ReportAtLevel(seqidpos,seqids,p_s,seqidpos)
speciesnumber,speciesseqnumber,count,genera=ReportAtLevel(p_s,seqids,p_g,seqidpos)
genusnumber,genusseqnumber,count,families=ReportAtLevel(p_g,seqids,p_f,seqidpos)
familynumber,familyseqnumber,count,orders=ReportAtLevel(p_f,seqids,p_o,seqidpos)
ordernumber,orderseqnumber,count,classes=ReportAtLevel(p_o,seqids,p_c,seqidpos)
classnumber,classseqnumber,count,phyla=ReportAtLevel(p_c,seqids,p_p,seqidpos)
phylumnumber,phylumseqnumber,count,kingdoms=ReportAtLevel(p_p,seqids,p_k,seqidpos)
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
SaveOverview("sequence",species,outputfilename + ".species")
print("The overview at the species level is saved in  file " + outputfilename + ".species")
SaveOverview("species",genera,outputfilename + ".genus")
print("The overview at the genus level is saved in  file " + outputfilename + ".genus")
SaveOverview("genus",families,outputfilename + ".family")
print("The overview at the family level is saved in  file " + outputfilename + ".family")
SaveOverview("family",orders,outputfilename + ".order")
print("The overview at the order level is saved in  file " + outputfilename + ".order")
SaveOverview("order",classes,outputfilename + ".class")
print("The overview at the class level is saved in  file " + outputfilename + ".class")
SaveOverview("class",phyla,outputfilename + ".phylum")
print("The overview at the phylum level is saved in  file " + outputfilename + ".phylum")

