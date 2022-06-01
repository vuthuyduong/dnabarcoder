#!/usr/bin/env python
# FILE: overview.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020

#import sys
import os, argparse
from Bio import SeqIO

parser=argparse.ArgumentParser(prog='overview.py',  
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -out outputname",
							   description='''The script that summarizes the number of taxonomic groups of the sequences of the given fasta file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', default="", help='the fasta input file.')
parser.add_argument('-c','--classification', default="", help='The taxonomic classification file.')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')

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

def LoadClassificationFromDescription(seqrecords):
	classificationdict={}
	for seqid in seqrecords.keys():
		description=seqrecords[seqid].description
		species=""
		genus=""
		family=""
		order=""
		bioclass=""
		phylum=""
		kingdom=""
		if " " in description:
			description=description.split(" ")[1]
		texts=description.split("|")
		for text in texts:
			text=text.rstrip()
			taxa=text.split(";")	
			for taxon in taxa:
				if taxon.startswith("k__"):
					kingdom=taxon.replace("k__","")
				elif taxon.startswith("p__"):
					phylum=taxon.replace("p__","")
				elif taxon.startswith("c__"):
					bioclass=taxon.replace("c__","")	
				elif taxon.startswith("o__"):
					order=taxon.replace("o__","")
				elif taxon.startswith("f__"):
					family=taxon.replace("f__","")	
				elif taxon.startswith("g__"):
					genus=taxon.replace("g__","")
				elif taxon.startswith("s__") and (" " in taxon.replace("s__","") or "_" in taxon.replace("s__","")):
					species=taxon.replace("s__","")
					species=species.replace("_"," ")
		classification=[kingdom,phylum,bioclass,order,family,genus,species]
		classificationdict.setdefault(seqid,classification)
	return classificationdict		

def LoadClassification(classificationfilename):
	classificationdict={}
	classificationfile = open(classificationfilename)
	header=next(classificationfile)
	words=header.split("\t")
	p_id=0
	p_s=-1
	p_g=-1
	p_f=-1
	p_o=-1
	p_c=-1
	p_p=-1
	p_k=-1
	i=0
	for word in words:
		if word.rstrip().lower()==args.idcolumnname.lower():
			p_id=i
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
	for line in classificationfile:
		seqid=""
		species=""
		genus=""
		family=""
		order=""
		bioclass=""
		phylum=""
		kingdom=""
		words=line.split("\t")
		if p_id >-1:
			seqid=words[p_id].rstrip()
		if p_s >-1:
			species=words[p_s].rstrip()	
		if p_g >-1:
			genus=words[p_g].rstrip()	
		if p_f >-1:
			family=words[p_f].rstrip()	
		if p_o >-1:
			order=words[p_o].rstrip()
		if p_c >-1:
			bioclass=words[p_c].rstrip()	
		if p_p >-1:
			phylum=words[p_p].rstrip()	
		if p_k >-1:
			kingdom=words[p_k].rstrip()	
		classification=[kingdom,phylum,bioclass,order,family,genus,species]
		if seqid!="":
			classificationdict.setdefault(seqid,classification)	
	classificationfile.close()	   
	return classificationdict
	
			
def ReportAtLevel(seqids,level,higherlevel,classificationdict):
	taxa=[]
	seqNo=0
	count=0
	highertaxa={}
	for seqid in classificationdict.keys():
		count=count+1
		if len(seqids) > 0 and (seqid not in seqids):
			continue
		classification=classificationdict[seqid]
		taxon=seqid
		if level >-1:
			taxon=classification[level].rstrip()
		highertaxon=classification[higherlevel].rstrip()
		if taxon=="" or ("unidentified" in taxon) or ("uncultured" in taxon) or ("_sp_" in taxon):
			continue
		seqNo=seqNo+1
		if not (taxon in taxa):
			taxa.append(taxon)
		if not( highertaxon=="" or ("unidentified" in highertaxon) or ("uncultured" in highertaxon)):
			if not (highertaxon in highertaxa.keys()):
				highertaxa.setdefault(highertaxon,[])
			if not (taxon in highertaxa[highertaxon]):
				highertaxa[highertaxon].append(taxon)
	numberoftaxa=len(taxa)		
	return numberoftaxa,seqNo,count,highertaxa

def SaveOverview(rank,taxa,outputname):
	outfile=open(outputname,"w")
	outfile.write("Taxon\t"+ rank + " number\n")
	for taxon in taxa.keys():
		outfile.write(taxon + "\t" + str(len(taxa[taxon])) + "\n")
	outfile.close()

######MAIN################################################################

outputfilename=""
seqids=[]
if fastafilename != "":
	seqrecords = SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	seqids=seqrecords.keys()
	outputfilename=GetWorkingBase(fastafilename) + ".overview"	
else:
	outputfilename=GetWorkingBase(classificationfilename) + ".overview"	
classificationdict={}
if classificationfilename!="":
	classificationdict=LoadClassification(classificationfilename)	
else:
	classificationdict=LoadClassificationFromDescription(seqrecords)	
outputfile=open(outputfilename,"w")
seqnumber,seqnumber,count,species=ReportAtLevel(seqids,-1,6,classificationdict)
speciesnumber,speciesseqnumber,count,genera=ReportAtLevel(seqids,6,5,classificationdict)
genusnumber,genusseqnumber,count,families=ReportAtLevel(seqids,5,4,classificationdict)
familynumber,familyseqnumber,count,orders=ReportAtLevel(seqids,4,3,classificationdict)
ordernumber,orderseqnumber,count,classes=ReportAtLevel(seqids,3,2,classificationdict)
classnumber,classseqnumber,count,phyla=ReportAtLevel(seqids,2,1,classificationdict)
phylumnumber,phylumseqnumber,count,kingdoms=ReportAtLevel(seqids,1,0,classificationdict)
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

