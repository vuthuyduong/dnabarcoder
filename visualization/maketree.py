#!/usr/bin/env python
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
from Bio import SeqIO
from Bio import Phylo
import pylab
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='maketree.py',  
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -ranks classificationranks",
							   description='''The script create a phylogenetic tree for the sequences given in the input fasta file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta input file.')
parser.add_argument('-c','--classification', default="", help='The taxonomic classification file.')
#parser.add_argument('-p','--classificationposition', default="", help='the classification positions for the prediction.')
parser.add_argument('-rank','--classificationranks', default="", help='the classification ranks for getting sequence dscriptions.')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence id in the classification file.')

args=parser.parse_args()
fastafilename= args.input
classificationfilename= args.classification
ranks=args.classificationranks
outputpath=args.out

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

def LoadClassification(seqids,classificationfilename,poslist,seqidpos):
	classificationfile= open(classificationfilename, errors="ignore")
	classificationdict={}
	for line in classificationfile:
		elements=line.rstrip().split("\t")
		seqid = elements[seqidpos].rstrip()
		if not (seqid in seqids):
			continue
		classification=""
		for pos in poslist:
			classification=classification + elements[pos].rstrip() + "|"
		classification=classification[:-1] 		
		classificationdict.setdefault(seqid,classification)
	classificationfile.close()
	return classificationdict

def LoadClassificationFromDescription(seqrecords,ranklist):
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
		classification=""
		if "kingdom" in ranklist:
			classification=classification + kingdom + "|"			
		if "phylum" in ranklist:
			classification=classification + phylum + "|"				
		if "class" in ranklist:
			classification=classification + bioclass + "|"	
		if "order" in ranklist:
			classification=classification + order + "|"	
		if "family" in ranklist:
			classification=classification + family + "|"	
		if "genus" in ranklist:
			classification=classification + genus + "|"	
		if "species" in ranklist:
			classification=classification + species + "|"	
		if classification!="":
			classification=classification[:-1]
		classificationdict.setdefault(seqid,classification)	
	return classificationdict	

def lookup_by_names(tree):
	names = {}
	for clade in tree.find_clades():
		if clade.name:
			if clade.name in names:
				raise ValueError("Duplicate key: %s" % clade.name)
			names[clade.name] = clade
	return names
#
#
def computeBranchLength(treefilename):
	tree = Phylo.read(treefilename, "newick")
	names = lookup_by_names(tree)
	for name in names:
		clade=names[name]
		print(clade.name + " " + str(clade.branch_length))


def PrintTree_ete(treefilename):
	svgfilename=GetWorkingBase(fastafilename) + ".tree.svg"
	from ete3 import Tree, TreeStyle, NodeStyle
	ts = TreeStyle()
	ts.show_leaf_name = True
	ts.show_branch_length = False
	ts.show_branch_support = True

	ns = NodeStyle()
	ns["size"] = 0
	ns["vt_line_width"] = 2
	ns["hz_line_width"] = 2

	# read tree
	t = Tree(treefilename)
	t.ladderize(1)
	for node in t.traverse():
		node.set_style(ns)
	# print final tree
	t.render(svgfilename, tree_style=ts)
	print("A iq-tree and its svg file are saved in " + treefilename + " and " + svgfilename + ".")	
	
def PrintTree(treefilename):
	svgfilename=GetWorkingBase(fastafilename) + ".tree.svg"
	tree = Phylo.read(treefilename, "newick")
	Phylo.draw(tree,do_show=True)	
	#Phylo.draw_ascii(tree,do_show=False)	
	pylab.savefig(svgfilename,format='svg', bbox_inches='tight', dpi=300)
	print("A iq-tree and its svg file are saved in " + treefilename + " and " + svgfilename + ".")	
	
def CreateTree(fastafilename):
	alignmentfilename =  GetWorkingBase(fastafilename) + ".aligned.fas"
	#make the alignment
	if not os.path.exists(alignmentfilename):
		command="clustalo -i " + fastafilename + " -o " + alignmentfilename
		print(command)
		os.system(command)
	#make tree
	treefilename=GetWorkingBase(fastafilename) + ".aligned.fas.treefile"
	if not os.path.exists(treefilename):
		#command="iqtree -pers 0.2 -s " + alignmentfilename + " -nt " +  str(nproc)
		command="iqtree version > dnabarcoder.iqtree.log"
		os.system(command)
		command="iqtree -pers 0.2 -n 500 -s " + alignmentfilename
		with open('dnabarcoder.iqtree.log') as f:
			if 'IQ-TREE multicore version' in f.read():
				command="iqtree -pers 0.2 -n 500 -s " + alignmentfilename + " -nt " +  str(nproc)
		os.system("rm dnabarcoder.iqtree.log") 
		#command="iqtree -s " + alignmentfilename
		#command="fasttree -nt -quote " + alignmentfilename + " > " + treefilename
		print(command)
		os.system(command)
	#print tree
	PrintTree(treefilename)
	return treefilename

def CreateFastaFileWithClassification(fastafilename,classificationdict):
	newfastafilename=GetWorkingBase(fastafilename) + ".classification.fasta"
	outputfile=open(newfastafilename,"w")	
	fastafile=open(fastafilename)
	for line in fastafile:
		if line.startswith(">"):
			seqid=line.rstrip().replace(">","")
			if " " in seqid:
				seqid=seqid.split(" ")[0]
			classification=""
			if seqid in classificationdict.keys():
				classification=unicode(classificationdict[seqid])
			header=">" + seqid + "|" + classification.replace(" ","_") + "\n"
			outputfile.write(header)	
		else:		
			outputfile.write(line)
	outputfile.close()		
	fastafile.close()
	return newfastafilename

def GetPositionList(classificationfilename,rankslist):
	positionlist=[]
	classificationfile=open(classificationfilename)
	header=classificationfile.readline()
	header=header.rstrip()
	classificationfile.close()
	texts=header.split("\t")
	isError=False
	i=0
	seqidpos=-1
	for text in texts:
		if text.lower()==args.idcolumnname.lower():
			seqidpos=i
		i=i+1
	if 	seqidpos==-1:
		print("Please specify the sequence id columnname by using -idcolumnname.")
		isError=True
	for rank in ranklist:
		if rank in texts:
			pos=texts.index(rank)
			positionlist.append(pos)
		else:
			print("The rank " + rank + " is not given in the classification." )
			isError=True
	return seqidpos,positionlist,isError

######MAIN
ranklist=[]	
if "," in ranks:
	ranklist=ranks.split(",")
elif ranks !="":
	ranklist.append(ranks)
classificationdict={}
seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
if os.path.exists(classificationfilename):
	seqidpos,positionlist,isError=GetPositionList(classificationfilename,ranks)
	if isError==True:
		sys.exit()
	classificationdict=LoadClassification(seqrecords.keys(),classificationfilename,positionlist,seqidpos)	
else:
	classificationdict=LoadClassificationFromDescription(seqrecords,ranklist)
newfastafilename=CreateFastaFileWithClassification(fastafilename,classificationdict)
treefilename=CreateTree(newfastafilename)

