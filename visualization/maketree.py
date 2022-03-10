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
							   usage="%(prog)s [options] -i fastafile -c classificationfilename -rank classification rank",
							   description='''The script create a phylogenetic tree for the sequences given in the input fasta file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta input file.')
parser.add_argument('-c','--classification', default="", help='The taxonomic classification file.')
#parser.add_argument('-p','--classificationposition', default="", help='the classification positions for the prediction.')
parser.add_argument('-ranks','--classificationranks', default="", help='the classification ranks for getting sequence dscriptions.')
parser.add_argument('-o','--out', default="dnabarcoder", help='The output folder.')

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

def LoadClassification(classificationfilename,poslist):
	classificationfile= open(classificationfilename, errors="ignore")
	seqids=[]
	classifications=[]
	numberoffeatures=0
	for line in classificationfile:
		elements=line.rstrip().split("\t")
		seqid = elements[0].replace(">","").rstrip()
		seqids.append(seqid)
		classification=(line[line.index("\t")+1:])
		texts=classification.split("\t")
		classification=""
		if len(poslist) >0:		
			for pos in poslist:
				pos=pos-1
				text=""
				if pos < len(texts):
					text=texts[pos]
				if text=="":
					text="unidentified"	
				text=text.replace(" ","_")	
				classification=classification + text + "|"
		else:
			for text in texts:
				text=text.rstrip()
				if text=="":
					text="unidentified"
				text=text.replace(" ","_")		
				classification=classification + text + "|"
		classification=classification[:-1] 		
		classifications.append(classification)
		if numberoffeatures==0:
			numberoffeatures=classification.count("|")
		else:
			numberoffeatures=min(numberoffeatures,classification.count("|"))
	classificationfile.close()
	return seqids,classifications,numberoffeatures

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
		#command="iqtree -pers 0.2 -s " + alignmentfilename
		command="iqtree -pers 0.2 -n 500 -s " + alignmentfilename
		#command="iqtree -s " + alignmentfilename
		#command="fasttree -nt -quote " + alignmentfilename + " > " + treefilename
		print(command)
		os.system(command)
	#print tree
	PrintTree(treefilename)
	return treefilename
def CreateFastaFileWithClassification(fastafilename,classificationfilename,positionlist):
	newfastafilename=GetWorkingBase(fastafilename) + ".classification.fasta"
	seqids,classifications,numberoffeatures=LoadClassification(classificationfilename,positionlist)
	outputfile=open(newfastafilename,"w")	
	fastafile=open(fastafilename)
	for line in fastafile:
		if line.startswith(">"):
			seqid=line.rstrip().replace(">","")
			if "|" in seqid:
				seqid=seqid[0:seqid.index("|")]
			classification=""
			if seqid in seqids:
				classification=unicode(classifications[seqids.index(seqid)])
			else:
				classification="|" * numberoffeatures
			header=">" + seqid + "|" + classification + "\n"
			outputfile.write(header)	
		else:		
			outputfile.write(line)
	outputfile.close()		
	fastafile.close()
	return newfastafilename

def GetPositionList(classificationfilename,ranks):
	ranklist=[]	
	if "," in ranks:
		ranklist=ranks.split(",")
	elif ranks !="":
		ranklist.append(ranks)
	positionlist=[]
	classificationfile=open(classificationfilename)
	header=classificationfile.readline()
	header=header.rstrip()
	classificationfile.close()
	texts=header.split("\t")
	isError=False
	for rank in ranklist:
		if rank in texts:
			pos=texts.index(rank)
			positionlist.append(pos)
		else:
			print("The rank " + rank + " is not given in the classification." )
			isError=True
	return positionlist,ranklist,isError

######MAIN
positionlist,ranklist,isError=GetPositionList(classificationfilename,ranks)
if isError==True:
	sys.exit()
newfastafilename=fastafilename
if os.path.exists(classificationfilename):
	newfastafilename=CreateFastaFileWithClassification(fastafilename,classificationfilename,positionlist)
treefilename=CreateTree(newfastafilename)

