#!/usr/bin/env python
# FILE: getSeqLengthDistribution.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020
import os
import sys, argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
plt.rc('font',size=6)
import numpy as np

parser=argparse.ArgumentParser(prog='getSeqLengthDistribution.py',  
							   usage="%(prog)s [options] -i fastafile -l intervallength -out outputname",
							   description='''Script that computes the number of the sequences in increasing intervals of sequence lengths. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta input file.')
parser.add_argument('-l','--intervallength', required=True, type=int, default=100, help='The length of the intervals.')
parser.add_argument('-o','--out',default="dnabarcoder", help='The output folder.')
parser.add_argument('-prefix','--prefix',default="", help='The prefix of the output files.')
parser.add_argument('-label','--label',default="", help='The label to display in the figure.')

args=parser.parse_args()
fastafilename= args.input
il= args.intervallength
outputpath=args.out
prefix=args.prefix
label=args.label

if not os.path.exists(outputpath):
	os.system("mkdir " + outputpath)
	
def GetBase(filename):
	if not ("." in filename):
		return filename
	return filename[:-(len(filename)-filename.rindex("."))] 

def GetWorkingBase(filename):
	basename=os.path.basename(filename)
	if "." in basename:
		basename=basename[:-(len(basename)-basename.rindex("."))] 
	path=outputpath + "/" + basename
	return path

def BarPlot(datasetname,labels,sums):
	x = np.arange(len(labels))  # the label locations
	#width = 0.35  # the width of the bars
	if len(labels) <50:
		fig, ax = plt.subplots(figsize=(3,2))
	else:
		fig, ax = plt.subplots(figsize=(6,2))
	ax.set_ylabel('Number of sequences')
	ax.set_title(datasetname + ": sequence length distribution")
	plt.bar(x,sums)
	ax.set_xticks(x)
	ax.set_xticklabels(labels,rotation=90)
	ax.legend()
	plt.tight_layout()
	plt.rcParams['font.size'] = 6.0
	plt.savefig(figoutput, dpi = 500)
	plt.show()
	
####MAIN
if prefix=="":
	prefix=GetBase(os.path.basename(fastafilename))
outputfilename=GetWorkingBase(prefix) + ".length.txt"
figoutput=GetBase(outputfilename)  + ".png"	
seqrecords = list(SeqIO.parse(fastafilename, "fasta"))
maxlength=0
minlength=0
for seqrecord in seqrecords:
	if len(seqrecord.seq) > maxlength:
		maxlength=len(seqrecord.seq)
	if minlength==0:
		minlength=len(seqrecord.seq) 
	elif len(seqrecord.seq)	< minlength:	
		minlength=len(seqrecord.seq)
outputfile=open(outputfilename,"w")
aver=0
n=int(maxlength/il) + 1
sums=[0] * n
for seqrecord in seqrecords:
	l=len(str(seqrecord.seq))
	i=int(l/il)
	sums[i]= sums[i] + 1
j=0
#save sequence length distribution
outputfile.write("Interval\tNumber of sequences\n")
intervals=[]
labels=[]
for i in range(0,n):
	outputfile.write("[" + str(j) + "," + str(j+il) + ")" + "\t" + str(sums[i]) + "\n")
	intervals.append(j)
	labels.append("[" + str(j) + "," + str(j+il) + ")")
	j=j+il
outputfile.close()
print("The minimum sequence length is " + str(minlength))
print("The maximum sequence length is " + str(maxlength))
print("The distribution and its figure are saved in files " + outputfilename + " and " + figoutput + ".")
#plot
if label=="":
	label=prefix
BarPlot(label,labels,sums)

