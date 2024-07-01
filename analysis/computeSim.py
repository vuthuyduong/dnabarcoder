#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 28 12:28:09 2021

@author: Duong Vu
"""
import sys
if sys.version_info[0] >= 3:
	unicode = str
import os, argparse
import multiprocessing
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font',size=6)
from matplotlib.patches import Polygon

nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='computeSim.py',  
							   usage="%(prog)s [options] -i fastafile -o output",
							   description='''Script that computes BLAST similarity matrix for sequences of a fasta file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the fasta file')
parser.add_argument('-o','--out',default="dnabarcoder", help='The output folder.') 
parser.add_argument('-ml','--minalignmentlength', type=int, default=400, help='Minimum sequence alignment length required for BLAST. For short barcode sequences like ITS2 (ITS1) sequences, minalignmentlength should be set to smaller, 50 for instance.')
parser.add_argument('-ms','--minsim', type=float, default=0, help='The minimum similarity score that will be saved for the output.')

args=parser.parse_args()
fastafilename= args.input
mincoverage=args.minalignmentlength
minsim=args.minsim
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


def BoxPlotAll(scorelist,filename):
	figoutput=GetWorkingBase(filename)+".boxplot.png"
	title=GetBase(os.path.basename(filename))
	data=[np.sort(np.array(scorelist))]
	minthreshold=round(float(np.min(np.array(scorelist))),4)
	threshold=round(float(np.median(np.array(scorelist))),4)
	maxthreshold=round(float(np.max(np.array(scorelist))),4)
	averthreshold=round(float(np.average(np.array(scorelist))),4)
#	fig, ax = plt.subplots(figsize=(10, 6))
#	#fig.canvas.set_window_title('Variation')
#	fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)
	fig, ax = plt.subplots(figsize=(3,3))
	box_colors = ['b']#['darkkhaki', 'royalblue']
	bp = ax.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
	plt.setp(bp['boxes'], color='black')
	plt.setp(bp['whiskers'], color='black')
	plt.setp(bp['fliers'], color='red', marker='+')
	
	ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)

	# Hide these grid behind plot objects
	ax.set_axisbelow(True)
	ax.set_title('Median similarity score of ' + title)
	#ax.set_xlabel('')
	ax.set_ylabel('Similarity score')
	num_boxes=len(data)
	medians=np.empty(num_boxes)
	for i in range(num_boxes):
		box=bp['boxes'][i]
		boxX=[]
		boxY=[]
		for j in range(5):
			boxX.append(box.get_xdata()[j])
			boxY.append(box.get_ydata()[j])
		box_coords = np.column_stack([boxX, boxY])	
		# Fill in the color
		ax.add_patch(Polygon(box_coords, facecolor=box_colors[i % 2]))
		med = bp['medians'][i]
		medianX=[]
		medianY=[]
		for j in range(2):
			medianX.append(med.get_xdata()[j])
			medianY.append(med.get_ydata()[j])
		medians[i] = medianY[0]
		#plot the average value
		ax.plot(np.average(med.get_xdata()), np.average(data[i]),color='w', marker='*', markeredgecolor='k')
	#add labels	
	ax.set_xticklabels(np.array(['']))	
	#add median values 
	upper_labels = [str(np.round(s, 4)) for s in medians]
	pos = np.arange(num_boxes) + 1
	k=0
	for tick, label in zip(range(num_boxes), ax.get_xticklabels()):
		ax.text(pos[tick], 0.97, upper_labels[tick], transform=ax.get_xaxis_transform(), horizontalalignment='center', size='x-small', color=box_colors[k])
		k=k+1
	#plt.legend()
	plt.tight_layout()
	plt.rcParams['font.size'] = 6.0
	plt.savefig(figoutput, dpi = 500)
# 	if displayed==True:
# 		if args.display=="yes":
#	plt.show()
	return figoutput,minthreshold,threshold,maxthreshold,averthreshold
			
def ComputeScoreList(simmatrix,output):
	scorelist=[]
	keys=list(simmatrix.keys())
	for i in range(0,len(keys)-2):
		try:
			list_i=simmatrix[keys[i]]
		except KeyError:
			continue
		for j in range(i+1,len(keys)-1):
			try:
				score=list_i[keys[j]]
				scorelist.append(score)
			except KeyError:
				continue
	figoutput,minthreshold,threshold,maxthreshold,averthreshold=BoxPlotAll(scorelist,output)
	print(len(scorelist))	
	return figoutput,minthreshold,threshold,maxthreshold,averthreshold	
				
def SaveSim(simmatrix,simfilename,ms):
	simfile=open(simfilename,"w")
	for i in simmatrix.keys():
		for j in simmatrix[i].keys():
			if simmatrix[i][j] >= ms:
				simfile.write(str(i) + " " + str(j) + " " + str(simmatrix[i][j]) + "\n")
	simfile.close()

def ComputeSim(fastafilename,mincoverage,minsim):
	#scorematrix = [[0 for x in range(len(seqrecords))] for y in range(len(seqrecords))] 
#	simfilename=GetBase(fastafilename) + ".sim"
#	simfile=open(simfilename,"w")
	#blast
	makedbcommand = "makeblastdb -in " + fastafilename + " -dbtype \'nucl\' " +  " -out db"
	os.system(makedbcommand)
	blastcommand = "blastn -query " + fastafilename + " -db  db -task blastn-short -outfmt 6 -out out.txt -num_threads " + str(nproc)
	if mincoverage >=400:
		blastcommand = "blastn -query " + fastafilename + " -db  db -outfmt 6 -out out.txt -num_threads " + str(nproc)
	os.system(blastcommand)
	
	seqrecords=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
	simmatrix={}
	for seqid in seqrecords.keys():
		simmatrix.setdefault(seqid,{})
		simmatrix[seqid][seqid]=1
	
	#read blast output
	blastoutputfile = open("out.txt")
	for line in blastoutputfile:
		if line.rstrip()=="":
			continue
		words = line.split("\t")
#		queryid = words[0].rstrip()
#		refid = words[1].rstrip()
		i = words[0].rstrip()
		j = words[1].rstrip()
#		i = seqids.index(queryid)
#		j = seqids.index(refid)
		pos1 = int(words[6])
		pos2 = int(words[7])
		iden = float(words[2]) 
		sim=float(iden)/100
		coverage=abs(pos2-pos1)
		score=sim
		if coverage < mincoverage:
			score=float(score * coverage)/mincoverage
		if j in simmatrix[i].keys():
			if simmatrix[i][j] < score:
				simmatrix[i][j]=round(score,4)
				simmatrix[j][i]=round(score,4)
		else:
			simmatrix[i][j]=round(score,4)
			simmatrix[j][i]=round(score,4)
#		if scorematrix[i][j] < score:
#			scorematrix[i][j]=score
#			scorematrix[j][i]=score
	os.system("rm out.txt")
	return simmatrix

if __name__ == "__main__":
	#load sequences
	base = sys.argv[0][0: sys.argv[0].rindex("/")+1]
	simmatrix={}
	print("Computing similarity matrix using BLAST..")
	simmatrix=ComputeSim(fastafilename,mincoverage,minsim)
	output=GetWorkingBase(fastafilename) + ".sim"
	#save the simmatrix
	SaveSim(simmatrix,output,minsim)
	print("The similarity file is saved in " + output + ".")	
	figoutput,minthreshold,threshold,maxthreshold,averthreshold=ComputeScoreList(simmatrix,output)
	print("The boxplot of similarity scores is saved in " + figoutput + ".")
	print("Min. similarity score: " + str(minthreshold) + ".")	
	print("Median similarity score: " + str(threshold) + ".")	
	print("Max. similarity score: " + str(maxthreshold) + ".")	
	print("Average similarity score: " + str(averthreshold) + ".")	
	
			
		
