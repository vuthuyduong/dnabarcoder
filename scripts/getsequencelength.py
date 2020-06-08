#!/usr/bin/env python
# FILE: getsequencelength.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import os
import sys
from Bio import SeqIO

fastafilename = sys.argv[1] # fasta file of unknown sequences
outputfilename = sys.argv[2] # output

seqrecords = list(SeqIO.parse(fastafilename, "fasta"))

outputfile=open(outputfilename,"w")
aver=0
intervals=[0] * 11
for seqrecord in seqrecords:
	l=len(str(seqrecord.seq))
	outputfile.write(str(l)+ "\n")
	aver=aver+l
	if l <100:
		intervals[0]= intervals[0] + 1
	elif l <200:
		intervals[1]= intervals[1] + 1
	elif l <300:
		intervals[2]= intervals[2] + 1
	elif l <400:
		intervals[3]= intervals[3] + 1
	elif l <500:
		intervals[4]= intervals[4] + 1
	elif l <600:
		intervals[5]= intervals[5] + 1
	elif l <700:
		intervals[6]= intervals[6] + 1
	elif l <800:
		intervals[7]= intervals[7] + 1
	elif l <900:
		intervals[8]= intervals[8] + 1	
	elif l <1000:
		intervals[9]= intervals[9] + 1
	else:
		intervals[10]= intervals[10] + 1
aver=float(aver)/float(len(seqrecords))
print(float(aver))
print("number of lengths 1-100: " + str(intervals[0]))
print("number of lengths 100-200: " + str(intervals[1]))
print("number of lengths 200-300: " + str(intervals[2]))
print("number of lengths 300-400: " + str(intervals[3]))
print("number of lengths 400-500: " + str(intervals[4]))
print("number of lengths 500-600: " + str(intervals[5]))
print("number of lengths 600-700: " + str(intervals[6]))
print("number of lengths 700-800: " + str(intervals[7]))
print("number of lengths 800-900: " + str(intervals[8]))
print("number of lengths 900-1000: " + str(intervals[9]))
print("number of lengths 1000-: " + str(intervals[10]))
outputfile.close()
