#!/usr/bin/env python
import sys
import os

fastafilename = sys.argv[1] # the name of the fasta file
resultfilename = sys.argv[2]
startingId = 0
if len(sys.argv) >3 :
	startingId = sys.argv[3] #starting id that should be added to sequence header


fastafile = open(fastafilename)
resultfile = open(resultfilename, "a")
i=int(startingId)
for line in fastafile:
	if line.startswith('>'):
		resultfile.write(">" + str(i) + "|" + line.rstrip()[1:] + "\n")
		i=i+1
	else:
		resultfile.write(line)    
      
