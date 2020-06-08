#!/usr/bin/env python
# FILE: convertVariationJson2Tab.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys
import os
import json

variationfilename=sys.argv[1]
output=sys.argv[2]

outputfile=open(output,"w")
variation={}
with open(variationfilename) as variation_file:
	variation = json.load(variation_file,encoding='latin1')	
outputfile.write("Taxonname\tThreshold\tMin threshold\tNumber of sequences\n")	
for classname in variation.keys():
	classname=classname.encode('ascii', 'ignore')
	#classname=classname.encode('ascii', 'ignore')
	#print(classname)
	if classname not in variation.keys():
		continue
	threshold=variation[classname]['Threshold']
	minthreshold=variation[classname]['MinThreshold']
	seqno=variation[classname]['NumberOfSequences']
	outputfile.write(classname + "\t" + str(threshold) + "\t" + str(minthreshold) + "\t" + str(seqno) + "\n")
outputfile.close()
  
