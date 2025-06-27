#!/usr/bin/env python
# FILE: combineCutoffs.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2020


import sys, argparse
import json
#from copy import deepcopy

parser=argparse.ArgumentParser(prog='combineCutoffs.py',  
							   usage="%(prog)s [options] -i listofjson_marker_specific_cutoffs_files -o outputname",
							   description='''The script that lists all cutoffs and confidence measures of different biomarkers for comparison. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='A list of cutoffs filenames separated by commas.')
parser.add_argument('-o','--out', help='The file lists cutoffs and confidence measures of different clades for comparisons.')


args=parser.parse_args()
dictionarylist= args.input
outputfilename=args.out


def SaveCutoffs(mergeddict,dictionarynames,outputfilename):
# 	#save mergeddict		
# 	with open(outputfilename,"w") as json_file:
# 		if sys.version_info[0] >= 3:
# 			json.dump(mergeddict,json_file,indent=2)
# 		else:
# 			json.dump(mergeddict,json_file,encoding='latin1',indent=2)			
	#save as tab. format
	textoutput=outputfilename+".txt"
	textfile=open(textoutput,"w")
	header="Rank\tDataset"
	for dictionary in dictionarynames:
		header=header + "\t" + dictionary + " cutoff" 
	for dictionary in dictionarynames:
		header=header + "\t" + dictionary + " confidence"
	textfile.write(header + "\n")
	cutoff=0
	for rank in mergeddict.keys():
		datasets=mergeddict[rank]
		for datasetname in datasets.keys():
			if len(list(datasets[datasetname].keys()))!=len(dictionaries):
				continue
			line=rank + "\t" + datasetname
			for dictionary in dictionaries:
				cutoff = datasets[datasetname][dictionary]["cut-off"]
				line=line + "\t" + str(cutoff)
			for dictionary in dictionaries:
				confidence = datasets[datasetname][dictionary]["confidence"]
				line=line +"\t" + str(confidence)
			textfile.write(line + "\n")
			
	textfile.close()
	print("The outputs are saved in " + outputfilename + ".")
####MAIN####
dictionaries=[]
if "," in dictionarylist:
	dictionaries=dictionarylist.split(",")
else:
	dictionaries.append(dictionarylist)	
mergeddict={}	
for dictionaryname in dictionaries:
	cutoffdict={}
	#load classes
	with open(dictionaryname,encoding='latin1') as json_file:
		cutoffdict = json.load(json_file)
	if cutoffdict!={}:	
		for rank in cutoffdict.keys():
			if not (rank in mergeddict.keys()):
				mergeddict.setdefault(rank,{})
			datasets=cutoffdict[rank]
			for datasetname in datasets.keys():
				dataset=datasets[datasetname]
				if not (datasetname in mergeddict[rank].keys()):
					mergeddict[rank].setdefault(datasetname,{})
				mergeddict[rank][datasetname].setdefault(dictionaryname,dataset)
SaveCutoffs(mergeddict,dictionaries,outputfilename)


