#!/usr/bin/env python
# FILE: replaceSeqID.py, requiring Python 2.7 if classification file contains non-ascii characters
# AUTHOR: Duong Vu
# CREATE DATE: 07 Oct 2024
import sys, argparse
import re
import os

parser=argparse.ArgumentParser(prog='replaceSeqID.py',  
							   usage="%(prog)s [options] -i inputfile -c classificationfile -o output",
							   description='''Script that replaces sequence ids with new ids given in the classification file ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the input file (fasta, tree files).')
parser.add_argument('-c','--classification', required=True, help='the tab. file containing sequence IDs and new IDS.')
parser.add_argument('-o','--out', help='The output file with sequence IDs being replaced with the new IDs.')
#parser.add_argument('-sep','--sep', help='The list of separators, separated by ",", to separate sequence IDs in the input file.')
parser.add_argument('-idcolumnname','--idcolumnname',default="ID", help='the column name of sequence ids in the newidfile.')
parser.add_argument('-newidcolumnname','--newidcolumnname',default="newID", help='the column name of new sequence ids in the newidfile.')

args=parser.parse_args()
inputfilename= args.input
output=args.out

def LoadIDDict(idfilename):
	iddict={}
	idfile=open(idfilename)
	header=next(idfile)
	header=header.rstrip()
	texts=header.split("\t")
	p_id=-1
	p_newid=-1
	i=0
	for text in texts:
		if text==args.idcolumnname:
			p_id=i
		if text==args.newidcolumnname:
			p_newid=i	
		i=i+1	
	for line in idfile:
		line=line.rstrip()
		texts=line.split("\t")
		iddict.setdefault(texts[p_id],texts[p_newid])		
	return iddict
	
###### MAIN####
iddict=LoadIDDict(args.classification)

with open(inputfilename, 'r') as inputfile:
    content = inputfile.read()
for old_id in iddict.keys():	
	new_id=iddict[old_id]
	if old_id==new_id:
		continue
	# Define the pattern to match text between > and " "
	pattern = rf"(?<=>){old_id}(?= )" #sequence headers in fasta files
	# Use re.sub to replace the matched pattern
	content = re.sub(pattern, new_id, content)	
	# Define the pattern to match text between > and \n
	pattern = rf"(?<=>){old_id}(?=\n)" #sequence IDs in fasta files
	# Use re.sub to replace the matched pattern
	content = re.sub(pattern, new_id, content)	
    # Define the pattern to match text between ( and :
	pattern = rf"(?<=\(){old_id}(?=:)" #sequence IDs in tree files
	# Use re.sub to replace the matched pattern
	content = re.sub(pattern, new_id, content)	
    # Define the pattern to match text between , and :
	pattern = rf"(?<=,){old_id}(?=:)" #sequence IDs in tree files
	# Use re.sub to replace the matched pattern
	content = re.sub(pattern, new_id, content)	    

with open(output, 'w') as outputfile:
    outputfile.write(content)
print("The new file with new IDs is saved in " + output + ".")
