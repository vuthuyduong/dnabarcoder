#!/usr/bin/env python
# FILE: fasta2tab.py
# AUTHOR: Duong Vu
# CREATE DATE: 07 June 2019
import sys, argparse
if sys.version_info[0] >= 3:
	unicode = str

from Bio import SeqIO
import multiprocessing
nproc=multiprocessing.cpu_count()

parser=argparse.ArgumentParser(prog='fasta2tab.py',  
							   usage="%(prog)s [options] -i fasta -o output",
							   description='''Script that creates a tab. file from a fasta file. ''',
							   epilog="""Written by Duong Vu duong.t.vu@gmail.com""",
   )

parser.add_argument('-i','--input', required=True, help='the tab. file.')
parser.add_argument('-o','--out', help='The fasta filename.') 


args=parser.parse_args()
fastafilename= args.input

output=args.out

def GetBase(filename):
	return filename[:-(len(filename)-filename.rindex("."))]

sequences=SeqIO.to_dict(SeqIO.parse(fastafilename, "fasta"))
outputfile=open(output,"w")
outputfile.write("ID\tSequence\n")
for seqid in sequences.keys():
	sequence=sequences[seqid]
	outputfile.write(seqid + "\t" + str(sequence.seq) + "\n")
outputfile.close()		
