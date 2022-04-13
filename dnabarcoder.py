#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 15:13:18 2021

@author: duong vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
"""
import sys, os, subprocess, inspect
path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#import lib.library as lib

#git_version = lib.git_version()
#base_version = '1.0.0'
#if git_version:
#    version = base_version+'-'+git_version
#else:
#    version = base_version

default_help = """

Usage:       dnabarcoder <command> <arguments>
version:     %s

Description: dnabarcoder is a tool for the analysis, visualization, classification, and predictions of similarity cut-offs for dna barcodes
    
Command:     overview     	                 Get an overview of the dataset 
             length                          Compute length distribution
             distribution                    Compute sequence distribution
             variation                       Compute sequence variation
             sim                             Compute similarity matrix	                          
             visualize                       Visualize the sequences
             tree                            Create a phylogenetic tree of the sequences			 
             cluster                         Cluster the sequences
             predict                         Predict similarity cut-offs for the sequences
             best                            Compute best similarity cut-offs for the sequences			 
             remove                          Remove similar sequences of the same complexes based on a give threshold
             search                          Search for best matches of the sequences against a file of reference sequences
             classify                        Classify the sequences to the group of their best match if the score is greater than the given cutoff
             verify                          Verify the assigned sequences based on the phylogenetic tree branch lengths			 
             krona                           Visualize classification results using Krona
             evaluate                        Compute accuracy for classification results
			       
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
        """ #% version
wrongcommand=False		
if len(sys.argv) > 1:
	if sys.argv[1] == 'overview':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script gives an overview about the barcode dataset
    
Arguments:   -i, --input             Fasta file of Dna sequences, required
             -c, --classification   The taxonomic classification file in tab delimited format, required
             -o, --out               The output folder, default= "dnabarcoder"
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'analysis', 'overview.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)	
	elif sys.argv[1] == 'length':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script computes the barcode length distribution
    
Arguments:   -i, --input             Fasta file of Dna sequences, required
             -l, --intervallength    The length of intervals, default=100
             -o, --out               The output folder, default= "dnabarcoder"
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'analysis', 'computeLengthDistribution.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)	
	elif sys.argv[1] == 'distribution':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script computes the distribution of the sequences
    
Arguments:   -i, --input             		Fasta file of Dna sequences, required
             -c, --classification    		The taxonomic classification file in tab delimited format, required
             -ranks, --classificationranks 	The classification ranks for computing sequence distribution, required			 
             -n, --numberofdisplayedlabels  The number of the largest taxa to be displayed in the figure, default=8	
             -method, --visualizationmethod The visualization method: krona or plot, default=plot		 			 
             -o, --out                      The output folder, default= "dnabarcoder"
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'analysis', 'computeDistribution.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)		
	elif sys.argv[1] == 'variation':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script computes the variation of the sequences
    
Arguments:   -i, --input             	    FASTA file of Dna sequences, required
             -c, --classification    	    The taxonomic classification file in tab delimited format, required
             -ranks, --classificationranks  The classification ranks for computing sequence variation, default=""		 
             -ml, --minalignmentlength      Minimum sequence alignment length required for BLAST, default=400. For short barcode sequences like ITS2 (ITS1) sequences, ml should be set to smaller, 50 for instance.	
             -m, --maxSeqNo                 The maximum number of randomly selected sequences to be computed in the case the groups are too big, default=0 (no maximum is given).		 	
             -plt, --plottype               The type of visualization: boxplot or plot, default=boxplot
             -o, --out                      The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'analysis', 'computeVariation.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)
	elif sys.argv[1] == 'sim':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script computes BLAST similarity matrix for sequences
    
Arguments:   -i, --input             	   Fasta file of Dna sequences, required
             -ml, --minalignmentlength      Minimum sequence alignment length required for BLAST, default=400. For short barcode sequences like ITS2 (ITS1) sequences, ml should be set to smaller, 50 for instance.	
             -ms, --minsim                 The minimum similarity that will be saved.
             -o, --out                     The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'analysis', 'computeSim.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)	
	elif sys.argv[1] == 'visualize':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script visualizes the sequences based on BLAST similarity 
    
Arguments:   -i, --input             	   Fasta file of Dna sequences, required
             -c, --classification    	   The taxonomic classification file in tab delimited format, required
             -rank, --classificationrank    The classification rank for coloring the sequences
             -coord, --coordinates          The coordinates file in json format if exists.
             -sim, --simfilename            The similarity file if exists.				 
             -ml, --minalignmentlength      Minimum sequence alignment length required for BLAST, default=400. For short barcode sequences like ITS2 (ITS1) sequences, ml should be set to smaller, 50 for instance.	
             -ms, --minsim                  The minimum similarity that will be saved for the similarity matrix, default=0.5			 
             -dim, --dimension              The dimension 2D or 3D for visualization, default=3
             -kneigh, --kneigbors           The k-neighbors number for visualization, default=150
             -size, --size                  The size of the dot
			 -method, --visualizationmethod The visualization method: DiVE and plot, default=plot			
             -o, --out                      The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'visualization', 'visualize.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)	
	elif sys.argv[1] == 'tree':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script visualizes the sequences based on BLAST similarity 
    
Arguments:   -i, --input             	   	     FASTA file of Dna sequences, required
             -c, --classification    	         The taxonomic classification file in tab delimited format 
             -ranks, --classificationranks       The classification ranks for getting sequence descriptions
             -o, --out                           The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'visualization', 'maketree.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)				
	elif sys.argv[1] == 'cluster':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script cluster the sequences with a given threshold (cutoff)
    
Arguments:   -i, --input             	        FASTA file of Dna sequences, required
             -t, --cutoff            	        The threshold (cutoff) used for clustering, default=0.97            
             -sim, --simfilename                The similarity file if exists				 
             -ml, --minalignmentlength          Minimum sequence alignment length required for BLAST in case simmatrix is not given, default=400. For short barcode sequences like ITS2 (ITS1) sequences, ml should probably be set to smaller, 50 for instance.	
             -c, --classification    	        The taxonomic classification file in tab delimited format, optional. This is used to compute the confidence of clustering.
             -rank, --classificationrank        The classification rank to compute the confidence measure for clustering, optional
             -o, --out                          The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'prediction', 'cluster.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)		
	elif sys.argv[1] == 'remove':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script removes sequences of the same complexes based on a given threshold (cutoff)
    
Arguments:   -i, --input             	        FASTA file of Dna sequences, required
             -t, --cutoff            	        The threshold (cutoff) used for clustering, default=0.97            
             -sim, --simfilename                The similarity matrix file if exists				 
             -ml, --minalignmentlength          Minimum sequence alignment length required for BLAST in case simmatrix is not given, default=400. For short barcode sequences like ITS2 (ITS1) sequences, ml should probably be set to smaller, 50 for instance.	
             -c, --classification    	        The taxonomic classification file in tab delimited format, optional. This is used to compute the confidence of clustering
             -rank, --classificationrank        The classification rank to remove complexes.
             -o, --out                          The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'prediction', 'removeComplexes.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)	
	elif sys.argv[1] == 'predict':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script predicts similarity cut-offs for sequence identification
    
Arguments:   -i, --input             	                   Fasta file of Dna sequences, required       
             -sim, --simfilename                           The similarity matrix file if exists				 
             -ml, --minalignmentlength                     Minimum sequence alignment length required for BLAST in case simmatrix is not given, default=400. For short barcode sequences like ITS2 (ITS1) sequences, ml should probably be set to smaller, 50 for instance.	
             -c, --classification    	                   The taxonomic classification file in tab delimited format 			 
             -ranks, --classificationranks                 The classification ranks for prediction, separated by ","
             -higherranks, --higherclassificationranks     The higher classification ranks for the prediction. If hp=="", the whole dataset will be used for prediction
             -st, --startingthreshold                      The starting threshold of the prediction
             -et, --endthreshold                           The ending threshold of the prediction		
             -s, --step       	           ,               The step for the prediction				  
             -minGroupNo,--minimumgroupnumber              The minimum number of groups for prediction, default=5
             -minSeqNo,--minimumsequencenumber             The minimum number of sequences for prediction, default=50	
             -redo,--redo                                  Redo the prediction if redo !="", default=""		
             -prefix,--prefix                              Prefix of all output files, default as the base of the input file				  
             -o, --out                                     The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl

		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'prediction', 'predict.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)		
	elif sys.argv[1] == 'best':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script predicts similarity cut-offs for sequence identification
    
Arguments:   -i, --input             	        the cutoffs file, required       
             -c, --classification    	        The taxonomic classification file in tab delimited format 			 		            
			 -prefix,--prefix                   Prefix of all output files, default as the base of the input file				  
             -o, --out                          The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl

		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'prediction', 'computeBestCutoffs.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)				
	elif sys.argv[1] == 'search':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script classifies a fasta file of sequences against a reference file of sequences
    
Arguments:   -i, --input             	        Fasta file of Dna sequences to be classified, required       
             -r, --reference                    Fasta file of reference sequences, required    			 
             -ml, --minalignmentlength          Minimum sequence alignment length required for BLAST in case simmatrix is not given, default=400. For short barcode sequences like ITS2 (ITS1) sequences, ml should probably be set to smaller, 50 for instance.	
             -prefix,--prefix                   Prefix of all output files, default as the base of the input file				  
             -o, --out                          The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'classification', 'search.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)				
	elif sys.argv[1] == 'classify':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script classifies the sequences to their best match based on local cutoffs or a given global cutoff. For a sequence, a BLAST similarity score to the sequences of the predicted group is computed. The sequence is assigned if the BLAST similarity score is greater or equal than the similarity cut-off
    
Arguments:   -i, --input             	        The file of classified sequences, required
             -f, --fasta                        The fasta file of classified sequences, required    
             -r, --reference                    The fasta file of reference sequences, required    			 
             -ml, --minalignmentlength          Minimum sequence alignment length required for BLAST in case simmatrix is not given, default=400. For short barcode sequences like ITS2 (ITS1) sequences, ml should probably be set to smaller, 50 for instance.	
             -mp, --minproba                    Only consider the classified sequences with a probability greater or equal minproba
             -m, --maxseqno                     Maximum number of the sequences of the predicted taxon name from the classification file will be selected for the comparison to find the best match. If it is not given, all the sequences will be selected, default=0
             -c, --classification    	        The taxonomic classification file in tab delimited format 			 
             -cutoffs, --cutoffs                The similarity cutoffs file predicted by dnabarcoder predict if exists
             -cutoff, --cutoff                  The similarity cutoff, default=0, only used if the similarity cutoffs file is not given
             -confidence, --confidence          The confidence of the similarity cutoff if exists
             -rank, --rank                      The rank to classify the sequences, default=""
             -prefix,--prefix                   Prefix of all output files, default as the base of the input file				  
             -minGroupNo,--minimumgroupnumber   The minimum number of groups for prediction, default=5
             -minSeqNo,--minimumsequencenumber  The minimum number of sequences for prediction, default=50	
             -o, --out                          The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
	
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'classification', 'classify.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)		
	elif sys.argv[1] == 'verify':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script assigns classified sequences based on local cutoffs or a given global cutoff. For a sequence, a BLAST similarity score to the sequences of the predicted group is computed. The sequence is assigned if the BLAST similarity score is greater or equal than the similarity cut-off
    
Arguments:   -i, --input             	        The file of classified sequences, required
             -f, --fasta                        The fasta file of classified sequences, required    
             -r, --reference                    The fasta file of reference sequences, required    			 
             -m, --maxseqno                     Maximum number of the sequences of the predicted taxon name from the classification file will be selected for the comparison to find the best match. If it is not given, all the sequences will be selected, default=0
             -c, --classification    	        The taxonomic classification file in tab delimited format 			 
             -rank, --rank                      The rank to classify the sequences, default=""
             -prefix,--prefix                   Prefix of all output files, default as the base of the input file				  
             -o, --out                          The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
	
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'classification', 'verify.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)				
	elif sys.argv[1] == 'krona':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script visualizes the classification/assignment results using Krona
    
Arguments:   -i, --input             	        The file of classified/assigned sequences, required
             -c, --classification               The fasta file of classified sequences, required 
             -o, --out                          The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'classification', 'visualizeClassification.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)	
	elif sys.argv[1] == 'evaluate':
		help = """
Usage:       dnabarcoder %s <arguments>
version:     %s

Description: The script compute classification metrices
    
Arguments:   -i, --input             	        The file of classified/assigned sequences, required
             -c, --queryclassification          The classification file of query sequences, required
             -r, --refclassification            The classification file of reference sequences, required
             -o, --out                          The output folder, default= "dnabarcoder"			 
Written by Duong Vu duong.t.vu@gmail.com/d.vu@wi.knaw.nl
		""" # % (sys.argv[1], version)		
		arguments = sys.argv[2:]
		if len(arguments) > 1:
			cmd = os.path.join(path, 'classification', 'evaluate.py')
			arguments.insert(0, cmd)
			exe = sys.executable
			arguments.insert(0, exe)
			subprocess.call(arguments)
		else:
			print(help)
			sys.exit(1)				
	else:
		wrongcommand=True
if len(sys.argv) == 1 or wrongcommand==True:
	print(default_help)
	