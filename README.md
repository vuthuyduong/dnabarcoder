# Dnabarcoder

We demonstrate that metabarcoding loses significant resolution and scientific explanatory power by relying on a single sequence similarity value for taxonomic assignment by presenting the dnabarcoder software. Dnabarcoder is a tool to <strong> PREDICT </strong> global and local similarity cut-offs for fungal sequence identification for a reference dataset, and <strong> CLASSIFY </strong> unidentified sequences based on the predicted similarity cutoffs. This reference dataset should come in the form of a FASTA file and should contain barcode sequences from as many species as possible. The classification or taxonomy of reference sequences can be given in the sequence headers of the fasta file (see [data/CBSITS_classification.fasta](https://raw.githubusercontent.com/vuthuyduong/dnabarcoder/master/data/CBSITS_classification.fasta) as an example), or an auxiliary file must contain their full taxonomic classification in a tab-delimited way (kingdom, phylum, class, and so on, see the [data/CBSITS.fasta](https://raw.githubusercontent.com/vuthuyduong/dnabarcoder/master/data/CBSITS.fasta) and [data/CBSITS.current.classification](https://raw.githubusercontent.com/vuthuyduong/dnabarcoder/master/data/CBSITS.current.classification) as examples). 

Dnabarcoder consists of five components namely analysis, visualization, prediction, classification, and verification (see the figure below). 

The analysis and visualization components are to analyze and give an overview about the reference sequences.

The prediction component is designed to predict local and global similarity cut-offs for a reference dataset for sequence identification. For big datasets like the UNITE database, it is practical to compute local similarity cut-offs for sequence identification at the immediate higher taxonomic level. To avoid memory constraints, it is practical to provide a maximum number (M) of sequences selected for each clade for the prediction (default =20000). If the number of the sequences exceeds this maximum number M, a number of M sequences will be selected randomly for the prediction.

The classification component is used to classify unidentified sequences (DNA barcodes, ASVs, or OTUs) provided in a FASTA file against the reference dataset with the predicted similarity cut-offs, while the verification component verifies the classification results. These components are described below. For every function in a component of dnabarcoder, a figure is generated automatically to aid in the interpretation of the results.

<img src="https://github.com/vuthuyduong/dnabarcoder/blob/master/images/dnbarcoder_flowchart.jpg" width="1000" height="300">

For every function of dnabarcoder, a figure is generated to interpret the result. An example of a complete workflow of dnabarcoder can be found in file [data/CBSITS2.sh](https://github.com/vuthuyduong/dnabarcoder/blob/master/data/CBSITS2.sh). The explanation for each command of the file is given below.

Although dnabarcoder was initially developed for fungi, it is applicable to any other organisms using DNA barcodes for identification.

## Contact person 

Duong Vu (d.vu@wi.knaw.nl)

## Citation

Duong Vu, R. Henrik Nilsson, Gerard J.M. Verkley (2022). dnabarcoder: an open-source software package for analyzing and predicting DNA sequence similarity cut-offs for fungal sequence identification.  Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13651

## Dependencies:

- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), used for DNA sequence comparisons
- [Biopython](https://biopython.org/), used for DNA sequence analyses
- [Matplot](https://matplotlib.org/), used for the visualization of the results
- [Krona](https://github.com/marbl/Krona/wiki), optional for visualizing classifications
- [LARGEVIS](https://github.com/rugantio/LargeVis-python3), optional for visualization
- [DiVE](https://github.com/NLeSC/DiVE), optional for visualization
- [IQtree](http://www.iqtree.org/), optional for verification
- [MAFFT](https://mafft.cbrc.jp/alignment/software/), optional for verification
- [Clustalo](http://www.clustal.org/omega/), opional for verification
- scikit-learn, numpy (version 1.16.2) and scipy (version 1.2.1), optional for evaluating classification results. If for higher versions of numpy and scipy, python 3.8 is required. 

We can use the following command to install scikit-learn,scipy, and matplotlib. For python 2.7:

 pip install -U scikit-learn scipy matplotlib  
 
For python3:  

pip3 install -U scikit-learn scipy matplotlib

## Installation:

Using github (preferable):

git clone https://github.com/vuthuyduong/dnabarcoder.git

Using pypi:

pip install dnabarcoder

Using conda:

conda install -c duong.t.vu dnabarcoder

## Inputs

As mentioned earlier, most of the functions dnabarcoder requires two files as inputs: a FASTA file containing reference sequences ([data/CBSITS.fasta](https://github.com/vuthuyduong/dnabarcoder/blob/master/data/CBSITS.fasta)) with unique sequence ids, and a classification file in a tab-delimited format containing the taxonomic classification of the sequences ([data/CBSITS.current.classification](https://github.com/vuthuyduong/dnabarcoder/blob/master/data/CBSITS.current.classification)) where the header contains the ranks of the sequences as follows:

id	kingdom	phylum	class	order	family	genus	species	strain number

MH854569	Fungi	Ascomycota	Dothideomycetes			Monodictys	Monodictys castaneae	CBS 100.07

MH854570	Fungi	Ascomycota	Sordariomycetes	Hypocreales	Nectriaceae	Fusarium	Fusarium equiseti	CBS 107.07

We can use [mkCOInr](https://github.com/meglecz/mkCOInr) to obtain this tab-delimited format for a reference database.

The taxonomic classification of the sequences can also be provided in the sequence headers as well where the taxonomic classification should have the following format: k__kingdom;p__phylum;c__class;o__order;f__family;s__species (see [data/CBSITS_classification.fasta](https://github.com/vuthuyduong/dnabarcoder/blob/master/data/CBSITS_classification.fasta))

## Tips for preparing reference datasets  to speed up the prediction and classification

We reduce the number of the reference sequences for the prediction and classification by selecting only unique and identified reference sequences. 

- We can use an aid script to select only unique sequences for the references:

aidscripts/selectsequences.py -i CBSITS.fasta -c CBSITS.current.classification -unique yes -o CBSITS.unique.fasta

- To predict cut-offs for species (for example) identification, we first select only unique sequences that are identified at the species level:

aidscripts/selectsequences.py -i CBSITS.fasta -c CBSITS.current.classification -rank species -unique yes -o CBSITS.species.fasta

## Outputs

Outputs of dnabarcoder will be saved in an output folder, specified by the user. If this output folder is not given, a folder namely dnabarcoder will be created automatically.

## Analysis and Visualization

The analyzation component was to get an overview, and to study the length, the distribution, and the similarity variation of the sequences at different taxonomic levels. All outputs will be saved in an output folder, dnabarcoder is given as default.

- To get an overview of the CBSITS.fasta dataset:

../../dnabarcoder.py overview -i CBSITS.fasta -c CBSITS.current.classification

If the taxonomic classifications of the sequences are given in sequence headers, we use the following command:

../../dnabarcoder.py overview -i CBSITS_classification.fasta

The output is given in the file dnabarcoder/CBSITS.overview. The overview at the species, genus, family, order, and class levels are given in the files dnabarcoder/CBSITS.overview.species, dnabarcoder/CBSITS.overview.genus, dnabarcoder/CBSITS.overview.family, dnabarcoder/CBSITS.overview.order, and dnabarcoder/CBSITS.overview.class 

- To see the distribution of the barcode lengths:

../../dnabarcoder.py length -i CBSITS.fasta -l 100

Here l is the interval length. Analyzing sequence lengths is important to decide the minimum BLAST alignment length ml. Next to a text file, a figure will be generated as an output of the function as follows:

<img src="https://github.com/vuthuyduong/dnabarcoder/blob/master/images/CBSITS.length.png" width="1000" height="300">

- To get the distribution of the sequences at different taxonomic level. In the following example, the distribution of the sequences is computed from the species to the class level:

../../dnabarcoder.py distribution -i CBSITS.fasta -c CBSITS.current.classification -rank class,order,family,genus,species            

Or:

../../dnabarcoder.py distribution -i CBSITS_classification.fasta -rank class,order,family,genus,species

where ranks are the classification ranks that we are interested in. Next to an output text file containing output information, a figure is generated as follows:

<img src="https://github.com/vuthuyduong/dnabarcoder/blob/master/images/CBSITS.distribution.png" width="300" height="300">

If we want to visualize the distribution of the sequences with Krona, then we can use the following command:

../../dnabarcoder.py  distribution -i CBSITS.fasta -c CBSITS.current.classification -rank class,order,family,genus,species -method krona

or:

../../dnabarcoder.py  distribution -i CBSITS_classification.fasta -rank class,order,family,genus,species -method krona


- To get sequence variation with different taxonomic groups:

../../dnabarcoder.py variation -i CBSITS.fasta -c CBSITS.current.classification -rank class,order,family,genus,species  -ml 400

Or:

../../dnabarcoder.py variation -i CBSITS_classification.fasta -rank class,order,family,genus,species  -ml 400

Here the minimum BLAST alignment length ml is set to 400 as 95% of the barcodes have a length of more than 400bp. For short sequences like ITS1 or ITS2, ml should be set to smaller such as 50. Next to an output text file, a figure is generated as follows:

<img src="https://github.com/vuthuyduong/dnabarcoder/blob/master/images/CBSITS.variation.png" width="500" height="300">

- To compute a similarity matrix for the barcodes:

../../dnabarcoder.py sim -i CBSITS.fasta -ml 400

The output is given in the file dnabarcoder/CBSITS.sim. 

## Visualization

The second component of dnabarcoder is to visualize the sequences-based 2D/3D “embeddings” using Matplotlib. Sequences’ coordinates are computed using LargeVis.
Together with sequence distribution and variation, visualization helped evaluate the predicted similarity cut-offs and classifi-cation results. 

- To visualize the sequences, use the following command:

../../dnabarcoder.py visualize -i CBSITS.fasta -c CBSITS.current.classification -rank class -ml 400 -sim dnabarcoder/CBSITS.sim

Note that the minimum BLAST alignment length ml for ITS sequences is set to 400. For short sequences like ITS1 or ITS2, ml should be set to smaller such as 50. 

Here the sequences of the same taxonomic class will have the same color. The output is given below:

<img src="https://github.com/vuthuyduong/dnabarcoder/blob/master/images/CBSITS.3.visualization.png" width="300" height="300">

Here the sequences are colored by on the taxa at the class level. 

If the simmatrix is not given, dnabarcoder will compute it and save it in the file dnabarcoder/CBSITS.sim. Note that if the computer cannot handle the complete similarity matrix, it is better to use [fMLC](https://github.com/vuthuyduong/fMLC/tree/master/Windows) for visualization.

- We can also visualize the sequences using [DIVE](https://github.com/NLeSC/DiVE). In this case, please download DiVE and place in the visualization folder, and use the following command:

../../dnabarcoder.py visualize -i CBSITS.fasta -c CBSITS.current.classification -rank class -ml 400 -method dive

## Prediction

The third component is to cluster and predict a similarity cut-off for sequence identification based on taxonomic classification. Given a taxonomic level, if higher taxonomic levels are not given, then whole dataset will be used for the prediction. The local similarity cut-offs assign more sequences than the global similarity cut-offs, and less computationally expensive to compute.

- To predict a <strong> global similarity cut-off </strong> for genus identification of the CBSITS dataset for example, use the followig command. Note that this action should not be used for a very large dataset as it might take quite sometime to finish.

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -rank genus -ml 400

Or:

../../dnabarcoder.py predict -i CBSITS_classification.fasta -st 0.7 -et 1 -s 0.001 -rank genus -ml 400

<strong> Note that the minimum BLAST alignment length ml for ITS sequences is set to 400. For short sequences like ITS1 or ITS2, ml should be set to smaller such as 50. </strong>

For this action, a complete similarity matrix will be computed if dnabarcoder/CBSITS.sim does not exist. The prediction is saved in the file dnabarcoder/CBSITS.predicted, and the predicted cutoffs are saved in a json format file dnabarcoder/CBSITS.cutoffs.json and a tab delimited format file dnabarcoder/CBSITS.cutoffs.json.txt. Note that if the file dnabarcoder/CBSITS.predicted exists, the soft will not recompute the existing predictions, and the new prediction will appended to the file. 

- <strong> If a similarity matrix of the sequences is pre-calculated, then we can use it for the prediction:</strong> 

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -rank genus -higherrank family,order,class,phylum -ml 400 <strong> -sim dnabarcoder/CBSITS.sim </strong>

The similarity matrix can be given in a tab delimited format file of the following form:

seqid1 seqid1 sim11

seqid1 seqid2 sim12

.

.

.


- If we wish to recompute the existing prediction, please use the following command:

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -rank genus -ml 400 -redo yes

- If we wish the prediction will appended to an existing file with a given prefix such as existingCBSITS for example, please use the following command:

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -rank genus -ml 400 -prefix existingCBSITS


- To predict <strong> local similarity cut-offs </strong> for genus identification of the CBSITS dataset for example, use the followig command:

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -rank genus -higherrank family,order,class,phylum -ml 400

- We can also predict the <strong> local cutoffs </strong> at a given rank individually:

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -rank genus -higherrank family -ml 400

- or in a given taxa:

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -rank genus -higherrank phylum -ml 400 -taxa Ascomycota

- To predict <strong> global similarity cut-offs </strong> for the CBSITS dataset for <strong> species identification </strong>, we first need to <strong> remove </strong> sequences of species complexes that are indistinguishable by ITS with 100% similarity score:

../../dnabarcoder.py remove -i CBSITS.fasta -c CBSITS.current.classification -rank species -sim dnabarcoder/CBSITS.sim -ml 400 -t 1

Here t is the threshold or similarity cut-off for removing sequences of the same complex. The results will be saved in dnabarcoder/CBSITS.species.diff.fasta dna dnabarcoder/CBSITS.similar

- To predict <strong> global similarity cut-off for species identification </strong> of the CBSITS dataset:

../../dnabarcoder.py predict -i dnabarcoder/CBSITS.species.diff.fasta -c CBSITS.current.classification -st 0.9 -et 1 -s 0.001 -rank species -ml 400 -sim dnabarcoder/CBSITS.sim -prefix CBSITS

- To predict <strong> local similarity cut-offs for species identification for the genera </strong> of the CBSITS dataset:

../../dnabarcoder.py predict -i dnabarcoder/CBSITS.species.diff.fasta -c CBSITS.current.classification -st 0.9 -et 1 -s 0.001 -rank species -higherrank genus -ml 400 -sim dnabarcoder/CBSITS.sim -prefix CBSITS 

- For <strong> higher taxonomic levels </strong>:

../../dnabarcoder.py predict -i dnabarcoder/CBSITS.species.diff.fasta -c CBSITS.current.classification -st 0.9 -et 1 -s 0.001 -rank species -higherrank family,order,class,phylum -ml 400 -sim dnabarcoder/CBSITS.sim -prefix CBSITS 

The prediction and cutoffs will be saved in the files dnabarcoder/filamentousfungalITS.predicted, dnabarcoder/CBSITS.cutoffs.json and dnabarcoder/CBSITS.cutoffs.json.txt.
 

- If the similarity matrix is not given because of a large dataset, we can just simply use the following command to remove species complexes and predict local similarity cut-offs for species identification for the genera of the CBSITS dataset:

../../dnabarcoder.py predict -i CBSITS.species.fasta -c CBSITS.current.classification -st 0.9 -et 1 -s 0.001 -rank species -higherrank genus -ml 400 -removecomplexes yes -prefix CBSITS 

 - For <strong> large datasets  such as UNITE datasets </strong>, we can remove species complexes during the prediction:
 
 ../../dnabarcoder.py predict -i dnabarcoder/CBSITS.species.fasta -c CBSITS.current.classification -st 0.9 -et 1 -s 0.001 -rank species -higherrank genus -ml 400 <strong> -removecomplexes yes </strong> -prefix CBSITS
 
 We can also set up a maximum number of sequences loaded for each clade for prediction (default is 20000)
 
  ../../dnabarcoder.py predict -i dnabarcoder/CBSITS.species.fasta -c CBSITS.current.classification -st 0.9 -et 1 -s 0.001 -rank species -higherrank genus -ml 400 -removecomplexes yes -prefix CBSITS <strong> -maxseqno 10000 </strong>

- To <strong> visualize </strong> the global prediction for all ranks, use the following command:

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -rank species,genus,family,order,class

The output is given below:

<img src="https://github.com/vuthuyduong/dnabarcoder/blob/master/images/CBSITS.global.png" width="600" height="300">


- To visualize the local prediction for species identification in the genera of the reference dataset, use the following command:

../../dnabarcoder.py predict -i CBSITS2.fasta -c ITS_20211006.classification -rank species -higherrank genus

The output is given below:

<img src="https://github.com/vuthuyduong/dnabarcoder/blob/master/images/CBSITS.7.local.png" width="600" height="300">

## Computing the best cut-offs 

- To compute <strong> the best cutoffs </strong>:

../../dnabarcoder.py best -i dnabarcoder/CBSITS.cutoffs.json -c CBSITS.current.classification

Or if the classifications are given in the sequence headers:

../../dnabarcoder.py best -i dnabarcoder/CBSITS.cutoffs.json -f CBSITS_classification.fasta

The best similarity cut-offs to classify the sequences at different taxonomic levels for the taxa given in the cut-offs file are saved in json and text format files [dnabarcoder/CBSITS.cutoffs.best.json](https://github.com/vuthuyduong/dnabarcoder/blob/master/data/CBSITS.cutoffs.json) and dnabarcoder/CBSITS.cutoffs.best.txt.
The best similarity cut-offs to assign sequences to the taxa given in the classification file are saved in json and text format files [dnabarcoder/moldITS.cutoffs.assign.json](https://raw.githubusercontent.com/vuthuyduong/dnabarcoder/master/data/CBSITS.cutoffs.assign.json) and dnabarcoder/CBSITS.cutoffs.assign.txt.


- To <strong> merge </strong> two or more similarity cut-offs files, use the following commands. For a taxonomic level and group, the output file will keep the similarity cut-off having the highest confidence:

../../dnabarcoder.py merge -i dnabarcoder/CBSITS.cutoffs.json,dnabarcoder/existing.cutoffs.json -o mergedcutoffs.json


## Classification 

The last component of dnabarcode is to classify a dataset against a reference/barcode dataset using a similarity cut-off or the local cut-offs given by the users or predicted by dnabarcoder for the reference dataset. The similarity cut-offs file should be in json format and have a structure similar to the structure of the file data/CBSITS.cutoffs.json. For classification, dnabarcoder will use the cut-off given in the 'cut-off' tag for classification.

- To <strong> search </strong> for the best match of the sequences in the UNITErelease.fasta file, use the following command:

../../dnabarcoder.py <strong> search </strong> -i UNITErelease.fasta -r CBSITS.fasta -ml 400

The result is saved in the file dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch

- To <strong> classify </strong> the UNITE sequences based on best matches with BLAST, using the following commands:

 - Globally, based on only one similarity cut-off:

../../dnabarcoder.py <strong> classify </strong> -i dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch -c CBSITS.current.classification -cutoff 0.994 -rank species -confidence 0.8334 

Or if classifications are given in sequence headers:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch -r CBSITS_classification.fasta -cutoff 0.994 -rank species -confidence 0.8334 

Here 0.994 is the global similarity cut-off for sequence identification at the species level. The result including unidentified sequences will be saved in dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.species.0994.classified. If we want to save only the classified sequences, use the following command:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch -r CBSITS_classification.fasta -cutoff 0.994 -rank species -confidence 0.8334 -saveclassifiedonly True

- Note that dnabarcoder could also take <strong> the output of BLAST (fmt 6) </strong> as the input for classification. Please use the following command:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.CBSITS_BLAST.blastoutput -c CBSITS.current.classification -cutoffs dnabarcoder/CBSITS.cutoffs.best.json 

- Locally, based on the similarity cut-off predicted for the best match:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch -c CBSITS.current.classification -cutoffs [dnabarcoder/CBSITS.cutoffs.best.json](https://github.com/vuthuyduong/dnabarcoder/blob/master/data/CBSITS.cutoffs.json) 

The result will be saved in dnabarcoder/UNITErelease.CBSITS_BLAST.classified. 

- Or we can also use the cut-offs file to assign the sequences:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch -c CBSITS.current.classification -cutoffs [dnabarcoder/CBSITS.cutoffs.assign.json](https://raw.githubusercontent.com/vuthuyduong/dnabarcoder/master/data/CBSITS.cutoffs.assign.json) 

The result will be saved in dnabarcoder/UNITErelease.CBSITS_BLAST.classified. 

- Only classify at the species level:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch -c fCBSITS.current.classification -cutoffs dnabarcoder/CBSITS.cutoffs  -rank species

The result will be saved in dnabarcoder/UNITErelease.CBSITS_BLAST.species.classified. 

- To compute <strong> classification/assigment accuracy and precision </strong>, use the following commands:

../../dnabarcoder.py accuracy -i dnabarcoder/UNITErelease.CBSITS_BLAST.species.classified -c UNITErelease.current.classification -r CBSITS.current.classification

- To visualize the classification/assignment results with Krona:

../../dnabarcoder.py krona -i dnabarcoder/UNITErelease.CBSITS_BLAST.classified -c CBSITS.current.classification

## Verification

- To <strong> verify </strong> the classification results based on phylogenic trees or cutoffs.

- To verify the classifications based on phylogenic trees at the species level:

../../dnabarcoder.py verify -i dnabarcoder/UNITErelease.CBSITS_BLAST.classified -c CBSITS.current.classification -r CBSITS.fasta -f UNITErelease.fasta -rank species -method tree

- To verify the classifications based on cutoffs:

../../dnabarcoder.py verify -i dnabarcoder/UNITErelease.CBSITS_BLAST.classified -c CBSITS.current.classification -r CBSITS.fasta -f UNITErelease.fasta -rank -cutoffs CBSITS.cutoffs.best.json -method cutoff

## Data

The CBSITS barcode dataset was released in Vu et al. (2019), while the UNITErelease.fasta dataset is the [UNITE general FASTA release](https://plutof.ut.ee/#/doi/10.15156/BIO/786368). The CBSITS.current.classification and UNITErelease.current.classification were updated from [Mycobank](https://www.mycobank.org/).


## Citation

Duong Vu, R. Henrik Nilsson, Gerard J.M. Verkley (2022). dnabarcoder: an open-source software package for analyzing and predicting DNA sequence similarity cut-offs for fungal sequence identification.  Molecular Ecology Resources. https://doi.org/10.1111/1755-0998.13651

## References

The data used in Vu et al. (2022) were released in the following paper:

Vu D. et al. (2019). Large-scale generation and analysis of filamentous fungal DNA barcodes boosts coverage for kingdom fungi and reveals thresholds for fungal species and higher taxon delimitation. Studies in Mycology 92, 135-154.


