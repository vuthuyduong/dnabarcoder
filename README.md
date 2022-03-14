# Dnabarcoder

Dnabarcoder is a tool to predict global and local similarity cut-offs for fungal sequence identification. It was implemented in Python which takes DNA barcodes in a fasta file  and their taxonomic classification (at the species, genus, family, order, class, etc. levels) in a tab delimited file as inputs (see data/filamentousfungalITS.fasta and data/filamentousfungalITS.current.classification for the format of the files). The output of dnabarcode will be saved in an folder given by the user. If this folder is not given, the folder namely dnabarcoder will be created. Dnabarcoder contains four components: analysis, visualization, prediction, and classification to help analyze and predict similarity cut-offs for a dataset of barcodes as well as its subclades, and to classify a dataset against the barcode dataset with the predicted cut-offs. For every function of dnabarcoder, a figure is generated to interpret the result. An example of a complete workflow of dnabarcoder can be found in file data/CBSITS2.sh.

Although dnabarcoder was initially developed for fungi, it is applicable to any other organisms using DNA barcodes for identification.

## Dependencies:

- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), used for DNA sequence comparisons
- [Matplot](https://matplotlib.org/), used for the visualization of the results
- [Krona](https://github.com/marbl/Krona/wiki), optional for visualizing classifications
- [LARGEVIS](https://github.com/rugantio/LargeVis-python3), optinal for visualization
- [DiVE](https://nlesc.github.io/DiVE), optional for visualization
- [IQtree](http://www.iqtree.org/), optional for verification
- [Clustalo](http://www.clustal.org/omega/), opional for verification

## Analysis and Visualization

The analyzation component was to get an overview, and to study the length, the distribution, and the similarity variation of the sequences at different taxonomic levels. All outputs will be saved in an output folder, dnabarcoder is given as default.

- To get an overview of the moldITS.fasta dataset:

../../dnabarcoder.py overview -i CBSITS.fasta -c CBSITS.current.classification

The output is given in the file dnabarcoder/CBSITS.overview. The overview at the species, genus, family, order, and class levels are given in the files dnabarcoder/CBSITS.overview.species, dnabarcoder/CBSITS.overview.genus, dnabarcoder/CBSITS.overview.family, dnabarcoder/CBSITS.overview.order, and dnabarcoder/CBSITS.overview.class 

- To see the distribution of the barcode lengths:

../../dnabarcoder.py length -i CBSITS.fasta -l 100

Here l is the interval length. Analyzing sequence lengths is important to decide the minimum BLAST alignment length ml. 

-To get the distribution of the sequences at different taxonomic level. In the following example, the distribution of the sequences is computed from the species to the class level:

../../dnabarcoder.py distribution -i CBSITS.fasta -c CBSITS.current.classification -ranks class,order,family,genus,species            

where ranks are the classification ranks that we are interested in.

If we want to visualize the distribution of the sequences with Krona, then we can use the following command:

../../dnabarcoder.py  distribution -i CBSITS.fasta -c CBSITS.current.classification -ranks class,order,family,genus,species -method krona

- To get sequence variation with different taxonomic groups:

../../dnabarcoder.py variation -i CBSITS.fasta -c CBSITS.current.classification -ranks class,order,family,genus,species  -ml 400

Here the minimum BLAST alignment length ml is set to 400 as 95% of the barcodes have a length of more than 400bp. For short sequences like ITS1 or ITS2, ml should be set to smaller such as 50.

- To compute a similarity matrix for the barcodes:

../../dnabarcoder.py sim -i CBSITS.fasta -ml 400

The output is given in the file dnabarcoder/CBSITS.sim

The second component of dnabarcoder is to visualize the sequences-based 2D/3D “embeddings” using Matplotlib. Sequences’ coordinates are computed using LargeVis.
Together with sequence distribution and variation, visu-alization helped evaluate the predicted similarity cut-offs and classifi-cation results. 

../../dnabarcoder.py visualize -i CBSITS.fasta -c CBSITS.current.classification -rank class -ml 400 -sim dnabarcoder/CBSITS.sim

If the simmatrix is not given, dnabarcoder will compute it and save it in the file dnabarcoder/CBSITS.sim.

Here the sequences are colored by on the taxa at the class level. 

We can also visualize the sequences using [DIVE](https://github.com/NLeSC/DiVE). In this case, please download DiVE and place in the visualization folder, and use the following command:

../../dnabarcoder.py visualize -i CBSITS.fasta -c CBSITS.current.classification -rank class -ml 400 -method dive

## Prediction

The third component is to cluster and predict a similarity cut-off for sequence identification based on taxonomic classification. Given a taxonomic level, if higher taxonomic levels are not given, then whole dataset will be used for the prediction.

- To predict a global similarity cut-off at the genus level of the moldITS dataset for example, use the followig command:

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -ranks genus -ml 400

The prediction is saved in the file dnabarcoder/CBSITS.predicted, and the predicted cutoffs are saved in a json format file dnabarcoder/CBSITS.cutoffs.json and a tab delimited format file dnabarcoder/CBSITS.cutoffs.json.txt

- To predict local similarity cut-offs at the genus level of the CBSITS dataset for example, use the followig command:

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -ranks genus -higherranks family,order,class,phylum -ml 400

We can also predict the local cutoffs at a given rank individually:

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -ranks genus -higherranks family -ml 400

or in a given taxa:

../../dnabarcoder.py predict -i CBSITS.fasta -c CBSITS.current.classification -st 0.7 -et 1 -s 0.001 -ranks genus -higherranks phylum -ml 400 -taxa Ascomycota


- To predict global and local similarity cut-offs for the CBSITS dataset at the species level, we first need to remove sequences of species complexes that are indistinguishable by ITS with 100% similarity score:

../../dnabarcoder.py remove -i CBSITS.fasta -c CBSITS.current.classification -ranks species -sim dnabarcoder/CBSITS.sim -ml 400 -t 1

Here t is the threshold or cut-off for removing sequences of the same complex. The results will be saved in dnabarcoder/CBSITS.species.diff.fasta dna dnabarcoder/CBSITS.similar

-To predict global similarity cut-off for species identification of the CBSITS dataset:

../../dnabarcoder.py predict -i dnabarcoder/CBSITS.species.diff.fasta -c CBSITS.current.classification -st 0.9 -et 1 -s 0.001 -ranks species -ml 400 -sim dnabarcoder/CBSITS.sim -prefix CBSITS

The prefix is to save the prediction and the predicted cutoffs in files dnabarcoder/CBSITS.predicted, dnabarcoder/CBSITS.cutoffs.json and dnabarcoder/CBSITS.cutoffs.json.txt.

-To predict local similarity cut-offs for species identification of the CBSITS dataset:

../../dnabarcoder.py predict -i dnabarcoder/CBSITS.species.diff.fasta -c CBSITS.current.classification -st 0.9 -et 1 -s 0.001 -ranks species -higherranks genus,family,order,class,phylum -ml 400 -sim dnabarcoder/CBSITS.sim -prefix CBSITS 

The prediction and cutoffs will be saved in the files dnabarcoder/filamentousfungalITS.predicted, dnabarcoder/CBSITS.cutoffs.json and dnabarcoder/CBSITS.cutoffs.json.txt.


To compute the best cutoffs:

../../dnabarcoder.py best -i dnabarcoder/CBSITS.cutoffs.json -c CBSITS.current.classification

The best similarity cut-offs are saved in json and text format files dnabarcoder/CBSITS.cutoffs.best.json and dnabarcoder/CBSITS.cutoffs.best.txt.


## Classification and Verification

The last component of dnabarcode is to classify a dataset against a reference/barcode dataset using a similarity cut-off or the local cut-offs predicted for the reference dataset.

- To search for the best match of the sequences in the UNITErelease.fasta file, use the following command:

../../dnabarcoder.py classify -i UNITErelease.fasta -r CBSITS.fasta -ml 400

The result is saved in the file dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch

To classify the UNITE sequences based on best matches, using the following commands:

 - Globally, based on only one similarity cut-off:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch -f UNITErelease.fasta -r CBSITS.fasta -c CBSITS.current.classification -cutoff 0.994 -rank species -confidence 0.8334 

Here 0.994 is the global similarity cut-off for sequence identification at the species level. The result will be saved in dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.species.0994.classified. 

- Locally, based on the similarity cut-off predicted for the best match:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch -f UNITErelease.fasta -r CBSITS.fasta -c CBSITS.current.classification -cutoffs dnabarcoder/CBSITS.cutoffs.best.json 

The result will be saved in dnabarcoder/UNITErelease.CBSITS_BLAST.classified. 

Only classify at the species level:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.CBSITS_BLAST.bestmatch -f UNITErelease.fasta -r CBSITS.fasta -c fCBSITS.current.classification -cutoffs dnabarcoder/CBSITS.cutoffs  -rank species

The result will be saved in dnabarcoder/UNITErelease.CBSITS_BLAST.species.classified. 

- To compute classification/assigment accuracy and precision, use the following commands:

../../dnabarcoder.py accuracy -i dnabarcoder/UNITErelease.CBSITS_BLAST.species.classified -c UNITErelease.current.classification -r CBSITS.current.classification

-To visualize the classification/assignment results with Krona:

../../dnabarcoder.py krona -i dnabarcoder/UNITErelease.CBSITS_BLAST.classified -c filamentousfungalITS.current.classification

- To verify the classification results based on phylogenic trees at the species level:

../../dnabarcoder.py verify -i dnabarcoder/UNITErelease.CBSITS_BLAST.classified -c CBSITS.current.classification -r CBSITS.fasta -f UNITErelease.fasta -rank species

## Data

The CBSITS barcode dataset was released in Vu et al. (2019), while the UNITErelease.fasta dataset is the [UNITE general FASTA release](https://plutof.ut.ee/#/doi/10.15156/BIO/786368). The CBSITS.current.classification and UNITErelease.current.classification were updated from [Mycobank](https://www.mycobank.org/).

## Contact person 

Duong Vu (d.vu@wi.knaw.nl)


## References

Vu D. et al. (2019). Large-scale generation and analysis of filamentous fungal DNA barcodes boosts coverage for kingdom fungi and reveals thresholds for fungal species and higher taxon delimitation. Studies in Mycology 92, 135-154.


