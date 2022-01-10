# Dnabarcoder

Dnabarcoder is a tool to predict global and local similarity cut-offs for fungal sequence identification. It was implemented in Python which takes DNA barcodes in a fasta file  and their taxonomic classification in a tab delimited file as inputs (see data/filamentousfungalITS.fasta and data/filamentousfungalITS.current.classification for the format of the files). The output of dnabarcode will be saved in an folder given by the user. If this folder is not given, the folder namely dnabarcoder will be created. Dnabarcoder contains four components: analysis, visualization, prediction, and classification to help analyze and predict similarity cut-offs for a dataset of barcodes as well as its subclades, and to classify a dataset against the barcode dataset with the predicted cut-offs. For every function of dnabarcoder, a figure is generated to interpret the result.

Although dnabarcoder was initially developed for fungi, it is applicable to any other organisms using DNA barcodes for identification.

## Dependencies:

- BLAST for DNA sequence comparison
- [LARGEVIS](https://github.com/rugantio/LargeVis-python3), only for DNA sequence visualization

## Analyzation

The analyzation component was to get an overview, and to study the length, the distribution, and the similarity variation of the sequences at different taxonomic levels. 

- To get an overview of the moldITS.fasta dataset:

../../dnabarcoder.py overview -i filamentousfungalITS.fasta -c filamentousfungalITS.current.classification

- To see the distribution of the barcode lengths:

../../dnabarcoder.py length -i filamentousfungalITS.fasta -l 100

Here l is the interval length. Analyzing sequence lengths is important to decide the minimum BLAST alignment length mc. 

-To get the distribution of the sequences at different taxonomic level. In the following example, the distribution of the sequences is computed from the species to the class level:

../../dnabarcoder.py distribution -i filamentousfungalITS.fasta -c filamentousfungalITS.current.classification -p 3,4,5,6,7            

where p is the position of the taxonomic (phylum, class, order, family, genus, species) level in the classification file.

If we want to visualize the distribution of the sequences with Krona, then we can use the following command:

../../dnabarcoder.py  distribution -i filamentousfungalITS.fasta -c filamentousfungalITS.current.classification -p 2,3,4,5,6,7 -method krona

- To get sequence variation with different taxonomic groups:

../../dnabarcoder.py variation -i moldITS.fasta -c moldITS.current.classification -p 3,4,5,6,7  -mc 400

Here the minimum BLAST alignment length mc is set to 400 as 95% of the barcodes have a length of more than 400bp.

- To compute a similarity matrix for the barcodes:

../../dnabarcoder.py sim -i filamentousfungalITS.fasta -mc 400

## Visualization

The second component of dnabarcoder is to visualize the sequences-based 2D/3D “embeddings” using Matplotlib. Sequences’ coordinates are computed using LargeVis.
Together with sequence distribution and variation, visu-alization helped evaluate the predicted similarity cut-offs and classifi-cation results. 

../../dnabarcoder.py visualize -i filamentousfungalITS.fasta -c filamentousfungalITS.current.classification -p 3 -mc 400

Here the sequences are colored by on the taxa at the position 3 (the class level) in the classification file. 

We can also visualize the sequences using [DIVE](https://github.com/NLeSC/DiVE). In this case, please download DiVE and place in the visualization folder, and use the following command:

../../dnabarcoder.py visualize -i filamentousfungalITS.fasta -c filamentousfungalITS.current.classification -p 3 -mc 400 -method dive

## Prediction

The third component is to cluster and predict a similarity cut-off for sequence identification based on taxonomic classification. Given a taxonomic level, if higher taxonomic levels are not given, then whole dataset will be used for the prediction.

- To predict a global similarity cut-off at the genus level of the moldITS dataset for example, use the followig command:

../../dnabarcoder.py predict -i filamentousfungalITS.fasta -c filamentousfungalITS.current.classification -st 0.7 -et 1 -s 0.001 -p 6 -mc 400

- To predict local similarity cut-offs at the genus level of the filamentousfungalITS dataset for example, use the followig command:

../../dnabarcoder.py predict -i filamentousfungalITS.fasta -c filamentousfungalITS.current.classification -st 0.7 -et 1 -s 0.001 -p 6 -hp 5,4,3,2 -mc 400

The prediction will be saved in the json file dnabarcoder/filamentousfungalITS.predicted and the cut-offs will be saved in the json file dnabarcoder/filamentousfungalITS.cutoffs

- To predict global and local similarity cut-offs for the moldITS dataset at the species level, we first need to remove sequences of species complexes that are indistinguishable by ITS with 100% similarity score:

../../dnabarcoder.py remove -i filamentousfungalITS.fasta -c filamentousfungalITS.current.classification -p 7 -sim dnabarcoder/moldITS.sim -mc 400 -t 1

Here t is the threshold or cut-off for removing sequences of the same complex. The results will be saved in dnabarcoder/filamentousfungalITS.diff.fasta dna dnabarcoder/filamentousfungalITS.similar

-To predict global similarity cut-off for species identification of the moldITS dataset:

../../dnabarcoder.py -i filamentousfungalITS.diff.fasta -c filamentousfungalITS.current.classification -st 0.9 -et 1 -s 0.001 -p 7 -mc 400 -sim dnabarcoder/filamentousfungalITS.sim -prefix filamentousfungalITS

-To predict local similarity cut-offs for species identification of the moldITS dataset:

../../dnabarcoder.py -i filamentousfungalITS.diff.fasta -c filamentousfungalITS.current.classification -st 0.9 -et 1 -s 0.001 -p 7 -hp 6,5,4,3,2 -mc 400 -sim dnabarcoder/filamentousfungalITS.sim -prefix moldITS 

The prediction will be saved in the json file dnabarcoder/filamentousfungalITS.predicted and the cut-offs will be saved in the json file dnabarcoder/filamentousfungalITS.cutoffs

## Classification/Assigment

The last component of dnabarcode is to classify a dataset against a reference/barcode dataset using a similarity cut-off or the local cut-offs predicted for the reference dataset.

- To search for the best match of the sequences in the UNITErelease.fasta file, use the following command:

../../dnabarcoder.py classify -i UNITErelease.fasta -r filamentousfungalITS.fasta -mc 400

The result is saved in the file dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.bestmatch

To classify the UNITE sequences based on best matches, using the following commands:

Globally:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.bestmatch -f UNITErelease.fasta -r filamentousfungalITS.fasta -c filamentousfungalITS.current.classification -cutoff 0.994 -rank species -confidence 0.8334 -mc 400

Here 0.994 is the global similarity cut-off for sequence identification at the species level. The result will be saved in dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.species.0994.classified. 

Locally:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.bestmatch -f UNITErelease.fasta -r filamentousfungalITS.fasta -c filamentousfungalITS.current.classification -cutoffs dnabarcoder/filamentousfungalITS.cutoffs -mc 400

The result will be saved in dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.classified. 

Only classify at the species level:

../../dnabarcoder.py classify -i dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.bestmatch -f UNITErelease.fasta -r filamentousfungalITS.fasta -c filamentousfungalITS.current.classification -cutoffs dnabarcoder/filamentousfungalITS.cutoffs -mc 400 -rank species

The result will be saved in dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.species.classified. 

- To compute classification/assigment accuracy and precision, use the following commands:

../../dnabarcoder.py accuracy -i dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.species.classified -c UNITErelease.current.classification -r filamentousfungalITS.current.classification

-To visualize the classification/assignment results with Krona:

../../dnabarcoder.py krona -i dnabarcoder/UNITErelease.filamentousfungalITS_BLAST.classified -c filamentousfungalITS.current.classification

-To verify the classification results based on phylogenic trees at the species level:

../../dnabarcoder.py verify -i dnabarcoder/NITErelease.filamentousfungalITS_BLAST.classified -c filamentousfungalITS.current.classification -r filamentousfungalITS.fasta -f UNITErelease.fasta -rank species

## Data

The filamentousfungalITS barcode dataset was released in Vu et al. (2019), while the UNITErelease.fasta dataset is the [UNITE general FASTA release](https://plutof.ut.ee/#/doi/10.15156/BIO/786368). The filamentousfungalITS.current.classification and UNITErelease.current.classification were updated from [Mycobank](https://www.mycobank.org/).

## Contact person 

Duong Vu (d.vu@wi.knaw.nl)


## References

Vu D. et al. (2019). Large-scale generation and analysis of filamentous fungal DNA barcodes boosts coverage for kingdom fungi and reveals thresholds for fungal species and higher taxon delimitation. Studies in Mycology 92, 135-154.


