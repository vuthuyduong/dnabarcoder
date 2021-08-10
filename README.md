# Dnabarcoder

Dnabarcoder is a tool to predict global and local similarity cut-offs for fungal sequence identification. It was implemented in Python which takes DNA barcodes in a fasta file  and their taxonomic classification in a tab delimited file as inputs (see data/moldITS.fasta and data/moldITS.current.classification for the format of the files). The output of dnabarcode will be saved in an folder given by the user. If this folder is not given, the folder namely dnabarcoder will be created. Dnabarcoder contains four components: analyzation, visualization, prediction, and classification to help analyze and predict similarity cut-offs for a dataset of barcodes as well as its subclades, and to classify a dataset against the barcode dataset with the predicted cut-offs. For every function of dnabarcoder, a figure is generated to interpret the result.

Although dnabarcoder was initially developed for fungi, it is applicable to any other organisms using DNA barcodes for identification.

## Dependencies:

- BLAST for DNA sequence comparison
- [LARGEVIS](https://github.com/rugantio/LargeVis-python3) for DNA sequence visualization

## Analyzation

The analyzation component was to get an overview, and to study the length, the distribution, and the similarity variation of the sequences at different taxonomic levels. 

- To get an overview of the moldITS.fasta dataset:

../../dnabarcoder.py overview -i moldITS.fasta -c moldITS.current.classification

- To see the distribution of the barcode lengths:

../../dnabarcoder.py length -i moldITS.fasta -l 100

Here l is the interval length. Analyzing sequence lengths is important to decide the minimum BLAST alignment length mc. 

-To get the distribution of the sequences at different taxonomic level. In the following example, the distribution of the sequences is computed from the species to the class level:

../../dnabarcoder.py distribution -i moldITS.fasta -c moldITS.current.classification -p 3,4,5,6,7            

where p is the position of the taxonomic level in the classification file.

If we want to visualize the distribution of the sequences with Krona, then we can use the following command:

../../dnabarcoder.py  distribution -i moldITS.fasta -c moldITS.current.classification -p 2,3,4,5,6,7 -method krona

- To get sequence variation with different taxonomic groups:

../../dnabarcoder.py variation -i moldITS.fasta -c moldITS.current.classification -p 3,4,5,6,7  -mc 400

Here the minimum BLAST alignment length mc is set to 400 as 95% of the barcodes have a length of more than 400bp.

- To compute a similarity matrix for the barcodes:

../../dnabarcoder.py sim -i moldITS.fasta -mc 400

## Visualization

The second component of dnabarcoder is to visualize the sequences-based 2D/3D “embeddings” using Matplotlib. Sequences’ coordinates are computed using LargeVis.
Together with sequence distribution and variation, visu-alization helped evaluate the predicted similarity cut-offs and classifi-cation results. 

../../dnabarcoder.py visualize -i moldITS.fasta -c moldITS.current.classification -p 3 -mc 400

Here the sequences are colored by on the taxa at the position 3 (the class level) in the classification file. 

We can also visualize the sequences using [DIVE](https://github.com/NLeSC/DiVE). In this case, please download DiVE and place in the visualization folder, and use the following command:

../../dnabarcoder.py visualize -i moldITS.fasta -c moldITS.current.classification -p 3 -mc 400 -method dive

## Prediction

The third component is to cluster and predict a similarity cut-off for sequence identification based on taxonomic classification. Given a taxonomic level, if higher taxonomic levels are not given, then whole dataset will be used for the prediction.

- To predict a global similarity cut-off at the genus level of the moldITS dataset for example, use the followig command:

../../dnabarcoder.py predict -i moldITS.fasta -c moldITS.current.classification -st 0.7 -et 1 -s 0.001 -p 6 -mc 400

- To predict local similarity cut-offs at the genus level of the moldITS dataset for example, use the followig command:

../../dnabarcoder.py predict -i moldITS.fasta -c moldITS.current.classification -st 0.7 -et 1 -s 0.001 -p 6 -hp 5,4,3,2 -mc 400

The prediction will be saved in the file dnabarcoder/moldITS.predicted and the cut-offs will be saved in the file dnabarcoder/moldITS.cutoffs

- To predict global and local similarity cut-offs for the moldITS dataset at the species level, we first need to remove sequences of species complexes that are indistinguishable by ITS with 100% similarity score:

../../dnabarcoder.py remove -i moldITS.fasta -c moldITS.current.classification -p 7 -sim dnabarcoder/moldITS.sim -mc 400 -t 1

Here t is the threshold or cut-off for removing sequences of the same complex. The results will be saved in dnabarcoder/moldITS.diff.fasta dna dnabarcoder/moldITS.similar

-To predict global similarity cut-off for species identification of the moldITS dataset:

../../dnabarcoder.py -i moldITS.diff.fasta -c moldITS.current.classification -st 0.9 -et 1 -s 0.001 -p 7 -mc 400 -sim dnabarcoder/moldITS.sim -prefix moldITS 

-To predict local similarity cut-offs for species identification of the moldITS dataset:

../../dnabarcoder.py -i moldITS.diff.fasta -c moldITS.current.classification -st 0.9 -et 1 -s 0.001 -p 7 -hp 6,5,4,3,2 -mc 400 -sim dnabarcoder/moldITS.sim -prefix moldITS 

The prediction will be saved in the file dnabarcoder/moldITS.predicted and the cut-offs will be saved in the file dnabarcoder/moldITS.cutoffs

## Classification/Assigment

The last component of dnabarcode is to classify a dataset against a reference/barcode dataset using a similarity cut-off or the local cut-offs predicted for the reference dataset.

- To classify the SH.fasta file, use the following command:

../../dnabarcoder.py classify -i SH.fasta -r moldITS.fasta -c moldITS.classification -mc 400

- Using a global similarity cutoff, use the following command:

../../dnabarcoder.py classify -i SH.fasta -r moldITS.fasta -c moldITS.classification -cutoff 0.97 -mc 400

- Using lobal similarity cutoffs, please use the following command:

../../dnabarcoder.py classify -i SH.fasta -r moldITS.fasta -c moldITS.classification -cutoffs dnabarcoder/moldITS.cutoffs -mc 400

- If there exists a file of classified sequences by BLAST or other tools, we can assign the classified sequences with the predicted cut-offs as follows:

../../dnabarcoder.py assign -i dnabarcoder/SH.moldITS_BLAST.classified -f SH.fasta -r moldITS.fasta -c moldITS.current.classification -cutoffs dnabarcoder/moldITS.cutoffs -mc 400

- We can also assign the sequences at a specific taxonomic level. For example, at the species level:

Globally:

../../dnabarcoder.py assign -i dnabarcoder/SH.moldITS_BLAST.classified -f SH.fasta -r moldITS.fasta -c moldITS.current.classification -cutoff 0.994 -rank species -confidence 0.8334 -mc 400

Locally:

../../dnabarcoder.py assign -i dnabarcoder/SH.moldITS_BLAST.classified -f SH.fasta -r moldITS.fasta -c moldITS.current.classification -cutoffs dnabarcoder/moldITS.cutoffs -rank species -mc 400

The result will be saved in dnabarcoder/SH.moldITS_BLAST.species.assigned. 

- To compute classification/assigment accuracy and precision, use the following commands:

../../dnabarcoder.py accuracy -i dnabarcoder/SH.moldITS_BLAST.species.assigned -c SH.current.classification -r moldITS.current.classification

-To visualize the classification/assignment results with Krona:

../../dnabarcoder.py krona -i dnabarcoder/SH.moldITS_BLAST.classified -c moldITS.current.classification

../../dnabarcoder.py krona -i dnabarcoder/SH.moldITS_BLAST.assigned -c moldITS.current.classification

## Data

The mold ITS barcode dataset was released in Vu et al. (2019), while the SH.fasta dataset is the [UNITE general FASTA release](https://plutof.ut.ee/#/doi/10.15156/BIO/786368). The moldITS.current.classification and SH.current.classification were updated from [Mycobank](https://www.mycobank.org/).

## Contact person 

Duong Vu (d.vu@wi.knaw.nl)


## References

Vu D. et al. (2019). Large-scale generation and analysis of filamentous fungal DNA barcodes boosts coverage for kingdom fungi and reveals thresholds for fungal species and higher taxon delimitation. Studies in Mycology 92, 135-154.


