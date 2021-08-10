# dnabarcoder

Dnabarcoder is a tool to predict global and local similarity cut-offs for fungal sequence identification. It was implemented in Python which takes DNA barcodes in a fasta file  and their taxonomic classification in a tab delimited file as inputs (see data/moldITS.fasta and data/moldITS.current.classification for the format of the files). The output of dnabarcode will be saved in an folder given by the user. If this folder is not given, the folder namely dnabarcoder will be created. Dnabarcoder contains four components: analyzation, visualization, prediction, and classification to help analyze and predict similarity cut-offs for a dataset of barcodes as well as its subclades, and to classify a dataset against the barcode dataset with the predicted cut-offs. For every function of dnabarcoder, a figure is generated to interpret the result.

Although dnabarcoder was initially developed for fungi, it is applicable to any other organisms using DNA barcodes for identification.

## Dependencies:

- BLAST for DNA sequence comparison
- LARGEVIS (https://github.com/rugantio/LargeVis-python3) for DNA sequence visualization

## analyzation

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

## visualization

The second component of dnabarcoder is to visualize the sequences-based 2D/3D “embeddings” using Matplotlib. Sequences’ coordinates are computed using LargeVis.
Together with sequence distribution and variation, visu-alization helped evaluate the predicted similarity cut-offs and classifi-cation results. 

../../dnabarcoder.py visualize -i moldITS.fasta -c moldITS.current.classification -p 3 -mc 400

Here the sequences are colored by on the taxa at the position 3 (the class level) in the classification file. 

We can also visualize the sequences using DiVE (https://github.com/NLeSC/DiVE). In this case, please download DiVE and place in the visualization folder, and use the following command:

../../dnabarcoder.py visualize -i moldITS.fasta -c moldITS.current.classification -p 3 -mc 400 -method dive

## prediction

The third component is to cluster and predict a similarity cut-off for sequence identification based on taxonomic classification. Given a taxonomic level, if higher taxonomic levels are not given, then whole dataset will be used for the prediction.

- To predict a global similarity cut-off at the genus level of the moldITS dataset for example, use the followig command:

../../prediction/predict.py -i moldITS.fasta -c moldITS.current.classification -st 0.7 -et 1 -s 0.001 -p 6

- To predict local similarity cut-offs at the genus level of the moldITS dataset for example, use the followig command:

../../prediction/predict.py -i moldITS.fasta -c moldITS.current.classification -st 0.7 -et 1 -s 0.001 -p 6 -hp 5,4,3,2

The prediction will be saved in the file dnabarcoder/moldITS.predicted and the cut-offs will be saved in the file dnabarcoder/moldITS.cutoffs

## classification


## Contact person 

Duong Vu (d.vu@wi.knaw.nl)


## References

Vu D. et al. (2019). Massive fungal biodiversity data re-annotation with multi-level clustering. Scientific Reports 4: 6837.


