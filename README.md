# dnabarcoder

Dnabarcoder is a tool to predict global and local similarity cut-offs for fungal sequence identification. It was implemented in Python which takes DNA barcodes in a fasta file  and their taxonomic classification in a tab delimited file as inputs (see data/moldITS.fasta and data/moldITS.current.classification). Dnabarcoder contains four components: analyzation, visualization, prediction, and classification to help analyze and predict similarity cut-offs for a dataset of barcodes as well as its subclades, and to classify a dataset against the barcode dataset with the predicted cut-offs.

Although dnabarcoder was initially developed for fungi, it is applicable to any other organisms using DNA barcodes for identification.

## Dependencies:

- BLAST for DNA sequence comparison
- LARGEVIS (https://github.com/rugantio/LargeVis-python3) for DNA sequence visualization

## analyzation

The analyzation component was to study the length, the distribution, and the similarity variation of the sequences at different taxonomic levels.
## visualization

## prediction

## classification


## Contact person 

Duong Vu (d.vu@wi.knaw.nl)


## References

Vu D. et al. (2019). Massive fungal biodiversity data re-annotation with multi-level clustering. Scientific Reports 4: 6837.


