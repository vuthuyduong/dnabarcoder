# dnabarcoder

dnabarcoder is a tool to predict global and local similarity cut-offs for fungal sequence identification 

## analyzation

[Windows](https://github.com/FastMLC/fMLC/tree/master/Windows)

[Linux](https://github.com/FastMLC/fMLC/tree/master/Linux)

## visualization
There are two datasets available as inputs for fMLC. The "small" dataset contains ~4000 ITS yeast sequences, checked and validated by the specialists at the Westerdijk Fungal Biodiversity Institute. This dataset were analyzed and released in [Vu D. et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5192050/). The "large" dataset contains ~350K ITS fungal sequences downloaded from GenBank (https://www.ncbi.nlm.nih.gov/) which was used in [Vu D. et al. 2014](https://www.nature.com/articles/srep06837) to evaluate the speed of MLC.

<!---[Download](http://www.westerdijkinstitute.nl/Download/SmallDatasetOf4KYeastITSSequences.zip) the small demo dataset.  -->

[Download](http://www.westerdijkinstitute.nl/Download/LargeDatasetOf350KITSSequences.zip) the large demo dataset. 

## prediction

After clustering the DNA sequences by fMLC, the groupings of the sequences can be saved as output of fMLC. A sparse (or complete) similarity matrix (in .sim format) can be saved in the folder where the dataset is given, to capture the similarity structure of the sequences. Based on this similarity matrix, the coordiates of the sequences can be computed and saved (in .outLargeVis format) using LargeVis. Finally, a json file containing the coordinates and metadata of the sequences is resided in the folder DiVE/data folder as an input of DiVE to visualize the data. This json file can be used for visualization by external applications as well.The clustering and visualization results of the two datasets can be found at https://github.com/FastMLC/fMLC/tree/master/data.

## classification

## Contact person 

Duong Vu (d.vu@westerdijkinstitute.nl)


## References

Vu D. et al. (2019). Massive fungal biodiversity data re-annotation with multi-level clustering. Scientific Reports 4: 6837.


