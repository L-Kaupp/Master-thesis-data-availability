# Master-thesis-data-availability
This repository contains all processing scripts for the snRNAseq analysis of the dataset published by Yaghmaeian Salmani et al.(2024) (DOI: https://doi.org/10.7554/eLife.89482) and the image analysis pipeline for Arivis 4.0.0. 

# For the snRNAseeq analysis of the dataset published by Yaghmaeian Salmani et al.(2024):

To start go on GEO (accession number GSE233866) and download the count matrix "GSE233866_untreated_counts.csv.gz". 
Be aware that the code only works properly when adding an initial "," at the first position of the "GSE233866_untreated_counts.csv.gz" file to grant the correct dataframe format. 

Then run the scripts in the following order:
1. QC_thesis.ipynb
2. Noramlization_thesis.r
3. Figures_thesis.r

# For the image analysis:

Open Arivis 4.0.0 and load in the image of your choice containing a neurobiotin-filled cell visualized with streptavidin - Pacific Blue according to the protocol described in the METHODS section of this thesis. Load in the file "Machine_Learning_Segmentation.pipeline" and run. 
