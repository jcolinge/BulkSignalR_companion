# BulkSignalR_companion

***Additional material for BulkSignalR package.***

In this repository, you will find example applications of `BulkSignalR` with bulk RNASeq and Spatial Transcriptomics (ST):

- SDC-application.R : Use case for Bulk RNASeq.
- ST-application.R :  Use case for ST using the Brain example (see below).

For ST, you need to extract the raw data before. We examplify this process
for the two datasets analyzed in our manuscript:
- data.extraction.st.Brain.R
- data.extraction.st.BreastCancer.R 

These scripts produce two files needed for ST analysis plus an optional tissue image file.
- X_count.export.tsv: Raw counts
- X_label.export.tsv: Labeled spots
- X_image.export.png (optional): Tissue image

_Note1 :_ R version 4.2 was used to process the data.  
_Note2 :_ By default, ST-application.R works with the brain dataset.  
To work with the Breast Cancer ST dataset, you need to change the input files.
