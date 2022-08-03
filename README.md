# BulkSignalR_companion

***Additional material for BulkSignalR package.***

In this repository, you will find :  

Application of `BulkSignalR` with bulk RNASeq and Spatial Transcriptomics (ST).  

- SDC-application.R : Use case for Bulk RNASeq.  
- SDC-application-st.R :  Use case for ST using the Brain example.  

For ST, you need to extract the raw data before.  
Two datasets described in the manuscript can be retrieved as follows.    
- data.extraction.st.Brain.R  :
- data.extraction.st.BreastCancer.R : 

They produce three files needed for ST analysis.
- X_count.export.tsv : Raw counts
- X_label.export.tsv : Labeled spots
- X_image.export.png (optional) : Tissue image

These files are needed to use SDC-application-st.R or
with the code related to the ST section of the
vignette providen inside [BulkSignalR package](https://github.com/jcolinge/BulkSignalR).  

_Note1 :_ R version 4.2 was used to process the data.  
_Note2 :_ By default, SDC-application-st.R works with the brain dataset.  
To work with the Breast Cancer ST dataset, you need to change the inputs accordingly.  