# BulkSignalR_companion
  
***Additional material for BulkSignalR package.***

In this repository, you will find :  

Application of BulkSignalR with bulk RNASeq and Spatial Transcriptomics (ST).  

- SDC-application.R : Use case for Bulk RNASEQ.  
- SDC-application-st.R :  Use case for Spatial Transcriptomics (ST).  

For Spatial Transcriptomics (ST), you need to extract the raw data before.  
Two datasets that are in the manuscript are given as examples.    

- data.extraction.st.Brain.R  :
- data.extraction.st.BreastCancer.R : 

They produce three outputs needed for ST analysis.
- X_count.export.png : Raw counts
- X_label.export.tsv : Labeled spots
- X_image.export.png (optional) : Tissue image

They can be used with SDC-application-st.R or
with the code related to the Spatial Transcriptomics section of the
vignette in [BulkSignalR package](https://github.com/jcolinge/BulkSignalR).