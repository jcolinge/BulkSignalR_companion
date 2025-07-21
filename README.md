# BulkSignalR_companion


## Preamble

We recently integrated `BulkSignalR` to Bioconductor, which required some
small adjustment of the S4 classes design. The main branch of this
companion repository has been adapted to the updated, Bioconductor
version of our library.

In case you need or want to refer to the companion repository as it was
for `BulkSignalR` original version, we created a branch named
**beforeBioconducor** that contains the original, now frozen companion.


## Contents

In this companion repository, we show case applications of the `BulkSignalR`
library to bulk RNA-seq and medium resolution 10x Genomics Visium :tm:
Spatial Transcriptomics (ST):

- SDC-application.R : Use case for Bulk RNA-seq.
- ST-application.R :  Use case for ST using a mouse brain example (see below).

For ST, you need to extract the raw data before. We exemplify this process
for the two data sets analyzed in our manuscript:

- data.extraction.st.Brain.R
- data.extraction.st.BreastCancer.R 

These scripts produce two files needed for ST analysis plus an optional tissue
image file:

- X_count.export.tsv: Raw counts
- X_label.export.tsv: Labeled spots
- X_image.export.png (optional): Tissue image

_Note1:_ R version 4.2 was used to process the data.  
_Note2:_ By default, ST-application.R works with the brain data set.  
To use it with the breast cancer ST data set, you need to edit the
names of the input files.
