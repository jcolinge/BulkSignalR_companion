library(glue)
library(data.table)
library(dplyr)
library(Matrix)
library(SpatialExperiment)
library(STexampleData)

#########################################################################################
# Script to export Visium human DLPFC data 
#
# For more details on this dataset see:
# Maynard and Collado-Torres et al. (2020): https://www.nature.com/articles/s41593-020-00787-0
#########################################################################################

# set up directory
bench.dir <- "/yourDirectory/"

spe <- Visium_humanDLPFC()

############  Export Image  ############################################################

img <- imgRaster(spe, sample_id = "sample_151673", image_id = "lowres")

scaleFactor <- imgData(spe)[imgData(spe)$image_id=="lowres","scaleFactor"]

png(glue("{bench.dir}/Brain_image.export.png"))
plot(img)
dev.off()

############ Prepare Count Matrix ######################################################

colData(spe)$idSpatial <- paste(colData(spe)$array_row,  colData(spe)$array_col,sep = "x")
spot   <- as.data.frame(colData(spe))
coords <- as.data.frame(spatialCoords(spe))

coords$barcode_id <- rownames(spatialCoords(spe))
spot    <- merge(spot, coords, by = "barcode_id",  all = TRUE)
rm(coords)
spot$cell_count[is.na(spot$cell_count)]<-0

# With Spatial Experiment pxl are given in full resolution 
# Using stored scaleFactor for scaling coordinates
spot$pxl_row_in_lowres  <- spot$pxl_row_in_fullres * scaleFactor   
spot$pxl_col_in_lowres <-  spot$pxl_col_in_fullres * scaleFactor

dataframe <-   as.data.frame(as.matrix(assays(spe)[[1]]))
dataframe <- dataframe[,rownames(colData(spe)[spe$in_tissue==1,])] # Keep only InTissue Slot
dataframe$gene_id <- rownames(dataframe)
rownames(dataframe) <- NULL

# We convert matrix count gene_id to symbol using rowData(spe)
dataframe <- inner_join(dataframe, as.data.frame(rowData(spe)[,c( "gene_id" ,  "gene_name")]), by = "gene_id")  
dataframe <- distinct(dataframe, gene_name, .keep_all = TRUE)
rownames(dataframe) <- dataframe$gene_name
dataframe$gene_name <- NULL
dataframe$gene_id   <- NULL

# Match and re-order rows 
ord <- match(colnames(dataframe), colData(spe)[spe$in_tissue==1,]$barcode_id) 
dataframe <- dataframe[,ord]
# Check barcode order
if (!all(colnames(dataframe) == colData(spe)[spe$in_tissue==1,]$barcode_id) ) 
        stop("BarcodID not ordered.", call. = FALSE)
# Change colname with idSpatial 
colnames(dataframe) <- colData(spe)[spe$in_tissue==1,]$idSpatial 

############ Export  labels and counts ####################################################

dataframe <- dataframe[rowSums(dataframe)>0, ] 
print(dim(dataframe))

write.table(spot ,file = glue("{bench.dir}/Brain_label.export.tsv"),
		 col.names=TRUE,row.names =  FALSE, quote = FALSE,sep = "\t")
write.table(dataframe ,file = glue("{bench.dir}/Brain_count.export.tsv"), 
		 col.names=NA,row.names =  TRUE, quote = FALSE,sep = "\t")
