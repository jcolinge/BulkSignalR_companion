library(glue)
library(data.table)
library(dplyr)
library(stringr)
library(Matrix)
library(Seurat)
library(SpatialExperiment)

#########################################################################################
# Script to export Visium human breast cancers data 
#
# For more details on this dataset see:
# Wu, S.Z., Al-Eryani, G., Roden, D.L. et al. A single-cell and spatially resolved atlas of human breast cancers. 
# Nat Genet 53, 1334–1347 (2021). https://doi.org/10.1038/s41588-021-00911-1
#
#########################################################################################
download <- TRUE

# set up directory
bench.dir <- "/"
dir <- "BreastVisiumNoH5"
bench.dir <- "yourDirectory"
outdir <- glue("{yourDirectory}/{dir}")

# You might need to do that in bash before.
# export HDF5_USE_FILE_LOCKING='FALSE'

dir.create(outdir, showWarnings = F)
setwd(outdir)

# -----------------------------------------------------------------------------------------------------------
# Download  data
# -----------------------------------------------------------------------------------------------------------
#       . <data-folder>
#       ├── ...
#       ├── <tissue-folder>
#       │   ├── filtered_feature_bc_matrix
#       │   │   ├── barcodes.tsv.gz
#       │   │   ├── features.tsv.gz
#       │   │   └── matrix.mtx.gz 
#       │   ├── spatial
#       │   │   ├── tissue_positions_list.csv
#       └── ...

if (download) {
  # raw feature-barcode matrix
  # note: raw matrix contains all spots; filtered matrix contains only spots over tissue
  url <- "https://zenodo.org/api/files/f3c708c6-3f17-4426-b4f2-74a95f788818/images.pdf"
  fn <- basename(url)

  download.file(url, file.path(outdir, fn))

  #system(glue("tar -C {outdir} -xvzf ", file.path(outdir, fn)))
  #system(glue("rm ", file.path(outdir, fn)))

  url <- "https://zenodo.org/api/files/f3c708c6-3f17-4426-b4f2-74a95f788818/spatial.tar.gz"
  fn <- basename(url)

  download.file(url, file.path(outdir, fn))
  system(glue("tar -C {outdir} -xvzf ", file.path(outdir, fn)))
  system(glue("rm ", file.path(outdir, fn)))

  url <- "https://zenodo.org/api/files/f3c708c6-3f17-4426-b4f2-74a95f788818/metadata.tar.gz"
  fn <- basename(url)

  download.file(url, file.path(outdir, fn))
  system(glue("tar -C {outdir} -xvzf ", file.path(outdir, fn)))
  system(glue("rm ", file.path(outdir, fn)))

  url <- "https://zenodo.org/api/files/f3c708c6-3f17-4426-b4f2-74a95f788818/filtered_count_matrices.tar.gz"
  fn <- basename(url)

  download.file(url, file.path(outdir, fn))
  system(glue("tar -C {outdir} -xvzf ", file.path(outdir, fn)))
  system(glue("rm ", file.path(outdir, fn)))

}

# -----------------------------------------------------------------------------------------------------------
# Extract data
# -----------------------------------------------------------------------------------------------------------

id <- "1160920F"

outdir_matrix_dir <- glue ("{outdir}/filtered_count_matrices/")

# barcodes
file_barcodes <- file.path(outdir_matrix_dir, glue("{id}_filtered_count_matrix"), "barcodes.tsv.gz")
df_barcodes <- read.csv(file_barcodes, sep = "\t", header = FALSE, 
                        col.names = c("barcode_id"))

# features
file_features <- file.path(outdir_matrix_dir, glue("{id}_filtered_count_matrix"), "features.tsv.gz")
df_features <- read.csv(file_features, sep = "\t", header = FALSE, 
                        col.names = c("gene_id", "gene_name", "feature_type"))

# counts
file_counts <- file.path(outdir_matrix_dir, glue("{id}_filtered_count_matrix"), "matrix.mtx.gz")
counts <- readMM(file = file_counts)

stopifnot(nrow(counts) == nrow(df_features))
stopifnot(ncol(counts) == nrow(df_barcodes))

image_spatial_dir <- glue ("{outdir}/spatial")

# spatial coordinates
file_tisspos <- file.path(outdir, glue("spatial/{id}_spatial"), "tissue_positions_list.csv")
tisspos <- read.csv(file_tisspos, header = FALSE, 
                       col.names=c("barcode_id", "in_tissue", "array_row", "array_col", 
                                   "pxl_row_in_fullres", "pxl_col_in_fullres"))
tisspos$idSpatial <- paste(tisspos$array_row, tisspos$array_col,sep = "x")

file_labels <- file.path(outdir, "metadata", glue("{id}_metadata.csv"))
df_labels    <- read.csv(file_labels, header = TRUE, 
                       col.names=c("barcode_id","nCount_RNA", "nFeature_RNA" , "patientid" ,"subtype","Classification"))

# -----------------------------------------------------------------------------------------------------------
# Create Seurat Object
# -----------------------------------------------------------------------------------------------------------
# Load10X_Spatial is a convenient wrapper when you have access h5.file. 
# But you don't need to use a h5 file necesseraly.
# Read10X & Read10X_Image can do the trick

 my_image <- Read10X_Image(
  as.character(dirname(file_tisspos)),
  image.name = "tissue_lowres_image.png",
  filter.matrix = TRUE
)

seurat.object <- CreateSeuratObject(
     counts = Read10X( data.dir = as.character(dirname(file_features)) ,  gene.column = 1),
     assay = 'Spatial'
)

# -----------------------------------------------------------------------------------------------------------
# Export Image
# -----------------------------------------------------------------------------------------------------------
image <- my_image[Cells(x = seurat.object)] 

DefaultAssay(object = image) <- 'Spatial'

seurat.object[['Slice1']] <- image

# Extract Raster Image from Seurat Object
img <- GetImage(seurat.object, mode = c("raster")) 

png(glue("{bench.dir}/spatial/BreastCancer_image.export.png"))
plot(img)
dev.off()

# -----------------------------------------------------------------------------------------------------------
# Prepare Count Matrix
# -----------------------------------------------------------------------------------------------------------

spot <-   as.data.frame(seurat.object@meta.data)
coords <- as.data.frame(GetTissueCoordinates(seurat.object)) # give lowRes

coords$barcode_id <- rownames(coords)
spot$barcode_id  <- rownames(spot)

spot    <- merge(spot, df_labels, by = "barcode_id",  all = TRUE)
spot    <- merge(spot, tisspos, by = "barcode_id",  all = TRUE)
spot    <- merge(spot, coords, by = "barcode_id",  all = TRUE)

# Export count
dataframe <-   as.data.frame(as.matrix(seurat.object@assays$Spatial@counts))

# Keep only InTissue Slot
tisspos <- tisspos[tisspos$in_tissue==1,] 

# match and re-order rows in df_tisspos
ord <- match(tisspos[tisspos$in_tissue==1,]$barcode_id, colnames(dataframe) ) 
dataframe <- dataframe[,ord]

# Check barcode order
if (!all(colnames(dataframe) == tisspos[tisspos$in_tissue==1,]$barcode_id))  
        stop("BarcodID not ordered.", call. = FALSE)
# Change colname with idSpatial
colnames(dataframe) <- tisspos[tisspos$in_tissue==1,]$idSpatial 

# -----------------------------------------------------------------------------------------------------------
# Export Count Matrix
# -----------------------------------------------------------------------------------------------------------

dataframe <- dataframe[rowSums(dataframe)>0, ] 

write.table(spot ,file = glue("{bench.dir}/spatial/BreastCancer_label.export.tsv"), col.names=TRUE,row.names =  FALSE, quote = FALSE,sep = "\t")
write.table(dataframe ,file = glue("{bench.dir}/spatial/BreastCancer_count.export.tsv"), col.names=NA,row.names =  TRUE, quote = FALSE,sep = "\t")

