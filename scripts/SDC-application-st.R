library(BulkSignalR)
library(foreach)
library(doMC)
library(glue)

n.proc <- 2
registerDoMC(n.proc)


# set up working directory
bench.dir <- "./"

# load data =================================================

counts <- read.csv(glue("{bench.dir}/Brain_count.export.tsv"), sep="\t", stringsAsFactors=FALSE)
areas <- read.csv(glue("{bench.dir}/Brain_label.export.tsv"), sep="\t", stringsAsFactors=FALSE)

areas$idSpatial <- paste0("X", areas$idSpatial)
areas           <- areas[areas$idSpatial %in% names(counts),]
areas$label <- areas$ground_truth
areas$label[is.na(areas$label)] <- "unclassified"
table(areas$label)


# prepare data =================================================

bsrdm <- prepareDataset(counts,symbol.col=1, 
                              min.count=1, 
                              method="TC")

bsrdm <- learnParameters(bsrdm, quick=FALSE,
                         plot.folder="./",
                         min.positive=2, 
                         verbose=TRUE)
bsrdm
# save(bsrdm,file="spatial2-bsrdm.rda")
load("spatial2-bsrdm.rda")

bsrinf <- initialInference(bsrdm, min.cor=-1)
bsrinf
# save(bsrinf,file="spatial2-bsrinf.rda")
load("spatial2-bsrinf.rda")


# spatial analysis ============================================

bsrinf.red <- reduceToBestPathway(bsrinf)
pairs.red <- LRinter(bsrinf.red)
thres <- 0.01
sum(pairs.red$qval<thres)
s.red <- getLRGeneSignatures(bsrinf.red, qval.thres=thres)
scores.red <- scoreLRGeneSignatures(bsrdm,s.red)


# plot one specific interaction
rownames(scores.red) <- gsub("->","/",rownames(scores.red))
inter <- gsub("\\}","",gsub("\\{","",rownames(scores.red)[20]))
spatialPlot(scores.red[20,], areas, inter, ref.plot=TRUE, dot.size=1)

# generate visual index
spatialIndexPlot(scores.red, areas, "bigplot-BSR.pdf")

# statistical association with tissue areas
assoc.bsr <- spatialAssociation(scores.red, areas)
spatialAssociationPlot(assoc.bsr)

# 2D-projection of score spatial distribution patterns
spatialDiversityPlot(scores.red, assoc.bsr, with.names=TRUE)