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
                              prop=0.01,
                              method="TC")

bsrdm <- learnParameters(bsrdm, quick=FALSE,
                         plot.folder="./",
                         min.positive=2, 
                         verbose=TRUE)
bsrdm
save(bsrdm,file="spatial2-bsrdm.rda")
load("spatial2-bsrdm.rda")

bsrinf <- initialInference(bsrdm, min.cor=-1)
bsrinf
save(bsrinf,file="spatial2-bsrinf.rda")
load("spatial2-bsrinf.rda")


# spatial analysis ============================================

bsrinf.red <- reduceToBestPathway(bsrinf)
pairs.red <- LRinter(bsrinf.red)
thres <- 0.01
min.corr <- 0.01
pairs.red <- pairs.red[pairs.red$qval < thres & pairs.red$LR.corr > min.corr,]
s.red <- getLRGeneSignatures(bsrinf.red, qval.thres=thres)
scores.red <- scoreLRGeneSignatures(bsrdm,s.red)


# plot one specific interaction
inter <- "{CALM1} / {GRM5}" # we have to follow the syntax with {} to be compatible with reduction operations
spatialPlot(scores.red[inter,], areas, inter, ref.plot=TRUE, dot.size=1)

# dissect one interaction
separatedLRPlot(scores.red, "CALM1", "GRM5", ncounts(bsrdm), areas)

# generate visual index
spatialIndexPlot(scores.red, areas, "bench2-plots/bigplot-BSR.pdf")

# statistical association with tissue areas based on a statistical test (Kruskal-Wallis by default)
assoc.bsr <- spatialAssociation(scores.red, areas)
spatialAssociationPlot(assoc.bsr)

# statistical association with tissue areas based on correlations
assoc.bsr.corr <- spatialAssociation(scores.red, areas, test="Spearman")
spatialAssociationPlot(assoc.bsr.corr)

# 2D-projection of score spatial distribution patterns (compatible with all associations)
spatialDiversityPlot(scores.red, assoc.bsr.corr, with.names=TRUE)
spatialDiversityPlot(scores.red, assoc.bsr.corr, with.names=TRUE, score.based=TRUE)
spatialDiversityPlot(scores.red, assoc.bsr.corr, with.names=TRUE, proj="tSNE", score.based=TRUE)

# use of spatial plot for other quantities
mtc <- grep("^MT",rownames(ncounts(bsrdm)))
csc <- colSums(ncounts(bsrdm)[mtc,])/colSums(ncounts(bsrdm))
spatialPlot(csc, areas, "MT contents", ref.plot=TRUE, dot.size=1)

# use of associations for genes instead of interactions
sel.genes <- unique(c(pairs.red[pairs.red$pw.name=="Protein-protein interactions at synapses", "L"],
                      pairs.red[pairs.red$pw.name=="Protein-protein interactions at synapses", "R"],
                      pairs.red$L[grep("^COL", pairs.red$L)]))
assoc.g <- spatialAssociation(ncounts(bsrdm)[sel.genes,], areas, test="Spearman")
spatialAssociationPlot(assoc.g)
