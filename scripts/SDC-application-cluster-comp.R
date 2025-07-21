# In this example script, we illustrate the use of BulkSignalR based on the
# differential analysis of gene or protein expression. We use SDC data and
# apply elementary cellular deconvolution to define a group or cluster
# or tumor samples that are immune-rich and another one that contains
# immune-poor samples.
#
# This is just an example. The clustering and the differential gene expression
# analyses are the responsibility of the user, and would most likely rely
# on classical, existing R packages or statistical tests.


library(BulkSignalR)
library(igraph)

# activate parallel computing for faster model training [optional]
library(doParallel)
n.proc <- 2
cl <- makeCluster(n.proc)
registerDoParallel(cl)

# prepare data
data(sdc,package="BulkSignalR")
normal <- grep("^N", names(sdc))
bsrdm <- BSRDataModel(sdc[, -normal])

# ===============================================================
# Sample clustering
# ===============================================================

# define clusters of SDC samples based on their contents in immune cell types
data(immune.signatures, package="BulkSignalR")
unique(immune.signatures$signature)
imm.scores <- scoreSignatures(bsrdm, immune.signatures)
simpleHeatmap(imm.scores)
d <- dist(t(imm.scores))
h <- hclust(d)
plot(h)
clusters <- cutree(h, 4) # 4 because of SDC19 that is an outlier
clusters

# identify the high and low immune content clusters
imm.contents <- foreach(c=1:3, .combine=c) %do% {
  mean(imm.scores[,which(clusters == c)])
}
low.cluster <- which.min(imm.contents)
high.cluster <- which.max(imm.contents)

# find the corresponding columns in the count matrix
colA <- as.integer(which(clusters == high.cluster))
colB <- as.integer(which(clusters == low.cluster))

# ============================================================================
# Differential gene expression (DGE) analysis
# ============================================================================

# This part can be substituted with any piece of code that is favored by
# the user. Here we provide an example using edgeR, but any other such
# library, e.g., DESeq2, would be equally fine).

# Note that we apply edgeR to the original count matrix (sdc) to fit its
# its statistical model requirements. Remember that data were normalized
# differently creating the object bsrdm. Also, some rows were removed since
# they were plenty of zeros. Thos row are also removed from the analyisi
# hereafter to maintain full compatibility with bsrdm contents.

library(edgeR)
conditions <- rep("medium", ncol(sdc))
conditions[colA] <- "high"
conditions[colB] <- "low"
cl <- factor(conditions)
design <- model.matrix(~0+cl)
colnames(design) <- gsub("^cl", "", colnames(design))
cm <- makeContrasts(high-low, levels=design)
dge <- DGEList(sdc[rownames(ncounts(bsrdm)),], genes=rownames(ncounts(bsrdm)))
rownames(dge$counts) <- rownames(ncounts(bsrdm))
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge, design, robust=T)
fit.dge <- glmFit(dge, design)
comparison <- "high - low"
lrt <- glmLRT(fit.dge, contrast=cm[,comparison])
sel.r <- topTags(lrt, n=nrow(dge$counts))
edger.stats <- sel.r$table

# extract the required data from edgeR analysis and the normalized counts
stats <- data.frame(pval=edger.stats$PValue,
                    logFC=edger.stats$logFC,
                    expr=rowMeans(ncounts(bsrdm)[,colA]))
rownames(stats) <- edger.stats$genes

# ===========================================================================


# create the BSRDataModelComp and BSRClusterComp objects for this comparison
bsrdm.comp <- as(bsrdm, "BSRDataModelComp")
bsrcc <- BSRClusterComp(bsrdm.comp, colA, colB, stats)
bsrcc
bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "high.versus.low.immune")
bsrdm.comp

# here, we illustrate adding another cluster comparison (we use the same
# one with a different anme for the sake of the example) and removing it
bsrdm.comp <- addClusterComp(bsrdm.comp, bsrcc, "again")
bsrdm.comp
bsrdm.comp <- removeClusterComp(bsrdm.comp, "again")
bsrdm.comp

# score ligand-receptor interactions
bsrinf.comp <- BSRInferenceComp(bsrdm.comp, "high.versus.low.immune")
bsrinf.comp

# save(bsrdm.comp,file="bsrdm-comp.rda")
# save(bsrinf.comp,file="bsrinf-comp.rda")
load("bsrdm-comp.rda")
load("bsrinf-comp.rda")

# basic pathway statistics
pwstat <- getPathwayStats(bsrinf.comp, qval.thres=0.01)
head(pwstat)

# examples of reducing the LR analysis --------------------------

# best pathway per LR pair
head(LRinter(bsrinf.comp),n=10)
bsrinf.red <- reduceToBestPathway(bsrinf.comp)
head(LRinter(bsrinf.red))

# reduction to the ligand
bsrinf.redL <- reduceToLigand(bsrinf.comp)
head(LRinter(bsrinf.redL))

# reduction to the receptor
bsrinf.redR <- reduceToReceptor(bsrinf.comp)
head(LRinter(bsrinf.redR))

# reduction to pathway
bsrinf.redP <- reduceToPathway(bsrinf.comp)
head(LRinter(bsrinf.redP))

# reduction to pathway followed by reduction to best pathway
# to avoid redundancy
bsrinf.redPP <- reduceToBestPathway(bsrinf.redP)
head(LRinter(bsrinf.redPP))

# gene signatures --------------------------------------------

# extract gene signatures to report combined ligand-receptor and
# receptor downstream pathway scores
sum(LRinter(bsrinf.redPP)$qval < 0.01) # number of significant interactions
sum(LRinter(bsrinf.red)$qval < 1e-6)
bsrsig.red <- BSRSignatureComp(bsrinf.red, qval.thres=1e-6)
bsrsig.red
scores.red <- scoreLRGeneSignatures(bsrdm.comp, bsrsig.red, name.by.pathway=TRUE, rownames.LRP=TRUE)

# on the screen
simpleHeatmap(scores.red, width=6, height=8, pointsize=4)

# in a PDF file
pdf(file="SDC-diff-LR-heatmap.pdf", width=6, height=8, pointsize=4,
    useDingbats=FALSE)
simpleHeatmap(scores.red, width=6, height=12, pointsize=4)
dev.off()


# networks ----------------------------------------------------------

# generate a ligand-receptor network and export it in .graphML
# for Cytoscape or similar tools
gLR <- getLRNetwork(bsrinf.red, qval.thres=1e-8)
gLR
write_graph(gLR,file="SDC-LR-network.graphml",format="graphml")

# play around with igraph functions as an alternative to Cytoscape
plot(gLR)
plot(gLR,
     edge.arrow.size=0.6,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)
# community detection
u.gLR <- as_undirected(gLR) # most algorithms work for undirected graphs only
comm <- cluster_edge_betweenness(u.gLR)
plot(comm,u.gLR,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)
# cohesive blocks
cb <- cohesive_blocks(u.gLR)
plot(cb,u.gLR,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75,
     edge.color="black")

# generate a ligand-receptor network complemented with intra-cellular,
# receptor downstream pathways [computations are a bit longer here]
gLRintra <- getLRIntracellNetwork(bsrinf.red, qval.thres=1e-8)
write_graph(gLRintra, file="SDC-LR-intracellular-network.graphml",
            format="graphml")
lay <- layout_with_kk(gLRintra)
plot(gLRintra,
     layout=lay,
     edge.arrow.size=0.6,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)

# reduce complexity by focusing on strongly targeted pathways
pairs <- LRinter(bsrinf.red)
top <- unique(pairs[pairs$pval<1e-10,c("pw.id","pw.name")])
top
gLRintra.res <- getLRIntracellNetwork(bsrinf.red, qval.thres=0.01,
                                      restrict.pw=top$pw.id)
lay <- layout_with_fr(gLRintra.res)
plot(gLRintra.res,
     layout=lay,
     edge.arrow.size=0.6,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)

# select target genes based on their regulation P-values instead of
# correlations with the receptor
gLRintra.res <- getLRIntracellNetwork(bsrinf.red, qval.thres=0.01,
                                      restrict.pw=top$pw.id, max.pval=1e-3)
lay <- layout_with_fr(gLRintra.res)
plot(gLRintra.res,
     layout=lay,
     edge.arrow.size=0.6,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.75)

# re-scoring ---------------------------------------------------------

# obtain inference significance considering less deep receptor
# downstream signaling

bsrinf.less <- rescoreInference(bsrinf.comp, param=param(bsrdm.comp), rank.p=0.75)
head(LRinter(bsrinf.comp), n=10)
head(LRinter(bsrinf.less), n=10)
plot(x=LRinter(bsrinf.comp)$qval, y=LRinter(bsrinf.less)$qval, log="xy", pch=20)
abline(a=0, b=1, col="orange")

# end ---------------------------------------------------------------

# stop cluster if parallel computation was used [optional]
stopCluster(cl)
