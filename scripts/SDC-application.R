library(BulkSignalR)
library(igraph)

# activate parallel computing for faster model training [optional]
library(doParallel)
n.proc <- 6
cl <- makeCluster(n.proc)
registerDoParallel(cl)

# prepare data and train the statistical model
data(sdc,package="BulkSignalR")
normal <- grep("^N", names(sdc))
bsrdm <- BSRDataModel(sdc[, -normal])
bsrdm <- learnParameters(bsrdm, quick=FALSE, plot.folder=".", verbose=TRUE)
bsrdm

# score ligand-receptor interactions
bsrinf <- BSRInference(bsrdm)
bsrinf

# save(bsrdm,file="bsrdm.rda")
# save(bsrinf,file="bsrinf.rda")
load("bsrdm.rda")
load("bsrinf.rda")

# basic pathway statistics
pwstat <- getPathwayStats(bsrinf, qval.thres=0.01)
head(pwstat)

# examples of reducing the LR analysis --------------------------

# best pathway per LR pair
head(LRinter(bsrinf),n=10)
bsrinf.red <- reduceToBestPathway(bsrinf)
head(LRinter(bsrinf.red))

# reduction to the ligand
bsrinf.redL <- reduceToLigand(bsrinf)
head(LRinter(bsrinf.redL))

# reduction to the receptor
bsrinf.redR <- reduceToReceptor(bsrinf)
head(LRinter(bsrinf.redR))

# reduction to pathway
bsrinf.redP <- reduceToPathway(bsrinf)
head(LRinter(bsrinf.redP))

# reduction to pathway followed by reduction to best pathway
# to avoid redundancy
bsrinf.redPP <- reduceToBestPathway(bsrinf.redP)
head(LRinter(bsrinf.redPP))

# gene signatures --------------------------------------------

# extract gene signatures to report combined ligand-receptor and
# receptor downstream pathway scores
sum(LRinter(bsrinf.redPP)$qval < 0.01) # number of significant interactions
sum(LRinter(bsrinf.red)$qval < 1e-8)
bsrsig.red <- BSRSignature(bsrinf.red, qval.thres=1e-8)
scores.red <- scoreLRGeneSignatures(bsrdm, bsrsig.red)

# on the screen
simpleHeatmap(scores.red, pointsize=4)

# in a PDF file
pdf(file="SDC-LR-heatmap.pdf", width=6, height=12, pointsize=4,
    useDingbats=FALSE)
simpleHeatmap(scores.red, pointsize=4)
dev.off()

# correlate with the immune microenvironment
data(immune.signatures, package="BulkSignalR")
imm.scores <- scoreSignatures(bsrdm, immune.signatures)

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
write.graph(gLRintra, file="SDC-LR-intracellular-network.graphml",
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
top <- unique(pairs[pairs$pval<1e-20,c("pw.id","pw.name")])
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

# re-scoring ---------------------------------------------------------

# obtain inference significance considering less deep receptor
# downstream signaling

bsrinf.less <- rescoreInference(bsrinf, param=parameters(bsrdm), rank.p=0.75)
head(LRinter(bsrinf), n=10)
head(LRinter(bsrinf.less), n=10)
plot(x=LRinter(bsrinf)$qval, y=LRinter(bsrinf.less)$qval, log="xy", pch=20)
abline(a=0, b=1, col="orange")

# using the full network instead of restricting to observed genes ------

inter <- LRinter(bsrinf)

bsrinf.full <- BSRInference(bsrdm, use.full.network=TRUE) # lengthy computation
inter.full <- LRinter(bsrinf.full)

# end ---------------------------------------------------------------

# stop cluster if parallel computation was used [optional]
stopCluster(cl)

