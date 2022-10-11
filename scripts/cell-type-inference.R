library(BulkSignalR)
library(foreach)
library(doMC)
library(igraph)

n.proc <- 4
registerDoMC(n.proc)

# SDC data
data(sdc,package='BulkSignalR')
bsrdm <- prepareDataset(counts = sdc)
bsrdm <- learnParameters(bsrdm)
bsrinf <- initialInference(bsrdm)

# Common TME cell type signatures
data(immune.signatures, package="BulkSignalR")
unique(immune.signatures$signature)
immune.signatures <- immune.signatures[immune.signatures$signature %in% c("B cells","Dentritic cells","Macrophages",
                                                                          "NK cells","T cells","T regulatory cells"),]
data("tme.signatures", package="BulkSignalR")
signatures <- rbind(immune.signatures,tme.signatures[tme.signatures$signature%in%c("Endothelial cells","Fibroblasts"),])
tme.scores <- scoreSignatures(bsrdm, signatures)

# assign cell types to interactions
lr2ct <- assignCellTypesToInteractions(bsrdm, bsrinf, tme.scores)
head(lr2ct)

# cellular network computation and plot
g.table <- cellularNetworkTable(lr2ct)

gCN <- cellularNetwork(g.table)
plot(gCN, edge.width=5*E(gCN)$score)

gSummary <- summarizedCellularNetwork(g.table)
plot(gSummary, edge.width=1+30*E(gSummary)$score)


# relationship with partial EMT------------------------------
# Should be tested HNSCC data instead of SDC!!

# find the ligands
data(p.EMT, package="BulkSignalR")
gs <- p.EMT$gene
triggers <- relateToGeneSet(bsrinf, gs)
triggers <- triggers[triggers$n.genes>1,] # at least 2 target genes in the gs
ligands.in.gs <- intersect(triggers$L, gs)
triggers <- triggers[!(triggers$L%in%ligands.in.gs),]
ligands <- unique(triggers$L)

# link to cell types
cf <- cellTypeFrequency(triggers, lr2ct, min.n.genes=2)
missing <- setdiff(rownames(tme.scores),names(cf$s))
cf$s[missing] <- 0
cf$t[missing] <- 0
op <- par(mar=c(2,10,2,2))
barplot(cf$s,col="lightgray",horiz=T,las=2)
par(op)


# random selections based on random gene sets
qval.thres <- 1e-3
inter <- LRinter(bsrinf)
tg <- tGenes(bsrinf)
tcor <- tgCorr(bsrinf)
good <- inter$qval <= qval.thres
inter <- inter[good,]
tg <- tg[good]
tcor <- tcor[good]
all.targets <- unique(unlist(tg))
r.cf <- list()
for (k in 1:100){ # should 1000 or more
  r.gs <- sample(all.targets,length(intersect(gs,all.targets)))
  r.triggers <- relateToGeneSet(bsrinf, r.gs, qval.thres=qval.thres)
  r.triggers <- r.triggers[r.triggers$n.genes>1,]
  r.ligands.in.gs <- intersect(r.triggers$L, r.gs)
  r.triggers <- r.triggers[!(r.triggers$L%in%r.ligands.in.gs),]
  r <- cellTypeFrequency(r.triggers, lr2ct, min.n.genes=2)
  missing <- setdiff(rownames(tme.scores),names(r$s))
  r$s[missing] <- 0
  r$t[missing] <- 0
  o <- order(names(r$t))
  r$s <- r$s[o]
  r$t <- r$t[o]
  r.cf <- c(r.cf,list(r))
}
r.m.s <- foreach(i=seq_len(length(r.cf)),.combine=rbind) %do% {
  r.cf[[i]]$s
}

# plot results
op <- par(mar=c(2,10,2,2))
boxplot(r.m.s,col="lightgray",horizontal=T,las=2)
pts <- data.frame(x=as.numeric(cf$s[colnames(r.m.s)]),cty=colnames(r.m.s))
stripchart(x~cty, data=pts,add=TRUE,pch=19,col="red")
par(op)
for (cty in rownames(tme.scores))
  cat(cty, ": P=", sum(r.m.s[,cty]>=cf$s[cty])/nrow(r.m.s), "\n", sep="")

