## An Introduction to Biostatistical Analysis with R
## R codes for lecture on 05/08/2020

#### 0 check and install (if necessary) packages required in this script
required.packages <- c('ape', 'gplots', 'fields', 'cluster')
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages, repos="http://cran.us.r-project.org")
}
# load required packages
require(ape)
require(gplots)
require(fields)
require(cluster)

## Summarize information existing in a number of variables
# Unsupervised and supervised classification

# this data set was analyzed in Zhao 2011 (Nature Communications 2:467)
line <- read.csv("RiceDiversityLine.csv")
pheno <- read.csv("RiceDiversityPheno.csv")
geno <- read.csv("RiceDiversityGeno.csv")
line.pheno <- merge(line, pheno, by.x = "NSFTV.ID", by.y = "NSFTVID")
alldata <- merge(line.pheno, geno, by.x = "NSFTV.ID", by.y = "NSFTVID")
rownames(alldata) <- alldata$NSFTV.ID

# check the name of columns
colnames(alldata)[1:60]

# Unsupervided classificatin: Cluster analysis

#### analysis of marker data
data.mk <- alldata[, 50:ncol(alldata)]
subpop <- alldata$Sub.population
dim(data.mk)

# calculate Euclid distance
d <- dist(data.mk)
head(d)
as.matrix(d)[1:6,1:6]

# cluster samples based on the complete linkage method
tre <- hclust(d)
tre

# draw dendrogram
plot(tre)

# convert to a phylo object defined in the ape package
phy <- as.phylo(tre)

# plot as a phylo object
plot(phy)

# add colors to edges
head(phy$edge)
head(subpop[phy$edge[,2]], 10)
col <- as.numeric(subpop[phy$edge[,2]])
edge.col <- ifelse(is.na(col), "gray", col)

# plot a dendrogram
plot(phy, edge.color = edge.col, show.tip.label = F)

# different types of dendrogram
pdf("fig4.pdf", width = 10, height = 10)
op <- par(mfrow = c(2, 2), mar = rep(0, 4))
plot(phy, edge.color = edge.col, type = "phylogram", show.tip.label = F)
plot(phy, edge.color = edge.col, type = "cladogram", show.tip.label = F)
plot(phy, edge.color = edge.col, type = "fan", show.tip.label = F)
plot(phy, edge.color = edge.col, type = "unrooted", show.tip.label = F)
par(op)
dev.off()

# create an own function
myplot <- function(tre, subpop, type = "unrooted", ...) {
	phy <- as.phylo(tre)
	col <- as.numeric(subpop[phy$edge[,2]])
	edge.col <- ifelse(is.na(col), "gray", col)
	plot(phy, edge.color = edge.col, type = type, show.tip.label = F, ...)
}

# use the function
d <- dist(data.mk)
tre <- hclust(d)
myplot(tre, subpop)
myplot(tre, subpop, type = "cladogram")

# try different methods for calculating distance
pdf("fig5.pdf", width = 10, height = 10)
op <- par(mfrow = c(2, 2), mar = rep(0, 4))
d <- dist(data.mk, method = "euclidean") # default method
myplot(hclust(d), subpop)
d <- dist(data.mk, method = "manhattan")
myplot(hclust(d), subpop)
d <- dist(data.mk, method = "minkowski", p = 1.5)
myplot(hclust(d), subpop)
d <- as.dist(1 - cor(t(data.mk)))
myplot(hclust(d), subpop)
par(op)
dev.off()

# try different clustering methods
pdf("fig6.pdf", width = 10, height = 10)
d <- dist(data.mk)
op <- par(mfrow = c(2, 3), mar = rep(0, 4))
tre <- hclust(d, method = "complete") # default method
myplot(tre, subpop)
tre <- hclust(d, method = "single")
myplot(tre, subpop)
tre <- hclust(d, method = "average")
myplot(tre, subpop)
tre <- hclust(d, method = "median")
myplot(tre, subpop)
tre <- hclust(d, method = "centroid")
myplot(tre, subpop)
tre <- hclust(d, method = "ward.D2")
myplot(tre, subpop)
par(op)
dev.off()

# focus on two clustering methods
op <- par(mfrow = c(1, 2), mar = rep(0, 4))
d <- dist(data.mk)
tre <- hclust(d, method = "complete")
myplot(tre, subpop, type = "phylogram")
tre <- hclust(d, method = "ward.D2")
myplot(tre, subpop, type = "phylogram")
par(op)

##### analysis of phenotypic data ####
# preparation of data
required.traits <- c("Flowering.time.at.Arkansas",
	"Flowering.time.at.Faridpur", "Flowering.time.at.Aberdeen",
	"Culm.habit", "Flag.leaf.length", "Flag.leaf.width",
	"Panicle.number.per.plant", "Plant.height", "Panicle.length",
	"Primary.panicle.branch.number", "Seed.number.per.panicle",
	"Florets.per.panicle", "Panicle.fertility", "Seed.length",
	"Seed.width","Brown.rice.seed.length", "Brown.rice.seed.width",
	"Straighthead.suseptability","Blast.resistance",
	"Amylose.content", "Alkali.spreading.value", "Protein.content")
data.tr <- alldata[, required.traits]
missing <- apply(is.na(data.tr), 1, sum) > 0
data.tr <- data.tr[!missing, ]
subpop.tr <- alldata$Sub.population[!missing]

# scaling
data.tr <- scale(data.tr)

# perform clusetering for both varieties and traits
d <- dist(data.tr)
tre.var <- hclust(d, method = "ward.D2")
d <- dist(t(data.tr))
tre.tra <- hclust(d, "ward.D2")
op <- par(mfrow = c(1, 2))
myplot(tre.var, subpop.tr, type = "phylogram")
plot(tre.tra, cex = 0.5)
par(op)

# perform clustring from both sides
pdf("fig9.pdf")
heatmap(data.tr, margins = c(12,2))
dev.off()

# this part is optional (plot with heatmap.2)
require(gplots)
pdf("fig9-2.pdf", height = 12)
heatmap.2(data.tr, margins = c(12,2), col=redgreen(256), trace = "none", lhei = c(2,10), cexRow = 0.3)
dev.off()

# perform clustering with appropriate methods
pdf("fig10.pdf")
heatmap(data.tr, Rowv = as.dendrogram(tre.var),
				Colv = as.dendrogram(tre.tra),
				RowSideColors = as.character(as.numeric(subpop.tr)),
				labRow = subpop.tr,
				margins = c(12, 2))
dev.off()

# this part is again optional
pdf("fig10-2.pdf", height = 12)
heatmap.2(data.tr, Rowv = as.dendrogram(tre.var),
				Colv = as.dendrogram(tre.tra),
				RowSideColors = as.character(as.numeric(subpop.tr)),
				labRow = subpop.tr,
				margins = c(12,2), col=redgreen(256), trace = "none", lhei = c(2,10), cexRow = 0.3)
dev.off()

# perform clustering based on marker genotypes for determining row order
data.mk2 <- data.mk[!missing, ]
d <- dist(data.mk2)
tre.mrk <- hclust(d, method = "ward.D2")
pdf("fig11.pdf")
heatmap(data.tr, Rowv = as.dendrogram(tre.mrk),
				Colv = as.dendrogram(tre.tra),
				RowSideColors = as.character(as.numeric(subpop.tr)),
				labRow = subpop.tr,
				margins = c(12, 2))
dev.off()

# this part is optional
pdf("fig11-2.pdf", height = 12)
heatmap.2(data.tr, Rowv = as.dendrogram(tre.mrk),
				Colv = as.dendrogram(tre.tra),
				RowSideColors = as.character(as.numeric(subpop.tr)),
				labRow = subpop.tr,
				margins = c(12,2), col=redgreen(256), trace = "none", lhei = c(2,10), cexRow = 0.3)
dev.off()

###### classification based on the hierarchical clustering result #####
# classify samples with the cutree function
d <- dist(data.mk)
tre <- hclust(d, method = "ward.D2")
cluster.id <- cutree(tre, k = 5)
cluster.id

# visualize the result
op <- par(mfrow = c(1,2), mar = rep(0, 4))
myplot(tre, cluster.id, type = "phylogram")
myplot(tre, subpop, type = "phylogram", direction = "leftwards")
par(op)

# obtain the cross table of classification
table(cluster.id, subpop)

# visualize the result in the PCA space
pca <- prcomp(data.mk)
op <- par(mfrow = c(1,2))
plot(pca$x[,1:2], pch = cluster.id, col = as.numeric(subpop))
plot(pca$x[,3:4], pch = cluster.id, col = as.numeric(subpop))
par(op)

# kmeans clustering
kms <- kmeans(data.mk, centers = 5)
kms

# repeat kmeans clustering
for(i in 1:5) {
	kms <- kmeans(data.mk, centers = 5)
	print(table(kms$cluster, subpop))
}

# start from multiple sets of initial points
for(i in 1:5) {
	kms <- kmeans(data.mk, centers = 5, nstart = 50)
	print(table(kms$cluster, subpop))
}

# compare results
table(kms$cluster, subpop)
table(cluster.id, subpop)
table(kms$cluster, cluster.id)

# match id between the results of kmeans and hclust
convert.table <- apply(table(kms$cluster, cluster.id), 1, which.max)
convert.table
cluster.id.kms <- convert.table[kms$cluster]
pdf("fig14.pdf", width = 8, height = 8)
op <- par(mfrow = c(2,2))
plot(pca$x[,1:2], pch = cluster.id, col = as.numeric(subpop), main = "hclust")
plot(pca$x[,3:4], pch = cluster.id, col = as.numeric(subpop), main = "hclust")
plot(pca$x[,1:2], pch = cluster.id.kms, col = as.numeric(subpop), main = "kmeans")
plot(pca$x[,3:4], pch = cluster.id.kms, col = as.numeric(subpop), main = "kmeans")
par(op)
dev.off()

#### decide the appropriate number of groups #####
# visualize within-group and between-groups sum of squares
n <- nrow(data.mk)
wss <- rep(NA, 10)
wss[1] <- (n - 1) * sum(apply(data.mk, 2, var))
for(i in 2:10) {
	print(i)
	res <- kmeans(data.mk, centers = i, nstart = 50)
	wss[i] <- sum(res$withinss)
}
plot(1:10, wss, type = "b", xlab = "Number of groups",
			 ylab = "Within groups sum of squares")

#### find samples which are not clearly assigned to one cluseter ####
# calculate distance to the centers of five groups
require(fields)
d2ctr <- rdist(kms$centers, data.mk)
d2ctr
apply(d2ctr, 2, which.min)
kms$cluster

# prepare own function find the second best
nth.min <- function(x, n) {
	sort(x)[n]
}
nth.min(-10:10, 3)
d.1st <- apply(d2ctr, 2, min)
d.2nd <- apply(d2ctr, 2, nth.min, n = 2)

# evaluate unclearness of classificatinos (shadow values)
shadow <- 2 * d.1st / (d.1st + d.2nd)
unclear <- shadow > 0.9

# visualize the result
cluster.id.kms[unclear] <- 20
op <- par(mfrow = c(1,2))
plot(pca$x[,1:2], pch = cluster.id.kms, col = as.numeric(subpop), main = "kmeans")
plot(pca$x[,3:4], pch = cluster.id.kms, col = as.numeric(subpop), main = "kmeans")
par(op)

#### select representative varieties/lines based on trait data ####
# select 48 varieties/lines (by classifying all samples into 48 groups)
require(cluster)
n.sel <- 48
kmed <- pam(data.tr, k = n.sel)
kmed
kmed$id.med

# look at the distribution of 48 varieties/lines selected as medoids
pca.tr <- prcomp(data.tr, scale = T)
kmed$id.med
pch <- rep(1, nrow(data.tr))
pch[kmed$id.med] <- 19
op <- par(mfrow = c(1,2))
plot(pca.tr$x[,1:2], col = as.numeric(subpop.tr), pch = pch)
plot(pca.tr$x[,3:4], col = as.numeric(subpop.tr), pch = pch)
par(op)

# compare histogram between three datasets (all, 48 selected k-medoids, 48 selected randomly)
op <- par(mfcol = c(2,4))
for(i in 1:4) {
	res <- hist(pca.tr$x[,i], main = paste("PC", i, "all"))
	hist(pca.tr$x[kmed$id.med, i], breaks = res$breaks, main = paste("PC", i, "k-medoids"))
}
par(op)

#### Preparation for homework ####
mydata <- data.frame(
  number = alldata$Panicle.number.per.plant,
  length = alldata$Panicle.length,
  branch = alldata$Primary.panicle.branch.number,
  seed = alldata$Seed.number.per.panicle,
  florets = alldata$Florets.per.panicle
)
missing <- apply(is.na(mydata), 1, sum) > 0
mydata <- mydata[!missing, ]
subpop <- alldata$Sub.population[!missing]

