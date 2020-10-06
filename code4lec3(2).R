#### An Introduction to Biostatistical Analysis with R ####
## 3rd Lecture: Summarize information existing in a number of variables
## R codes for lecture on 05/01/2020


#### 0. Required packages ####
# required packages for this code
required.packages <- c("plotly")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages, repos="http://cran.us.r-project.org")
}
# load required packages
require(plotly)


#### 1. Introduction to principal component analysis (PCA) ####
## this data set was analyzed in Zhao 2011 (Nature Communications 2:467)
line <- read.csv("RiceDiversityLine.csv")
pheno <- read.csv("RiceDiversityPheno.csv")
geno <- read.csv("RiceDiversityGeno.csv")
line.pheno <- merge(line, pheno, by.x = "NSFTV.ID", by.y = "NSFTVID")
alldata <- merge(line.pheno, geno, by.x = "NSFTV.ID", by.y = "NSFTVID")

## extract panicle length and flag leaf length
mydata <- data.frame(
	panicle.length  = alldata$Panicle.length,
	leaf.length = alldata$Flag.leaf.length
	)
mydata
missing <- apply(is.na(mydata), 1, sum) > 0
mydata <- mydata[!missing, ]
mydata

## look at the relationship between two variables
plot(mydata)
lim <- range(mydata)
plot(mydata, xlim = lim, ylim = lim)

## statistics for measuring the relationship between two variables
cov(mydata)
cor(mydata)

## subtract the mean from each column to shift the center of data to the origin  
mydata <- sweep(mydata, 2, apply(mydata, 2, mean))
summary(mydata)
cov(mydata)
cor(mydata)
lim <- range(mydata)
plot(mydata, xlim = lim, ylim = lim)
abline(h = 0, v = 0)

## perform principal component analysis (PCA) and draw a scatterplot
res <- prcomp(mydata)
lim <- range(res$x)
plot(res$x, xlim = lim, ylim = lim)
abline(h = 0, v = 0)

## show graphs side by side
op <- par(mfrow = c(1,2))
lim <- range(mydata)
plot(mydata, xlim = lim, ylim = lim)
abline(h = 0, v = 0)
lim <- range(res$x)
plot(res$x, xlim = lim, ylim = lim)
abline(h = 0, v = 0)
par(op)

## show the result of PCA
summary(res)

## look into the details of the result
res
res$sdev
res$rotation

## draw graphs for PCA
op <- par(mfrow = c(1,2))
plot(res)
biplot(res)
par(op)

## plot again
lim <- range(mydata)
plot(mydata, xlim = lim, ylim = lim)
abline(h = 0, v = 0)

## arbitrary line
u.temp <- c(1 / sqrt(2), 1 / sqrt(2))
abline(0, u.temp[2] / u.temp[1], col = "red")

## draw scores
score.temp <- as.matrix(mydata) %*% u.temp
x <- score.temp * u.temp[1]
y <- score.temp * u.temp[2]
segments(x, y, mydata$panicle.length, mydata$leaf.length, col = "gray")
points(x, y, pch = 4, col = "green")

## focus on one sample
id <- which.max(mydata$leaf.length)
arrows(0, 0, mydata$panicle.length[id], mydata$leaf.length[id], col = "purple")
arrows(x[id], y[id], mydata$panicle.length[id], mydata$leaf.length[id], col = "pink")
arrows(0, 0, x[id], y[id], col = "blue")

#### 2. PCA without using prcomp ####
# calculate covariance matrix
cov <- var(mydata)
cov

# eigenvalue decomposition
eig <- eigen(cov)
eig

# compare results
res <- prcomp(mydata)
res
sqrt(eig$values)

# calculate principal component scores
mydata[1,]
eig$vectors[,1]
mydata[1,1] * eig$vectors[1,1] + mydata[1,2] * eig$vectors[2,1]
res$x[1,1]

score <- as.matrix(mydata) %*% eig$vectors
head(score)
head(res$x)

# variance of scores = eigenvalues 
var(score)
eig$values

# sum of variance 
sum(eig$values)
sum(diag(cov))

# contribution
eig$values / sum(eig$values)
cumsum(eig$values) / sum(eig$values)
summary(res)

#### 3. PCA based on a correlation matrix ####
# extract panicle length and florets per panicle
mydata <- data.frame(
	panicle.length  = alldata$Panicle.length,
	panicle.florets = alldata$Florets.per.panicle
	)
missing <- apply(is.na(mydata), 1, sum) > 0
mydata <- mydata[!missing, ]
# look at the relationship between two variables
plot(mydata)

# the following analysis is wrong
res <- prcomp(mydata)
res

# if panicle length is measured in meter unit
mydata$panicle.length <- mydata$panicle.length / 100
res.2 <- prcomp(mydata)
res.2

# scaling
mydata.scaled <- scale(mydata)
var(mydata.scaled)
res.scaled <- prcomp(mydata.scaled)
res.scaled

# cov and cor
eigen(cov(mydata.scaled))
eigen(cor(mydata))

# perform principal component analysis on scaled data
res.scaled.2 <- prcomp(mydata, scale = T)
res.scaled.2
res.scaled

#### 4. PCA with higher dimensional data ####
# multivariate (>3) analysis
mydata <- data.frame(
	leaf.length = alldata$Flag.leaf.length,
	leaf.width  = alldata$Flag.leaf.width,
	plant.height = alldata$Plant.height,
	panicle.number = alldata$Panicle.number,
	panicle.length = alldata$Panicle.length,
	seed.length = alldata$Seed.length,
	seed.width = alldata$Seed.width
	)
missing <- apply(is.na(mydata), 1, sum) > 0
mydata <- mydata[!missing, ]

# PCA based on a correlation matrix
res <- prcomp(mydata, scale = T)
summary(res)
plot(res)

# scatter plot principal component scores
subpop <- alldata$Sub.population[!missing]
op <- par(mfrow = c(1,2))
plot(res$x[,1:2], col = as.numeric(subpop))
legend("topleft", levels(subpop), col = 1:nlevels(subpop), pch = 1)
plot(res$x[,3:4], col = as.numeric(subpop))
par(op)
df <- data.frame(subpop = subpop, res$x[,1:3])
plot_ly(data = df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~subpop, type = "scatter3d", mode = "markers")

# understand the meaning of components (eigen vectors)
res$rotation[,1:4]

# understand the meaning of components (biplot)
op <- par(mfrow = c(1,2))
biplot(res, choices = 1:2)
biplot(res, choices = 3:4)
par(op)

# calculate factor loadings
# meaning of factor loadings (for PCA based on correlation matrix)
factor.loadings <- cor(mydata, res$x[,1:4])
factor.loadings

theta <- 2 * pi * (0:100 / 100)
x <- cos(theta)
y <- sin(theta)
op <- par(mfrow = c(1,2))
plot(factor.loadings[,1:2], xlim = c(-1,1), ylim = c(-1,1), pch = 4)
text(factor.loadings[,1:2], rownames(factor.loadings), col = "red")
lines(x, y, col = "gray")
abline(v = 0, h = 0)
plot(factor.loadings[,3:4], xlim = c(-1,1), ylim = c(-1,1), pch = 4)
text(factor.loadings[,3:4], rownames(factor.loadings), col = "red")
lines(x, y, col = "gray")
abline(v = 0, h = 0)
par(op)

factor.loadings <- t(res$sdev * t(res$rotation))[,1:4]
factor.loadings

#### 6. PCA with high dimensional data ####
# prepare multivariate data
mydata <- alldata[, 50:ncol(alldata)]
dim(mydata)
head(mydata)[,1:10]

# perform PCA
res.pca <- prcomp(mydata)
summary(res.pca)
plot(res.pca)

# plot principal component scores
subpop <- alldata$Sub.population
op <- par(mfrow = c(1,2))
plot(res.pca$x[,1:2], col = as.numeric(subpop))
plot(res.pca$x[,3:4], col = as.numeric(subpop))
legend(-10, 20, levels(subpop), col = 1:nlevels(subpop), pch = 1, cex = 0.5)
par(op)

# plot them with plotly
df <- data.frame(subpop = subpop, res.pca$x[,1:4])
plot_ly(data = df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~subpop, type = "scatter3d", mode = "markers")
plot_ly(data = df, x = ~PC2, y = ~PC3, z = ~PC4, color = ~subpop, type = "scatter3d", mode = "markers")

# correlation between PC1-4 in alldata on one hand and PC1-4 just calculated
cor(alldata[,c("PC1","PC2","PC3","PC4")], res.pca$x[,1:4])

#### 7. Multi-dimensional scaling (MDS) ####

# extract marker data
mydata <- alldata[, 50:ncol(alldata)]
D <- dist(mydata)

# perform MDS
res.mds <- cmdscale(D, k = 10, eig = T)
res.mds

# eigenvalues and contributions
res.mds$eig[1:10]
res.mds$eig[1:10] / sum(res.mds$eig)
cumsum(res.mds$eig[1:10]) / sum(res.mds$eig)
barplot(res.mds$eig[1:10])

# draw the result of MDS
subpop <- alldata$Sub.population
op <- par(mfrow = c(1,2))
plot(res.mds$points[,1:2], col = as.numeric(subpop))
plot(res.mds$points[,3:4], col = as.numeric(subpop))
legend(5, -10, levels(subpop), col = 1:nlevels(subpop), pch = 1, cex = 0.5)
par(op)

# draw the result of MDS with plotly
df <- data.frame(subpop = subpop, res.mds$points[,1:4])
plot_ly(data = df, x = ~X1, y = ~X2, z = ~X3, color = ~subpop, type = "scatter3d", mode = "markers")
plot_ly(data = df, x = ~X2, y = ~X3, z = ~X4, color = ~subpop, type = "scatter3d", mode = "markers")

## correlation between PC1-4 and scores in MDS
cor(res.pca$x[,1:4], res.mds$points[,1:4])

# prepare data
mydata <- data.frame(
	plant.height = alldata$Plant.height,
	res.mds$points[,1:4]
	)
mydata <- na.omit(mydata)
# analyze data
model <- lm(plant.height ~ ., data = mydata)
anova(model)
plot(mydata$plant.height, predict(model))

#### 8. MDS without using cmdscale ####
#prepare data again
mydata <- alldata[, 50:ncol(alldata)]
D <- dist(mydata)

#obtain B matrix
D2 <- as.matrix(D^2)
D2i. <- apply(D2, 1, mean)
D2.j <- apply(D2, 2, mean)
D2.. <- mean(D2)
B <- - 0.5 * (sweep(sweep(D2, 1, D2i.), 2, D2.j) + D2..)

#eigenvalue decomposition of B matrix
eig <- eigen(B)
eval <- eig$values[1:10]
evec <- eig$vectors[,1:10]
points <- evec * rep(sqrt(eval), each = nrow(evec))

# compare results
head(points, 4)
head(res.mds$points, 4)

# draw graph
subpop <- alldata$Sub.population
op <- par(mfrow = c(1,2))
plot(points[,1:2], col = as.numeric(subpop))
plot(points[,3:4], col = as.numeric(subpop))
legend(5, -10, levels(subpop), col = 1:nlevels(subpop), pch = 1, cex = 0.5)
par(op)

#### 9. Preparation for homework ####
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

