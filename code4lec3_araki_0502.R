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
save(res, subpop, file = "ts_araki_0502.rda")
plot(res$x[,1:2], col = as.numeric(subpop))
legend("topleft", levels(subpop), col = 1:nlevels(subpop), pch = 1)
plot(res$x[,3:4], col = as.numeric(subpop))
par(op)
df <- data.frame(subpop = subpop, res$x[,1:3])
plot_ly(data = df, x = ~PC1, y = ~PC2, z = ~PC3, color = ~subpop, type = "scatter3d", mode = "markers")