## An Introduction to Biostatistical Analysis with R
## R codes for lecture on 04/24/2020

# prepare required packages
required.packages <- c("plotly")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  install.packages(new.packages, repos="http://cran.us.r-project.org")
}
# load required packages
require(plotly)

## 2rd step: Find relationship between variables using regression

# this data set was analyzed in Zhao 2011 (Nature Communications 2:467)
pheno <- read.csv("RiceDiversityPheno.csv")
line <- read.csv("RiceDiversityLine.csv")
line.pheno <- merge(line, pheno, by.x = "NSFTV.ID", by.y = "NSFTVID")
head(line.pheno)[,1:12]

# extract variables for regression analysis
data <- data.frame(
	height = line.pheno$Plant.height,
	flower = line.pheno$Flowering.time.at.Arkansas,
	PC1 = line.pheno$PC1,
	PC2 = line.pheno$PC2,
	PC3 = line.pheno$PC3,
	PC4 = line.pheno$PC4)
data <- na.omit(data)
head(data)

# look at the relationship between plant height and flowering time
plot(data$height ~ data$flower)

# perform single linear regression
model <- lm(height ~ flower, data = data)
summary(model)

# again, plot the two variables
plot(data$height ~ data$flower)
abline(model, col = "red")

# calculate fitted values
height.fit <- fitted(model)
points(data$flower, height.fit, pch = 3, col = "green")

# plot residuals
segments(data$flower, height.fit,
		 data$flower, height.fit + resid(model), col = "gray")
# height = height.fit + resid(model)

# predict unknown data
height.pred <- predict(model, data.frame(flower = seq(60, 140, 20)))
points(seq(60, 140, 20), height.pred, pch = 2, col = "blue")

######
# visualize the plane for optimization
x <- data$flower
y <- data$height
mu <- seq(0, 100, 1)
beta <- seq(0, 2, 0.02)
sse <- matrix(NA, length(mu), length(beta))
for(i in 1:length(mu)) {
	for(j in 1:length(beta)) {
		sse[i, j] <- sum((y - mu[i] - beta[j] * x)^2)
	}

}
persp(mu, beta, sse, col = "green")

# draw the figure with plotly
df <- data.frame(mu, beta, sse)
plot_ly(data = df, x = ~mu, y = ~beta, z = ~sse) %>% add_surface()

# calculate sum of squares (ss) of x and ss of xy
n <- length(x)
ssx <- sum(x^2) - n * mean(x)^2
ssxy <- sum(x * y) - n * mean(x) * mean(y)

# calculate b
b <- ssxy / ssx
b

# calculate m
m <- mean(y) - b * mean(x)
m

# draw scatter plot and regression line 
plot(y ~ x)
abline(m, b)

# calculate fitted values
y.hat <- m + b * x
lim <- range(c(y, y.hat))
plot(y, y.hat, xlab = "Observed", ylab = "Fitted", xlim = lim, ylim = lim)
abline(0, 1)

# calculate correlation between observed and fitted values
cor(y, y.hat)

# compare the square of the correlation and R2 
cor(y, y.hat)^2
summary(model)

###### 
model <- lm(height ~ flower, data = data)

# analysis of variance of regression
anova(model)

### ANOVA--ANOVA--ANOVA...
# calculate sum of squares of regression and error
ssr <- b * ssxy
ssr
ssy <- sum(y^2) - n * mean(y)^2
sse <- ssy - ssr
sse

# calculate mean squares of residuals
msr <- ssr / 1
msr
mse <- sse / (n - 2)
mse

# calculate F value
f.value <- msr / mse
f.value
# calculate p value for the F value
1 - pf(f.value, 1, n - 2)

# check the summary of the result of regression analysis		 
summary(model)

# square root of mse
sqrt(mse)

# R squared
ssr / ssy

# adjusted R squared
(ssy / (n - 1) - mse) / (ssy / (n - 1))

# confidnence interval of parameters

# test beta = 0
t.value <- (b - 0) / sqrt(mse/ssx)
t.value

# visualize the t distribution under H0
s <- seq(-10, 10, 0.2)
plot(s, dt(s, n - 2), type = "l")
abline(v = t.value, col = "green")
abline(v = - t.value, col = "gray")

# perform t test
2 * (1 - pt(abs(t.value), n - 2))

# check the summary of the model
summary(model)

# test beta = 0.5
t.value <- (b - 0.5) / sqrt(mse/ssx)
t.value
2 * (1 - pt(abs(t.value), n - 2))

# test mu = 0
t.value <- (m  - 0) / sqrt(mse * (1/n + mean(x)^2 / ssx))
t.value
2 * (1 - pt(abs(t.value), n - 2))

# check the summary of the model again
summary(model)

# test mu = 70
t.value <- (m  - 70) / sqrt(mse * (1/n + mean(x)^2 / ssx))
t.value
2 * (1 - pt(abs(t.value), n - 2))

# check the summary of the model
summary(model)

# visualize the t distribution under H0
s <- seq(-5, 5, 0.1)
plot(s, dt(s, n - 2), type = "l")
abline(v = t.value, col = "green")
abline(v = - t.value, col = "gray")

## confidence interval of estimates of y
# fitted values
pred <- predict(model)
head(pred)
head(fitted(model))

# calculate confidence interval
pred <- predict(model, interval = "confidence", level = 0.95)
head(pred)

# draw confidence bands
pred <- data.frame(flower = 50:160)
pc <- predict(model, interval = "confidence", newdata = pred)
plot(data$height ~ data$flower)
matlines(pred$flower, pc, lty = c(1, 2, 2), col = "red")

# draw confidence bands
pc <- predict(model, interval= "prediction", newdata = pred)
plot(data$height ~ data$flower)
matlines(pred$flower, pc, lty = c(1, 2, 2), col = "green")

# estimate the confidence intervals for the estimate and prediction of y 
pred <- data.frame(flower = 120)
predict(model, interval = "confidence", newdata = pred, level = 0.99)
predict(model, interval = "prediction", newdata = pred, level = 0.99)


# polynomial regression
model.quad <- lm(height ~ flower + I(flower^2), data = data)
summary(model.quad)

# plot(data$height ~ data$flower)
pred <- data.frame(flower = 50:160)
pc <- predict(model.quad, int = "c", newdata = pred)
plot(data$height ~ data$flower)
matlines(pred$flower, pc, lty = c(1, 2, 2), col = "red")

# compare predicted and observed values
lim <- range(c(data$height, fitted(model), fitted(model.quad)))
plot(data$height, fitted(model), 
		xlab = "Observed", ylab = "Expected",
	 	xlim = lim, ylim = lim)
points(data$height, fitted(model.quad), col = "red")
abline(0, 1)

# compare error variance between two models
anova(model, model.quad)

# extend polynormial regression model to a higher dimensional one...
model.cube <- lm(height ~ flower + I(flower^2) + I(flower^3), data = data)
summary(model.cube)
# compare error variance between two models
anova(model.quad, model.cube)

# multi-linear regression with genetic background 
model.wgb <- lm(height ~ PC1 + PC2 + PC3 + PC4, data = data)
summary(model.wgb)
anova(model.wgb)

# multi-linear regression with all information
model.all <- lm(height ~ flower + I(flower^2) + PC1 + PC2 + PC3 + PC4, data = data)
summary(model.all)
# compare error variance
anova(model.all, model.wgb)

# compare between the simplest and final models
lim <- range(data$height, fitted(model), fitted(model.all))
plot(data$height, fitted(model), xlab = "Observed", ylab = "Fitted", xlim = lim, ylim = lim)
points(data$height, fitted(model.all), col = "red")
abline(0,1)


## simulations of a yield experiment with randomized complete block design
# set a seed for random number generation
set.seed(12)

# simulate the condition of an experimental field with 4 x 4 plots
# The blocks have unequal fertility among them
field.cond <- matrix(rep(c(4,2,-2,-4), each = 4), nrow = 4)
field.cond

# set block to consider the heterogeneity of field condition
block <- c("I", "II", "III", "IV")
blomat <- matrix(rep(block, each = 4), nrow = 4)
blomat

# assume that there are four varieties
variety <- c("A", "B", "C", "D")
# sample the order of the four varieties randomly
sample(variety)
sample(variety)

# allocate the varieties randomly to each column of the field
varmat <- matrix(c(sample(variety), sample(variety), 
			sample(variety), sample(variety)), nrow = 4)
varmat 

# simulate genetic ability of the varieties
g.value <- matrix(NA, 4, 4)
g.value[varmat == "A"] <- 4
g.value[varmat == "B"] <- 2
g.value[varmat == "C"] <- -2
g.value[varmat == "D"] <- -4
g.value

# simulate error variance (variation due to the heterogeneity of local environment)
e.value <- matrix(rnorm(16, sd = 2.5), 4, 4)
e.value

# simulate phenotypic values
grand.mean <- 50
simyield <- grand.mean + field.cond + g.value + e.value
simyield

# unfold a matrix to a vector
as.vector(simyield)
as.vector(varmat)
as.vector(blomat)

# create a dataframe for the analysis of variance
simdata <- data.frame(variety = as.vector(varmat), block = as.vector(blomat), yield = as.vector(simyield))
simdata

# draw interaction plot
interaction.plot(simdata$block, simdata$variety, simdata$yield)

# perform the analysis of variance (ANOVA) with simulated data
res <- aov(yield ~ block + variety, data = simdata)
summary(res)

# perform ANOVA with a linear model
res <- lm(yield ~ block + variety, data = simdata)
anova(res)

## simulations of a yield experiment with completely randomized design (CRD)
# completely randomized the plots of each variety in a field
varmat.crd <- matrix(sample(varmat), nrow = 4)
varmat.crd

# simulate genetic ability of the varieties
g.value.crd <- matrix(NA, 4, 4)
g.value.crd[varmat.crd == "A"] <- 4
g.value.crd[varmat.crd == "B"] <- 2
g.value.crd[varmat.crd == "C"] <- -2
g.value.crd[varmat.crd == "D"] <- -4
g.value.crd

# simulate phenotypic values
simyield.crd <- grand.mean + g.value.crd + field.cond + e.value
simyield.crd

# create a dataframe for the analysis of variance
simdata.crd <- data.frame(variety = as.vector(varmat.crd), 
							yield = as.vector(simyield.crd))
simdata.crd

# perform ANOVA
res <- lm(yield ~ variety, data = simdata.crd)
anova(res)
summary(res)

# perform multiple simulations
n.rep <- 100
p.rbd <- rep(NA, n.rep)
p.crd <- rep(NA, n.rep)
for(i in 1:n.rep) {
	# experiment with randomized block design
	varmat <- matrix(c(sample(variety), sample(variety), 
			sample(variety), sample(variety)), nrow = 4)
	g.value <- matrix(NA, 4, 4)
	g.value[varmat == "A"] <- 4
	g.value[varmat == "B"] <- 2
	g.value[varmat == "C"] <- -2
	g.value[varmat == "D"] <- -4
	e.value <- matrix(rnorm(16, sd = 2.5), 4, 4)
	simyield <- grand.mean + field.cond + g.value + e.value
	simdata <- data.frame(variety = as.vector(varmat), 
			block = as.vector(blomat), yield = as.vector(simyield))
	res <- lm(yield ~ block + variety, data = simdata)
	p.rbd[i] <- anova(res)$Pr[2]
	
	# experiment with completed randomized design
	varmat.crd <- matrix(sample(varmat), nrow = 4)
	g.value.crd <- matrix(NA, 4, 4)
	g.value.crd[varmat.crd == "A"] <- 4
	g.value.crd[varmat.crd == "B"] <- 2
	g.value.crd[varmat.crd == "C"] <- -2
	g.value.crd[varmat.crd == "D"] <- -4
	simyield.crd <- grand.mean + g.value.crd + field.cond + e.value
	simdata.crd <- data.frame(variety = as.vector(varmat.crd), 
							yield = as.vector(simyield.crd))
	res <- lm(yield ~ variety, data = simdata.crd)
	p.crd[i] <- anova(res)$Pr[1]
}
sum(p.rbd < 0.05) / n.rep
sum(p.crd < 0.05) / n.rep
sum(p.rbd < 0.01) / n.rep
sum(p.crd < 0.01) / n.rep


############ Preparation of data for homework ############
# read two datasets and merge them
pheno <- read.csv("RiceDiversityPheno.csv")
line <- read.csv("RiceDiversityLine.csv")
line.pheno <- merge(line, pheno, by.x = "NSFTV.ID", by.y = "NSFTVID")

# extract variables for regression analysis
data <- data.frame(
	seed = line.pheno$Seed.number.per.panicle,
	length = line.pheno$Panicle.length)
data <- na.omit(data)
head(data)

