## An Introduction to Biostatistical Analysis with R
## R codes for lecture on 04/19/2019

## 1st step: Begin with the ABCs of R

# simple calculation
3 + 5 * 3

# assign value to a varible
x <- 1 + 2
x

# calculation with the value assigned variable
x + 5 * x

# calculation with various functions
abs(x)
sin(x)
atan(x)
log(x)
log10(x)

# calculate PDF of normal distribution with mean mu and sd sigma
mu <- 3
s2 <- 2
x <- 5
1 / sqrt(2 * pi * s2) * exp(- (x - mu)^2 / (2 * s2)) 

# calculate PDF of normal distribution using the R fanction dnorm
dnorm(x, mu, sqrt(s2))

# vectorized arithmetic
length <- c(8.1, 7.7, 8.2, 9.7, 7.1, 7.3)    # mm scale
length

# calculate length-width ratio
width <- c(3.7, 3.0, 2.9, 2.4, 3.3, 2.5)
ratio <- length / width
ratio

# calculate mean of seed length
sum(ratio)
length(ratio)
sum(ratio) / length(ratio)

# calculate mean using the R function mean
mean(ratio)

# calculate variance
xbar <- mean(ratio)
(ratio - xbar)^2
sum((ratio - xbar)^2)
sum((ratio - xbar)^2) / (length(ratio) - 1)

# calculate standard deviation using the R function sd
var(ratio)

# calculate covariance
xbar <- mean(length)
ybar <- mean(width)
sum((length - xbar) * (width - ybar)) / (length(length) - 1)

# calculate covariance using the R function cov
cov(length, width)

# calculate correlation
s12 <- sum((length - xbar) * (width - ybar))
s1 <- sum((length - xbar)^2)
s2 <- sum((width - ybar)^2)
s12 / (sqrt(s1) * sqrt(s2))

# calculate correlation with cov and sd
cov(length, width) / (sd(length) * sd(width))

# calculate correlation using the R function cor
cor(length, width)

# combine two vectors to make a 6x2 matrix
x <- cbind(length, width)
x

# calculate column means
m <- apply(x, 2, mean)
m

# subtract column mean from each column
z <- sweep(x, 2, m)
z

# calculate a covariance matrix
t(z) %*% z / (nrow(z) - 1)


# read rice phenotype data
# this data set was analyzed in Zhao 2011 (Nature Communications 2:467)
pheno <- read.csv("RiceDiversityPheno.csv")

# check the size of data
dim(pheno)

# look at the first 6 lines of data
head(pheno)

# convine the dataset with another dataset
line <- read.csv("RiceDiversityLine.csv")
head(line)

# merge two datasets
data <- merge(line, pheno, by.x = "NSFTV.ID", by.y = "NSFTVID")
head(data)


# calculate ratio
ratio <- data$Seed.length / data$Seed.width
mean(ratio)
ratio

# calculate mean again
mean(ratio, na.rm = T)

# calculate mean of each column
sapply(data, mean, na.rm = T)

# summary of the pheno dataset 
summary(data)

# calculate correation between seed length and widht
cor(data$Seed.length, data$Seed.width)

# calculate correation again with an option for trearing missing values
cor(data$Seed.length, data$Seed.width, use = "pair")

## controll statement

# Fibonacci sequence
fibSeq <- function(n) {
	fs <- c(0, 1, rep(NA, n - 1))
	for(i in 3:(n + 1)) {
		fs[i] <- fs[i - 2] + fs[i - 1]
	}
	return(fs)
}

fibSeq(15)

# Fibonacci again
fib <- function(n) {
	if(n < 2) {
		return(n)
	} else {
		return(fib(n - 1) + fib(n - 2)) # calculate recursively
	}
}
for(i in 0:15) {
	print(fib(i))
}

# prime number
isPrime <- function(n) {
	if(n == 2) {
		return(FALSE)	
	}
	
	for(i in 2:(n-1)) {
		if(n %% i == 0) {
			return(FALSE)
		}
	}
	return(TRUE)
}

isPrime(11)
isPrime(111)


## 2nd step: Visualization
# attach the database to the R search path
attach(data)

## 1d plot
# draw histogram
hist(Plant.height)

# draw stem-and-leaf plot 
stem(Plant.height)

# draw boxplot
boxplot(Plant.height)

# draw histogram again for Blast resistance
hist(Blast.resistance)

# create a contigency table for Blast resistance 
t <- table(Blast.resistance)
t

# draw a barplot for table object
plot(t)

# change the title of axes
plot(t, xlab = "Blast resistance scores", ylab = "Frequency")

# draw another type of barplot
barplot(t)

# draw pie chart
pie(t)

# add the main title
pie(t, main = "Blast resistance")


## 2d plot
# scatter plot
plot(Plant.height, Panicle.length)

# draw regression line
abline(lm(Panicle.length ~ Plant.height))

# draw rug for both axes
rug(Plant.height, side = 1)
rug(Panicle.length, side = 2)

# draw Scatterplot with boxplots showing the marginal distributions
def.par <- par(no.readonly = T)	
layout(matrix(c(2, 0, 1, 3), nrow = 2, byrow = T), widths = c(2, 1), heights = c(1, 2), respect = T)
plot(Plant.height, Panicle.length)
boxplot(Plant.height, horizontal = T)
boxplot(Panicle.length)
par(def.par)

## 3d plot
# draw scatterplot with estimated bivariate density
require("KernSmooth")
x <- data.frame(Plant.height, Panicle.length)
x <- na.omit(x)
d <- bkde2D(x, bandwidth = 4)
plot(x)
contour(d$x1, d$x2, d$fhat, add = T)

# draw perspective plot of estimated bivariate density
persp(d$x1, d$x2, d$fhat, xlab = "Plant.height", ylab = "Panicle.length", zlab = "density", theta = -30, phi = 30)

## look at the relationship between genetic background and phenotype
pop.id <- as.numeric(Sub.population)
plot(Plant.height, Panicle.length, col = pop.id)
levels(Sub.population)
legend("bottomright", levels(Sub.population), col = 1:nlevels(Sub.population), pch = 1)

# draw boxplot for each sub population
boxplot(Plant.height ~ Sub.population)
# add color to the boxplot
boxplot(Plant.height ~ Sub.population, border = 1:nlevels(Sub.population))

# draw Scatterplot with boxplots showing the marginal distributions
def.par <- par(no.readonly = T)	
layout(matrix(c(2, 0, 1, 3), nrow = 2, byrow = T), widths = c(2, 1), heights = c(1, 2), respect = T)
plot(Plant.height, Panicle.length, col = pop.id)
legend("bottomright", levels(Sub.population), col = 1:nlevels(Sub.population), pch = 1)
boxplot(Plant.height ~ Sub.population, border = 1:nlevels(Sub.population), horizontal = T)
boxplot(Panicle.length ~ Sub.population, border = 1:nlevels(Sub.population))
par(def.par)

## relationships among more than three variables
# draw bubble plot
symbols(Plant.height, Panicle.length, circles = Flag.leaf.length, inches = 0.1, fg = pop.id)

# draw a scatter plot matrix
x <- data.frame(Plant.height, Panicle.length, Flag.leaf.length)
pairs(x, col = pop.id)

# draw a bit complex scater plot matrix
pairs(x, panel = function(x, y, ...) {
	points(x, y, ...)
	abline(lm(y ~ x), col = "gray")
}, col = pop.id)

# draw 3d scatter plot using scatterplot3d
require(scatterplot3d)
scatterplot3d(Plant.height, Panicle.length, Flag.leaf.length, color = pop.id)

## draw world map and ovarlay origins of germplasm
require(maps)
require(mapdata)
map('worldHires')
points(Longitude, Latitude, col = pop.id)
legend("bottomleft", levels(Sub.population), col = 1:nlevels(Sub.population), pch = 1)

map('worldHires')
points(jitter(Longitude, 200), Latitude, col = pop.id)
legend("bottomleft", levels(Sub.population), col = 1:nlevels(Sub.population), pch = 1)

# output figure to a PDF file
pdf("map.pdf")
map('worldHires')
points(jitter(Longitude, 200), Latitude, col = pop.id)
legend("bottomleft", levels(Sub.population), col = 1:nlevels(Sub.population), pch = 1)
dev.off()

# output figure to a PDF file
pdf("map_large.pdf", width = 20, height = 10)
map('worldHires')
points(jitter(Longitude, 200), Latitude, col = pop.id)
legend("bottomleft", levels(Sub.population), col = 1:nlevels(Sub.population), pch = 1)
dev.off()

## Use plotly
## plotly is a package for drawing interactive graphs
# let's draw a histogram with the plot_ly function
plot_ly(x = Plant.height, type = "histogram")

# you can easily draw a horizontal histotram
plot_ly(y = Plant.height, type = "histogram")

# draw a boxplot
plot_ly(y = Plant.height, type = "box")

# draw a boxplot separatedly for each subpopulaiton 
plot_ly(y = Plant.height, color = Sub.population, type = "box")

# draw a barplot based on table data
t <- table(Blast.resistance) # do it again, just in case
plot_ly(x = names(t), y = t, type = "bar")

# draw a pie chart
plot_ly(labels = names(t), values = t, type = "pie")
# add title to the chart
plot_ly(labels = names(t), values = t, type = "pie") %>% 
layout(title = "Blast resisntance")

# draw a scatter plot for two variables
plot_ly(x = Plant.height, y = Panicle.length, type = "scatter", mode = "markers")
# add colors to plots
plot_ly(x = Plant.height, y = Panicle.length, color = Sub.population, type = "scatter", mode = "markers")
# change the size of plots
plot_ly(x = Plant.height, y = Panicle.length, color = Sub.population, size = Flag.leaf.length, type = "scatter", mode = "markers")

# in the case using data.frame
df <- data.frame(Sub.population, Plant.height, Panicle.length, Flag.leaf.length)
df <- na.omit(df)
plot_ly(data = df, x = ~Plant.height, y = ~Panicle.length, color = ~Sub.population, size = ~Flag.leaf.length, type = "scatter", mode = "markers")

# draw a 3D smoothed density surface
plot_ly(data = data.frame(d), x = d$x1, y = d$x2, z = d$fhat) %>%
add_surface()

# draw a 3D scatter plot
plot_ly(data = df, x = ~Plant.height, y = ~Panicle.length, z = ~Flag.leaf.length, color = ~Sub.population, type = "scatter3d", mode = "markers")

