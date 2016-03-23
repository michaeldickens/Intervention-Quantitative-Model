require(fitdistrplus)

DF <- read.csv("data.csv")

basicPlot <- function() {
    plot(sort(DF[["Value"]]), col="blue", log="")
}

basicHistogram <- function() {
    hist(DF[["Value"]], breaks=100)
}

fit <- function() {
    ## Distribution names: norm, lnorm, unif
    res = fitdist(DF[["Value"]], "norm")
    print(summary(res))
    plot(res)
}

fit()
