rm(list=ls())
library(LaplacesDemon)
library(bayesm)
library(data.table)
library(knitr)
library(flextable)
library(tidyverse)
library(tmvtnorm)
library(car)
library(rgl)
# Reject-Accept Algorithm for (standard) truncated normal.
# 1) draw X* from g(X) rnorm(mean=0,sd=1)
# 2) accept X* if X*>=lower truncation and X*<=upper truncation, else reject

sigma <- matrix(c(1,0.5,0.5,1), ncol=2)
x <- rtmvnorm(n=10000, mean=c(0,0), sigma=sigma, upper=c(Inf, 0), lower = c(0,-Inf), algorithm="gibbs", burn.in.samples=100)
plot(x, main="samples from truncated bivariate normal distribution",
     xlim=c(-6,6), ylim=c(-6,6), 
     xlab=expression(x[1]), ylab=expression(x[2]))
abline(v=1, lty=3, lwd=2, col="red")
abline(h=0, lty=3, lwd=2, col="blue")

x <-data.frame(x)
first <- density(x$X1)
plot(first)
second <- density(x$X2)
plot(second)
joint.density.plot(x$X1, x$X2, Title=NULL, contour=TRUE, color=FALSE, Trace=NULL)

