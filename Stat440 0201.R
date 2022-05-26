---
title: "Stat 440 0201"
author: "Sangjun Pyo"
date: "5/26/2022"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Monte Carlo Standard Error

we illustrate methods for quantifying errors from MC integral
estimates based on the central limit theorem. By the end of this lesson, you
should be able to:

• Estimate uncertainty in MC integration estimates.
• Perform this uncertainty estimation in R.
• Reduce the uncertainty in MC integration estimates using Rao–Blackwellization.

## Example 1: Trigonometric Integral Standard Error

```{r }
library(ggplot2)
ggplot(data=data.frame(x=seq(-.01, 1.01, .001),
                        y=(cos(50*seq(-.01, 1.01, .001)) +
                        sin(20*seq(-.01, 1.01, .001)))**2),
      mapping=aes(x=x,y=y)) +
  geom_line() +
  geom_area(fill="blue", alpha=.36)
```

```{r}
# Define the function h(x)
h=function(x){(cos(50*x)+sin(20*x))ˆ2}

# Calculate the true value of the integral
# (value from Wolfram Mathematica)
true_integral = (
  (
    8240 - 105 * sin(40) + 42 * sin(100) +
    280 * cos(30) - 120 * cos(70)
  ) / 8400
)

# Simulate 10,000 random numbers from U(0,1)
x=h(runif(1e4,min=0,max=1))

# Estimate E[h(X)] using a naive Monte Carlo
estint=cumsum(x)/(1:10ˆ4)

# Estimate standard error, i.e., sqrt(Var)
esterr=sqrt(cumsum((x-estint)ˆ2))/(1:10ˆ4)

ggplot(data=data.frame(samps=1:10ˆ4,
                      ests=estint,
                      errs=esterr)) +

  # integral estimate
  geom_line(aes(x=samps, y=ests)) +

  # integral estimate plus 2 standard error upper bound
  geom_line(aes(x=samps, y=ests+2*errs), color="blue") +

  # integral estimate plus 2 standard error lower bound
  geom_line(aes(x=samps, y=ests-2*errs), color="red") +
  geom_hline(yintercept=true_integral,
              linetype="dashed", alpha=.5) +
              xlab("Number of MC samples") +
              ylab("Integral estimate +/- 2*SE")

```

## Rao–Blackwellization

This theorem also has bearings on computational methodology since it gives a
generic approach to reducing the variance of a Monte Carlo estimator, which is
to use conditioning. This technique is sometimes called Rao–Blackwellization,
although the conditioning is not always in terms of sufficient statistics in
Monte Carlo settings. It basically states that using conditional expectations—
that can be computed—in Monte Carlo expressions brings an improvement in
the variability of those Monte Carlo estimators while it does not perturb their
inherent unbiasedness

## Example 2: Bivariate Normal Distribution
```{r}
# library for historical statistics data sets
library(HistData)

# parameters for MVN for Galton data
model_mu = colMeans(Galton)
model_sigma = cov(Galton)

# function from MVN lecture to generate iid samples
mvrnorm_chol<-function(n,mu,Sigma){
  #' sample from N_d(mu, Sigma) with Cholesky decomposition
  #' @param n sample size
  #' @param mu mean vector
  #' @param Sigma covariance matrix
  
  # apply Cholesky decomposition to Sigma
  Sig_chol<-chol(Sigma)

  # determine the dimension
  p<-length(mu)

  # sample from N_d(mu, Sigma)
  X <- mu + t(Sig_chol)%*%matrix(rnorm(n*p),nrow=p,ncol=n)

  # return the sample
  return(t(X))
}

galton_naive_mc = function(n) {
  #' Naive MC approach for estimating p(X <Y)
  #' @param n MC sample size

  # generate iid MVN samples
  sample = mvrnorm_chol(n, model_mu, model_sigma)

  # estimate probability as mean of indicators
  mean(sample[, 1] < sample[, 2])
}

galton_rb_mc = function(n) {
  #' Rao-Blackwellized approach for estimating p(X <Y)
  #' @param n MC sample size

  # sample from child heights
  y = rnorm(n, model_mu[2], sqrt(model_sigma[2, 2]))

  # averaging over all samples...
  mean(
    # calculate the conditional probabilities using pnorm
    pnorm(
      y,
      # conditional mean
      mean=(
        model_mu[1] +
          model_sigma[2, 1] / model_sigma[2, 2] * (y - model_mu[2])
      ),
      # conditional SD
      sd=sqrt(
        (1 - model_sigma[1, 2]ˆ2 /
          (model_sigma[1,1] * model_sigma[2,2])) *
        model_sigma[1, 1]
      )
    )
  )
}

naive_mc_ests = replicate(10000, galton_naive_mc(1000))
rb_mc_ests = replicate(10000, galton_rb_mc(1000))

ggplot(data=data.frame(
  method=c(rep("NaiveMC", 10000), rep("RaoBlackwell", 10000)),
  sample=c(naive_mc_ests, rb_mc_ests)
)) +
  geom_histogram(aes(x=sample, y=..density.., fill=method), bins=30, alpha=.64)

var(naive_mc_ests)

var(rb_mc_ests)

```

## Example 3: Dickey’s Decomposition

```{r}
set.seed(999)
# Set true parameter values
nu=5; mu=3; sigma=0.5

# Set the number of simulations
Nsim=10e3

# Sample Y 10,000 times
y=sqrt(rchisq(Nsim,df=nu)/nu)

# Sample X given Y 10,000 times
x=rnorm(Nsim,mu,sigma/y)

# Estimate E[h(X)] by the empirical average
mean(exp(-xˆ2))

# Estimate E[h(X)] by the Rao–Blackwellized average
mean(exp(-muˆ2/(1+2*(sigma/y)ˆ2))/sqrt(1+2*(sigma/y)ˆ2))

xt <- mu + sigma * rt(Nsim, df=nu)
mean(exp(-xtˆ2))

# Estimate the variance of the empirical average
var(exp(-xˆ2))

# Estimate the variance of the Rao–Blackwellized average
var(exp(-muˆ2/(1+2*(sigma/y)ˆ2))/sqrt(1+2*(sigma/y)ˆ2))

# Visualize the convergence of two estimators for E[h(X)]
avg_est=cumsum(exp(-xˆ2))/(1:Nsim)
rba_est=cumsum(exp(-muˆ2/(1+2*(sigma/y)ˆ2))/sqrt(1+2*(sigma/y)ˆ2))/(1:Nsim)

library(ggplot2)

ggplot(data=data.frame(n=1:Nsim,avg_est=avg_est,rba_est=rba_est)) +
  geom_line(mapping=aes(x=n, y=avg_est), color="red", linetype="dotted") +
  geom_line(mapping=aes(x=n, y=rba_est), color="blue") +
  geom_hline(yintercept=0.00771, linetype="dashed") +
  xlab("Number of simulations") + ylab("Estimated E[h(X)]")

```
