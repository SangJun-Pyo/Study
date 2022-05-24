---
title: "Stat 440 0127"
author: "Sangjun Pyo"
date: "5/23/2022"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Importance Sampling

In this lecture, we introduce importance sampling, an alternative approach to
MC integration that relies on a clever intepretation of the MC integrand. By the
end of this lesson, you should be able to:

• Understand how to evaluate an integral using importance sampling.
• Choose a reasonable proposal distribution for importance sampling.
• Write numerically stable codes to address the potential underflow issue.

## Importance Sampling

Suppose we want to compute the expectation of a function h(X) of a random
variable X with respect to a distribution with density f(x)

where the expectation is evaluated over the density f. We could approximate I
by a “naive simulation”

where $X_1, . . . , X_M$ are independently sampled from the density f(x). From
our last lecture we have learned that $I^{MC}$ converges to I as $M → ∞$ by the law of large numbers. The main limitation of using this naive approach is that we may not know how to sample from f in many cases.

Let g(x) denote any other density function that we know how to sample
from. We further assume that that g(x) is non-zero whenever f(x) is non-zero
(We need this condition to avoid dividing by 0 in what follows).

where $X^`_1, . . . , X^`_M$ are independently sampled from the density g(x). This
new simulation is called “importance sampling” and g is called the “proposal
distribution”. Using the law of large numbers again, we can show that $I^{IS}$ also
converges to I as $M → ∞.$

## Choice of Proposal

The key idea behind importance sampling is that if the proposal distribution
g is well chosen then the approximation $\hat I^{IS}$ to I will be better than the naive
approximation $I^MC$. So the question then is how to choose a “good” proposal
distribution g. Here we provide some guiding principles and a theoretical result
that will help with this process.

First, we need to choose a density g with heavier tails than f. 
by “heavier tails”

As an example, if f(x) is the standard normal density and g(x) is the standard
Laplace density, then g has “heavier tails” than f because:

$$ $$


We can also show this visually:
```{r }
xseq <- seq(-5, 5, .01)

library(ggplot2)

ggplot() +
  geom_line(data=data.frame(x=xseq, y=dnorm(xseq)),
            mapping=aes(x=x,y=y),
            color="blue",
            linetype="dashed") +
  geom_line(data=data.frame(x=xseq, y=1/2*exp(-abs(xseq))),
            mapping=aes(x=x,y=y),
            color="red",
            linetype="dashed") +
ylab("density")

```

Above, we see that the normal distribution density (blue) decays to 0 faster
than the Laplace distribution density (red).
  Choosing the proposal distribution with this property is quite important.
If this is not satisfied, $\hat I^{IS}$ might have an infinite variance. 
To see this, let’s compute the following second moment:

$$ $$

If g has lighter tails than f, then this integral might be infinite and so is the
variance of $\hat I^{IS}$

Second, we need to choose g to be similar in shape to f. To see this, suppose that g(x) is small over some set $A ⊆ R$ where f(x) is large. Again, the ratio of f(x)/g(x) could be large, leading to a large variance. The first two principles together suggest that a good choice for an importance sampling
density g should be similar to f but with heavier tails, as visualized below
(source: https://pure.tue.nl/ws/portalfiles/portal/2895107/657045.pdf).

As a theoretical exercise, we might ask what choice of density g explicitly
minimizes the variance of our importance sampling estimator. It turns out this
choice of g has a closed form expression:

$$ $$

## Example 1: Tail Probability of Standard Normal

Suppose $X ∼ N (0, 1)$, and we want to compute $P(X > z)$ for z = 10. That
is, $h(x) = I(x > 10)$ and f(x) = φ(x) is the density of the standard normal
distribution.
  Let’s try naive simulation $\hat I^{MC}$, and compare it with the “truth”, as ascertained by the R function pnorm.

```{r }
set.seed(999)

# naive simulation
X = rnorm(100000)
mean(1*(X>10))

# ground truth
pnorm(10,mean=0,sd=1,lower.tail=FALSE)

```

  You can see the problem with our naive simulation: all of our simulated
samples are less than 10 (where f(x) = 0), so we have an imprecise estimator
of P(X > z).
  Instead, we use importance sampling. Here we code the general case for
any value of z, setting the proposal distribution g to be N (z, 1). Note that
because of this choice of g, many of the samples are > z (where f is non-zero)
and we hope to get better precision.

```{r}
pnorm.IS = function(z, nsamp=100000){
  #' Importance sampling from upper tail of standard normal
  #' @param z upper tail value
  #' @param nsamp number of proposal samples

  # sample from g: N(z,1)
  y = rnorm(nsamp,z,1)

  # compute the ratio f/g
  # f/g = exp(log(f/g)) = exp[log(f)-log(g)]
  w = exp(dnorm(y,mean=0,sd=1,log=TRUE)
          - dnorm(y,mean=z,sd=1,log=TRUE))

  # compute the importance sampling estimate
  return(mean(w*(y>z)))
}

# importance sampling
pnorm.IS(10)

# ground truth
pnorm(10,mean=0,sd=1,lower.tail=FALSE)
```

  As expected, our importance sampling method does a better job of calculating this probability. Methods like these are very common in the study of statistically rare events, such as those arising in climate science and economics.
  
## Example 2: Computing Means on Log Scale

We now push Example 1 a bit further, to illustrate a numerical issue that can
arise quite generally (not just for importance sampling). Let’s try the above
with z = 100.

```{r}

pnorm.IS(100)
pnorm(100,lower.tail=FALSE)

```
When z = 100 both the naive approach and importance sampling produce
exactly zero, whereas from our probability class we learned that P(X >
100) 6= 0. This numerical issue is known as underflow (source: https://
en.wikipedia.org/wiki/Arithmetic_underflow), where R mistakenly
represents an extremely small non-zero number like P(X > 100) as an exact
zero.
The trick to solving this underflow issue is doing things on log scale. But
the importance sampling involves averaging, and we have to do the averaging
on the raw scale, not the log scale. The following function allows us to do this
computation:

```{r}
lmean = function(lx) {
  m = max(lx)
  return(m + log(mean(exp(lx-m))))
}

```
  Exploiting this trick, we can now perform importance sampling for z = 100
in a more stable way:

```{r}
lpnorm.IS = function(z,nsamp=100000){
  #' Importance sampling from upper tail of standard normal
  #' using the log-scale computation
  #' @param z upper tail value
  #' @param nsamp number of proposal samples

  # sample from g: N(z,1)
  y = rnorm(nsamp,z,1)

  # compute the ratio f/g on the log scale
  # i.e., log(f/g) = log(f) - log(g)
  w = dnorm(y,0,1,log=TRUE) - dnorm(y,z,1,log=TRUE)

  # compute the importance sampling estimate
  # note that mean(w*(y>z)) = mean(y>z) * mean(w[y>z])
  return(log(mean(y>z)) + lmean(w[y>z]))
}

# importance sampling
lpnorm.IS(100)

# ground truth
pnorm(100,mean=0,sd=1,lower.tail=FALSE,log=TRUE)
```

