---
title: "Stat 440 0125"
author: "Sangjun Pyo"
date: "5/22/2022"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Monte Carlo Integration

Monte Carlo (MC) methods are a broad collection of techniques for performing numerical computations using simulated random numbers. Now that we’ve covered some fundamental simulation ideas, we can introduce our first major application area: numerical integration with MC. By the end of this lesson, you should be able to:

• Understand the theoretical justification behind MC methods.
• Implement MC estimation for integrals.

## Law of Large Numbers

Calculating integrals is of the most common computational problems in mathematics and statistics. While there are plenty of tools for analytically calculating integrals, many integrals do not admit nice analytic solutions. This is where simulation methods can help!

To understand how simulation can help, we need to review some connections between expected values and integrals.

## Monte Carlo Integral

Algorithm:
1. Choose a random variable with density f(x) that allows for easy sample simulation.
2. Choose h(x) = t(x)/f(x). Note that this expression must be finite for all x ∈ R.
3. Generate n iid samples with density f(x), X1, . . . , Xn.
4. Approximate

## Example 1: Uniform Sampling for Generic Functions
```{r }
# generate uniform random samples
set.seed(123)


ns = 1:1000
us = runif(1000)
hs = exp(-us)

# calculate Monte Carlo estimate at each sample size
mc_ests = cumsum(hs) / ns

# plot results
library(ggplot2)

ggplot(data=data.frame(n=ns, est=mc_ests)) +
  geom_line(mapping=aes(x=n, y=est), color="blue") +
  geom_hline(yintercept=(1 - exp(-1)), linetype="dashed") +
  xlab("Number of Monte Carlo Samples") +
  ylab("Integral approximation")

```

  The dashed line above represents our true integral value. Above, we see that
the Monte Carlo approximation to the integral converges on the true value, as
expected by the law of large numbers.

  We can use R to define a generic function for performing this kind of integration. Here we will use a nice property of R, in that functions can be used as inputs to other functions:

```{r }
bounded_mc_estimate = function(h, n, a, b) {
  #'
  #' Estimate the function h integrated from a to b using MC
  #'
  #' @param h function
  #' @param n number of MC samples
  #' @param a integral lower bound
  #' @param b integral upper bound
  
  # generate uniform samples
  us = runif(n, min=a, max=b)

  # return MC estimate
  (b - a) * sum(h(us)) / n
}

```

  Let’s use this function to estimate π numerically using MC. You can verify from calculus that:

$$ \int _{-1}^{1} 2 \sqrt{1-x^2}dx = 2 arcsin(1) = π $$
```{r}
bounded_mc_estimate(
  h=function(x) { 2 * sqrt(1 - x**2) },
  n=10000,
  a=-1,
  b=1
)
```

## Example 2: CDFs

In statistics, there are many distributions for which the PDF is easy to calculate, but the CDF does not have analytic closed-form equations. Examples include the Normal distribution, Gamma distribution (including χ^2 distribution), and Beta distribution, to name just a small few with this problem.

Alternatively, we can estimate the CDF by using MC integration by observing
that, for a random variable X with density f(x):
$$ $$

Therefore, if X1, . . . , Xn are iid with density f, then the MC approximation for
the CDF at a particular value x takes the form:

$$ $$

Not only does this work for CDFs, but also for the probability of any particular
event. Suppose there’s a set A ⊆ R for which P(X ∈ A) exists. Then:

$$ $$

Let’s use this to estimate the standard normal CDF at x = 1, i.e.:

$$ $$

```{r}
ns = 1:10000
xs = rnorm(10000)
hs = ifelse(xs <= 1, 1, 0)

# calculate Monte Carlo estimate at each sample size
mc_ests = cumsum(hs) / ns

# plot results
ggplot(data=data.frame(n=ns, est=mc_ests)) +
  geom_line(mapping=aes(x=n, y=est), color="blue") +
  geom_hline(yintercept=pnorm(1), linetype="dashed") +
  xlab("Number of Monte Carlo Samples") +
  ylab("Integral approximation")


```


  Let’s also use this to estimate P(X ∈ A), where A = {x ∈ R : |x| ≥ 2}:

$$P(|X| ≥ 2) ≈ .0455 $$

This corresponds to the probability that X is at least two standard deviations
away from the mean for any normal random variable.

```{r}
ns = 1:10000
xs = rnorm(10000)
hs = ifelse(abs(xs) >= 2, 1, 0)

# calculate Monte Carlo estimate at each sample size
mc_ests = cumsum(hs) / ns
true_value = 1 - (pnorm(2) - pnorm(-2))

# plot results
ggplot(data=data.frame(n=ns, est=mc_ests)) +
  geom_line(mapping=aes(x=n, y=est), color="blue") +
  geom_hline(yintercept=true_value, linetype="dashed") +
  xlab("Number of Monte Carlo Samples") +
  ylab("Integral approximation")

```

One last note: sometimes, we can convert integrals over unbounded domains to integrals over bounded domains. For example, when x ≥ 0 from our standard normal example:
$$ $$

When working with integrating functions like these, exploiting useful function information, such as symmetry, can make computation easier in practice.
