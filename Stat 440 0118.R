---
title: "Stat 440 0118"
author: "Sangjun Pyo"
date: "1/19/2022"
output: pdf_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

Univariate Normal Distribution

```{r}
# set up a seed for reproducibility

# try a different seed value for a different result

set.seed(123)

# sample from N(0,1) using rnorm
sims <- data.frame(x=rnorm(n=1000, mean=0, sd=1))

# compute densities of N(0,1) using dnorm
xset <- seq(-3, 3, .001)
zdat <- data.frame(x=xset,y=dnorm(xset, mean=0, sd=1))

# plot the results
library(ggplot2)

ggplot() +
  geom_histogram(data=sims,
                mapping=aes(x=x, y=..density..),
                bins=30, fill="white", color="black") +

  geom_line(data=zdat, color="blue", mapping=aes(x=x, y=y)) +
  ggtitle("N(0,1): simulation vs truth")

```

Box-Muller Transformation

```{r}
box_muller = function(n) {
  #' Apply Box-Muller to generate N(0,1) samples
  #' @param n output sample size
  
  # simulate independent U1 and U2 from Unif(0,1)
  u1 <- runif(n, min=0, max=1)
  u2 <- runif(n, min=0, max=1)

    # apply Box-Muller
  r <- sqrt(-2 * log(u1))
  t <- 2 * pi * u2

  x1 <- r * cos(t)
  x2 <- r * sin(t)

  # return independent X1 and X2 from N(0,1)
  return(list(x1,x2))
}

# run algorithm and plot results
box_muller_results = box_muller(n=10000)

ggplot() +
  geom_histogram(data=data.frame(x=box_muller_results[[1]]),
                mapping=aes(x=x, y=..density..),
                bins=30, fill="white", color="black") +
  geom_line(data=zdat, color="blue", mapping=aes(x=x, y=y)) +
  ggtitle("X1: simulation vs truth")

ggplot() +
  geom_histogram(data=data.frame(x=box_muller_results[[2]]),
                mapping=aes(x=x, y=..density..),
                bins=30, fill="white", color="black") +
geom_line(data=zdat, color="blue", mapping=aes(x=x, y=y)) +
ggtitle("X2: simulation vs truth")

ggplot(data=data.frame(x=box_muller_results[[1]],
                       y=box_muller_results[[2]]),
       mapping=aes(x=x,y=y)) +
  geom_bin2d() +
  geom_density_2d(color="white") +
  coord_fixed() +
  xlab("x1") +
  ylab("x2")

```

