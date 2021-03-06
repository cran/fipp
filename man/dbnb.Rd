% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{dbnb}
\alias{dbnb}
\title{Probability density function of the BNB distribution}
\usage{
dbnb(x, mu, a, b, log = FALSE)
}
\arguments{
\item{x}{vector of quantiles.}

\item{mu}{mean parameter.}

\item{a}{1st shape parameter.}

\item{b}{2nd shape parameter.}

\item{log}{logical; if TRUE, density values p are given as log(p).}
}
\value{
Numeric vector of density values.
}
\description{
Evaluates the probability density function of the
beta-negative-binomial (BNB) distribution with a mean parameter and
two shape parameters.
}
\details{
The BNB distribution has density

\deqn{
  f(x) = \frac{\Gamma(\mu + x) B(\mu + a, x + b)}{\Gamma(\mu) \Gamma(x + 1) B(a, b)},
}{
  f(x) = (\Gamma(\mu + x) B(\mu + a, x + b)) / (Gamma(\mu) Gamma(x + 1) B(a, b)),
}
where \eqn{\mu} is the mean parameter and \eqn{a} and \eqn{b} are the first and
second shape parameter.
}
\examples{
## Similar to other d+DISTRIBUTION_NAME functions such as dnorm, it
## evaluates the density of a distribution (in this case the BNB distri)
## at point x
##
## Let's try with the density of x = 1 for BNB(1,4,3)
x <- 1
dbnb(x, mu = 1, a = 4, b = 3)

## The primary use of this function is in the closures returned from
## fipp() or nCluststers() as a prior on K-1
pmf <- nClusters(Kplus = 1:15, N = 100, type = "static",
gamma = 1, maxK = 150)

## Now evaluate above when K-1 ~ BNB(1,4,3)
pmf(priorK = dbnb, priorKparams = list(mu = 1, a = 4, b = 3))

## Compare the result with the case when K-1 ~ Pois(1)
pmf(priorK = dpois, priorKparams = list(lambda = 1))

## Although both BNB(1,4,3) and Pois(1) have 1 as their mean, the former
## has a fatter rhs tail. We see that it is reflected in the induced prior 
## on K+ as well

}
\references{
Frühwirth-Schnatter, S., Malsiner-Walli, G., and Grün, B. (2020)
Generalized mixtures of finite mixtures and telescoping sampling
\url{https://arxiv.org/abs/2005.09918}
}
