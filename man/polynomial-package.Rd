\name{polynomial-package}
\alias{polynomial-package}
\docType{package}
\title{
\packageTitle{polynomial}
}
\description{
\packageDescription{polynomial}
}
\details{
A set of functions for some common polynomials have been provided, they can be
found in \code{\link{CommonPolynomials}}.

The most important functions are:

\code{\link{as.polynomial}} to turn an object into a polynomial.

\code{\link[polynomial:polynomial-class]{coef}}, \code{\link{coef<-}},
\code{\link{exponents}} and \code{\link{exponents<-}} for extracting and
replacing polynomial coefficients and exponents.

\code{\link{reduce}} for simplifying polynomials.

\code{\link[polynomial:polynomial-class]{as.function}} to turn a polynomial into
a function.

\code{\link[polynomial:plot.polynomial]{plot}} for visualizing a polynomial.

Arithmetic is defined for polynomials (unary and binary \code{+} and \code{-},
\code{*}, \code{/}, \code{^}), but not in the same way as the base arithmetic
operators.
}
\author{
\packageAuthor{polynomial}

Maintainer: \packageMaintainer{polynomial}
}
\keyword{ package }
