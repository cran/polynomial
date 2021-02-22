\name{exponents}
\alias{exponents}
\alias{exponents<-}
\title{
Polynomial Exponents
}
\description{
Functions to get or set the exponents of a polynomial.
}
\usage{
exponents(x)
exponents(x) <- value
}
\arguments{
  \item{x}{
a polynomial, or an object which can be coerced by \code{\link{as.polynomial}}
to a polynomial.
}
  \item{value}{
an integer vector of the same length as \code{x}, or \code{NULL}.
}
}
\details{
\code{exponents} and \code{exponents<-} get and set the \code{"exponents"}
attribute of a polynomial.
}
\value{
For \code{exponents}, an integer vector of the same length as \code{x}.
(\code{NULL} indicates that the polynomial has default exponents, given by
\code{seq(0L, along.with = x)}).

For \code{exponents<-}, the updated object. (Note that the value of
\code{exponents(x) <- value} is that of the assignment, \code{value}, not the
return value from the left-hand side.)
}
\examples{
x <- 1:10


## 'x' has NULL exponents, representing the default exponents
exponents(x)


## using 'exponents<-' on an object that is not
## a polynomial coerces the object to polynomial
exponents(x) <- 9:0
is.polynomial(x)
x
}
