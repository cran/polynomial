\name{as.body}
\alias{as.body}
\title{
Create a Function Body From a Polynomial
}
\description{
Coerce a polynomial to class \code{"\linkS4class{call}"} such that it can be
used in function bodies and for labels in figures.
}
\usage{
as.body(x, var = "x", ...)
}
\arguments{
  \item{x}{
a polynomial, or an object which can be coerced by \code{\link{as.polynomial}}
to a polynomial.
  }
  \item{var}{
a character string; the name of the variable to use for the polynomial.
Optionally a \code{\link{call}} or \code{\link{expression}}, see
\strong{Examples}.
  }
  \item{\dots}{
further arguments passed to \code{\link{reduce}}.
  }
}
\details{
A zero length polynomial will return one of \code{numeric(length(var))} or
\code{complex(length(var))}, depending on the mode of the coefficients.

Simplifications are made where possible including not multiplying a variable by
\code{1}, not raising a variable to a power of \code{1}, not multiplying by a
variable raised to a power of \code{0}, pulling signs from trailing
coefficients, and exchanging complex coefficients for real coefficients. These
guidelines are not always followed for edge cases, these are for general use.
}
\value{
A class \code{"\linkS4class{call}"} object.
}
\examples{
## when supplied with exactly one coefficient, most simplifactions are ignored.
as.body(1+0i, strict = FALSE)


## zero valued coefficients are removed by 'reduce', this can be changed.
as.body(0)
as.body(0, zero.rm = FALSE)  # zero valued coefficient is retained.


as.body(Hermite(10), decreasing = TRUE)
as.body(Cyclotomic(105), decreasing = TRUE)


## 'var' can be a call to result in a more interesting result.
## the Chebyshev polynomials of the first kind are based on 'x = cos(theta)'
## and we can use 'var' to simulate this.
as.body(Chebyshev1(5), decreasing = TRUE,
    var = call("cos", as.symbol("theta")))


## we can also specify 'var' as an expresion to get the same thing.
as.body(Chebyshev1(5), decreasing = TRUE,
    var = expression(cos(theta)))


## you can combine this with 'as.function' like this:
as.function(reduce(Chebyshev1(5), decreasing = TRUE),
    xname = "theta", var = expression(cos(theta)))
}
