\name{polynomial}
\alias{polynomial}
\alias{as.polynomial}
\alias{is.polynomial}
\title{
Polynomials
}
\description{
Functions to create, coerce and check for polynomials.
}
\usage{
polynomial(length = 0)
as.polynomial(x, exponents = NULL, ...)
is.polynomial(x)
}
\arguments{
  \item{length}{
A non-negative integer specifying the desired length. Double values will be
coerced to integer: supplying an argument of length other than one is an error.
  }
  \item{x}{
object to be coerced or tested.
  }
  \item{exponents}{
\code{NULL} or an integer vector of the same length as \code{x}.
  }
  \item{\dots}{
further arguments passed to \code{\link[essentials:numbers]{as.numbers}}.
  }
}
\details{
\code{polynomial} creates a double-precision vector of the specified length with
each coefficient equal to 0 and exponents equal to
\code{seq(0L, length.out = length)}.

\code{as.polynomial} returns \code{x} if \code{x} is a polynomial. If
not, it attempts to coerce its argument to be of double or complex type: most
attributes will be stripped (excluding \code{\link{dim}}, \code{\link{dimnames}}
and \code{\link{names}}). Finally, the \code{"class"} and \code{"exponents"}
attibutes are populated.
}
\value{
\code{polynomial} and \code{as.polynomial} both return polynomial objects.

\code{is.polynomial} returns \code{TRUE} or \code{FALSE} depending on whether
its argument is a polynomial or not.
}
\examples{
as.polynomial(c(0, 4, 2, -1))


## 'as.polynomial' reduces complex numbers to real numbers if appropriate.
## You can use 'strict = FALSE' to stop this conversion.
as.polynomial(complex(real = c(0, 4, 2, -1)))
as.polynomial(complex(real = c(0, 4, 2, -1)), strict = FALSE)
}
