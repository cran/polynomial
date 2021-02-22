\name{reduce}
\alias{reduce}
\title{
Polynomial Simplification
}
\description{
Simplify a polynomial by removing unwanted coefficients, combining coefficients
corresponding to duplicate exponents, and sorting exponents.
}
\usage{
reduce(x, zero.rm = TRUE, na.rm = FALSE, finite = FALSE, dflt.exp = FALSE,
    fix.dup.exp = is.na(dflt.exp) || dflt.exp, decreasing = NA,
    ...)
}
\arguments{
  \item{x}{
a polynomial, or an object which can be coerced by \code{\link{as.polynomial}}
to a polynomial.
  }
  \item{zero.rm}{
logical. Should zero valued coefficients be removed?
  }
  \item{na.rm}{
logical. Should missing valued coefficients (including \code{NaN}) be removed?
  }
  \item{finite}{
logical. Should non-finite valued coefficients be removed?
  }
  \item{dflt.exp}{
logical. Should the coefficients be rearranged such that the exponents take on
their default values?
  }
  \item{fix.dup.exp}{
logical. Should coefficients corresponding to duplicate exponents be summed?
  }
  \item{decreasing}{
logical. Should the sort order of the exponents be increasing or decreasing?
  }
  \item{\dots}{
further arguments passed to \code{\link{as.polynomial}}.
  }
}
\details{
If \code{zero.rm} is \code{TRUE}, \code{0} valued coefficients and their
corresponding exponents will be removed from the polynomial.

If \code{na.rm} is \code{FALSE}, \code{NA} and \code{NaN} valued coefficients
will be retained, otherwise \code{NA} valued coefficients and their
corresponding exponents will be removed from the polynomial.

If \code{finite} is \code{TRUE}, non-finite valued coefficients and their
corresponding exponents will be removed from the polynomial, i.e., \code{finite = TRUE}
\emph{includes} \code{na.rm = TRUE}.

If \code{fix.dup.exp} is \code{TRUE}, coefficients corresponding to duplicated
exponents will be summed.

If \code{dflt.exp} is \code{TRUE}, the coefficients are rearranged such that the
exponents take on their default values, given by \code{seq(0L, along.with = x)}.
Note that this \emph{requires} the exponents to be unique (or \code{fix.dup.exp
= TRUE}).

If \code{decreasing} is \code{FALSE}, the polynomial is reordered such that the
exponents are increasing. If \code{decreasing} is \code{TRUE}, the polynomial is
reordered such that the exponents are decreasing. If \code{decreasing} is
\code{NA}, the polynomial is left as is (the default).
}
\value{
A polynomial object.
}
\examples{
x <- as.polynomial(
    c(1, NA, Inf, 6, 4, 3, 0, 1, 9),
    c(5,  7,   4, 4, 5, 8, 2, 3, 8)
)


reduce(x)                   # 0 valued coefficient is removed
reduce(x, zero.rm = FALSE)  # 0 valued coefficient is not removed
reduce(x, na.rm = TRUE)     # missing valued coefficient is removed
reduce(x, finite = TRUE)    # non-finite valued coefficients are removed


## 'dflt.exp = NA' means that 'reduce' will decide whether the coefficients
## should be rearranged such that the exponents take on their default values.
## In this case, it decides 'dflt.exp = TRUE'
reduce(x, dflt.exp = NA)
reduce(x, dflt.exp = TRUE)  # exponents take on their default values (forced)


## coefficients corresponding to duplicate exponents are summed
reduce(x, fix.dup.exp = TRUE)


reduce(x, decreasing = FALSE)  # increasing exponents
reduce(x, decreasing = TRUE)   # decreasing exponents
}
