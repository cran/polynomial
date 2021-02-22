\name{polynomial-class}
\Rdversion{1.1}
\docType{class}
\alias{polynomial-class}
\alias{coerce,ANY,polynomial-method}
\alias{degree}
\alias{integral}

\title{Class \code{"polynomial"}}
\description{
A class to hold the coefficients and exponents of a polynomial.
}
\section{Objects from the Class}{
Objects can be created by calls to the function \code{\link{as.polynomial}}.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{
Object of class \code{"\linkS4class{numbers}"}: numeric or complex coefficients.
Despite this, the functions \code{\link{is.numeric}} and
\code{\link{is.complex}} will \emph{always} return \code{FALSE} because
arithmetic on them does not make sense (inconsistent with the base arithmetic
operators).
    }
    \item{\code{exponents}:}{
Object of class \code{"integer"}: exponents must be not \code{NA} and non-negative.
    }
  }
}
\section{Extends}{
Class \code{"\linkS4class{numbers}"}, from data part.
}
\section{Methods}{
  \describe{
    \item{-}{
Negate a polynomial or subtract two polynomials.
    }
    \item{*}{
Multiply two polynomials.
    }
    \item{/}{
Divide two polynomials. A warning is issued for non-zero remainder.
    }
    \item{[}{
Subset a polynomial. Matrix and array indexing is allowed. Special care has been
taken to ensure that in such cases, the exponents are subsetted correctly.
    }
    \item{[[}{
Subset a polynomial. Matrix and array indexing is allowed. Special care has been
taken to ensure that in such cases, the exponents are subsetted correctly.
    }
    \item{^}{
Raise a polynomial to a non-negative integer power. As with ordinary \code{^}, a
power of \code{0} returns \code{1} (as a polynomial, however), and a power of
\code{1} returns the unmodified polynomial. A zero length power returns a zero
length polynomial. A negative or missing valued power returns the same
polynomial with all missing coefficients.
    }
    \item{+}{
Unary \code{+} or add two polynomials.
    }
    \item{as.function}{
Coerce a polynomial to a univariate function. The name of the variable is
determined by \code{xname} and the body of the function is determined by
\code{\link{as.body}}.
    }
    \item{as.list}{
Coerce a polynomial to a list. Names are preserved.
    }
    \item{coef, exponents}{
Select the coefficients or exponents part only of a polynomial.
    }
    \item{coerce}{
Coerce any object to a polynomial.
    }
    \item{degree}{
Calculate the degree of a polynomial, the maximum coefficient.
    }
    \item{deriv, integral}{
Calculate the derivative or integral of a polynomial.
    }
    \item{plot, points, lines}{
Plot a polynomial. See \code{\link{plot.polynomial}}.
    }
    \item{predict}{
Calculate the value of a polynomial at a set of points.
    }
  }
}
\seealso{
\code{\link{CommonPolynomials}}, \code{\link{reduce}}, \code{\link{as.body}}.
}
\keyword{classes}
