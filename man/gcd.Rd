\name{gcd}
\alias{gcd}
\alias{pgcd}
\alias{lcm}
\alias{plcm}
\alias{coprime}
\alias{pcoprime}
\title{
Greatest Common Divisor and Least Common Multiple
}
\description{
\code{gcd} and \code{lcm} return the greatest common divisor and least common
multiple of all the values present in its arguments.

\code{pgcd} and \code{plcm} return the \strong{p}arallel greatest common divisor
and \strong{p}arallel least common multiple of the input values. They take any
number of vectors as arguments, recycle them to common length and return a
single vector.
}
\usage{
gcd(..., na.rm = FALSE)
pgcd(..., na.rm = FALSE)

lcm(..., na.rm = FALSE)
plcm(..., na.rm = FALSE)

coprime(...)
pcoprime(...)
}
\arguments{
  \item{\dots}{
integer arguments. \code{NULL} is accepted as equivalent to \code{integer(0)}.
  }
  \item{na.rm}{
a logical indicating whether missing values should be removed.
  }
}
\details{
The greatest common divisor \code{d} of two integers \code{x} and \code{y} (at
least one of which is non-zero) is the greatest number for which \code{x/d} and
\code{y/d} are both integers. When both \code{x} and \code{y} are zero, the
greatest common divisor is also zero.

The least common multiple \code{m} of two non-zero integers \code{x} and
\code{y} is the least number for which \code{m/x} and \code{m/y} are both
integers. When one or both of \code{x} and \code{y} are zero, the least common
multiple is also zero.

If \code{x} or \code{y} is \code{NA}, the greatest common divisor and least
common multiple are also \code{NA}.

It is guaranteed that \code{d * m == abs(x * y)}

\code{x} is coprime with \code{y} if their greatest common divisor is \code{1}.
}
\value{
For \code{gcd} and \code{lcm} an integer vector of length 1. For \code{coprime}
a logical vector of length 1.

For \code{pgcd} and \code{plcm} an integer vector. For \code{pcoprime} a logical
vector. If any of the input values is a zero-length vector the result has length
zero. Otherwise, the result has length equal to the length of the longest
vector. The rules for determining the attributes of the result are rather
complicated. Attributes are only copied from input values whose lengths are
equal to the length of the result. If any such input values have a
\code{\link{dim}} attribute, the first \code{dim} attribute is copied to the
result. \code{\link{dimnames}} are copied in a similar manner (but only
\emph{after} the result has a \code{dim} attribute). If any such input values
have a \emph{conformable} \code{dimnames} attribute, the first conformable
\code{dimnames} attribute is copied to the result. If a \code{dim} attribute has
\emph{not} been assigned to the result, then finally \code{\link{names}} are
copied in a similar manner. If any such input values have a \code{names}
attribute, the first \code{names} attribute is copied to the result. A result
can have a \code{dim} attribute, a \code{names} attribute, neither, but cannot
have both. \code{dim} has priority over \code{names} (similar to
\code{\link{Arithmetic}} operators).
}
\examples{
y <- 105
x <- -y:y
d <- pgcd(x, y)
m <- plcm(x, y)


graphics::plot(
    x = x, y = d,
    ylab = substitute("Greatest Common Divisor of x and " ~ y, list(y = y))
)
graphics::plot(
    x = x, y = m,
    ylab = substitute("Least Common Multiple of x and " ~ y, list(y = y))
)
d * m == abs(x * y)


pcoprime(x, 105)
}
