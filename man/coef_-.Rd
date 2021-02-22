\name{coef<-}
\alias{coef<-}
\alias{coef<-.default}
\alias{coefficients<-}
\title{
Replace Coefficients
}
\description{
\code{coef<-} is a generic function which replaces the coefficients of a
polynomial. \code{coefficients<-} is an \emph{alias} for it.
}
\usage{
coef(object, ...) <- value
coefficients(object, ...) <- value

\S3method{coef}{default}(object, \dots) <- value
}
\arguments{
  \item{object}{
any object.
  }
  \item{\dots}{
further arguments passed to or from other methods. The default method passes
these arguments to \code{\link{as.numbers}}.
  }
  \item{value}{
a numeric or complex vector, or an object which can be coerced by
\code{\link{as.numbers}} to a numeric or complex vector.
  }
}
\details{
All attributes of \code{object} are retained except for the \code{class}
attribute which is set to \code{"polynomial"}.
}
\value{
For \code{coef<-}, the updated object. (Note that the value of
\code{coef(object, ...) <- value} is that of the assignment, \code{value}, not
the return value from the left-hand side.)
}
\examples{
x <- as.polynomial(Hermite(5))
names(x) <- c("a", "b", "c", "d", "e", "f")
print(x)

coef(x) <- 1:6
print(x)  # names are preserved

coef(x)[2] <- 0
print(x)
}
