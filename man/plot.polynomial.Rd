\name{plot.polynomial}
\alias{plot.polynomial}
\alias{points.polynomial}
\alias{lines.polynomial}
\title{
Plot Method for Polynomials
}
\description{
\code{plot.polynomial}, a method for the
\code{\link[graphics:plot.default]{plot}} generic. It draws a curve
corresponding to a polynomial over the interval \code{[from, to]}.
}
\usage{
\S3method{plot}{polynomial}(x, y = NULL, to = NULL, from = y, n = 101, add = FALSE,
    type = "l", xname = "x", xlab = xname, ylab = NULL,
    log = NULL, xlim = NULL, ...)

\S3method{points}{polynomial}(x, y = NULL, to = NULL, from = y, n = 101, add = FALSE,
    type = "p", ...)

\S3method{lines}{polynomial}(x, y = NULL, to = NULL, from = y, n = 101, add = FALSE,
    type = "l", ...)
}
\arguments{
  \item{x}{
a polynomial, or an object which can be coerced by \code{\link{as.polynomial}}
to a polynomial.
  }
  \item{y}{
alias for \code{from} for compatibility with \code{plot}.
  }
  \item{from, to}{
the range over which the polynomial will be plotted.
  }
  \item{n}{
integer; the number of x values at which to evaluate.
  }
  \item{add}{
logical; if \code{TRUE} add to an already existing plot; if \code{NA} start a
new plot taking the defaults for the limits and log-scaling of the x-axis from
the previous plot. Taken as \code{FALSE} (with a warning if a different value is
supplied) if no graphics device is open.
  }
  \item{xlim}{
\code{NULL} or a numeric vector of length 2; if non-\code{NULL} it provides the
defaults for \code{c(from, to)} and, unless \code{add = TRUE}, selects the
x-limits of the plot.
  }
  \item{type}{
plot type: see \code{\link[graphics]{plot.default}}.
  }
  \item{xname}{
character string giving the name to be used for the x axis. Optionally a
\code{\link{call}} or \code{\link{expression}}, see \code{\link{as.body}}.
  }
  \item{xlab, ylab, log, \dots}{
labels and \code{\link[graphics:par]{graphical parameters}} can also be
specified as arguments. See \strong{Details} for the interpretation of the
default for \code{log}.
  }
}
\details{
The polynomial \code{x} is evaluated at \code{n} points equally spaced over the
range \code{[from, to]}. The points determined in this way are then plotted.

If either \code{from} or \code{to} is \code{NULL}, it defaults to the
corresponding element of \code{xlim} if that is not \code{NULL}.

When neither \code{from/to} nor \code{xlim} specifies both x-limits, the
x-limits are inferred from \code{add}. When \code{add = FALSE} the defaults are
\code{(0, 1)}. For \code{add = NA} and \code{add = TRUE} the defaults are taken
from the x-limits used for the previous plot.

The value of \code{log} is used both to specify the plot axes (unless \code{add
= TRUE}) and how \sQuote{equally spaced} is interpreted: if the x component
indicates log-scaling, the points at which the polynomial is plotted are equally
spaced on the log scale.

The default value of \code{log} is taken from the current plot when \code{add =
TRUE}, whereas if \code{add = NA} the x component is taken from the existing
plot (if any) and the y component defaults to linear. For \code{add = FALSE} the
default is \code{""}.
}
\note{
For \code{points.polynomial} and \code{lines.polynomial}, \code{add} is ignored
and is taken to be \code{TRUE}. If the \sQuote{active} device is the "null
device", an error is raised. See \code{\link[grDevices:dev]{dev.cur}}.
}
\value{
A list with components \code{x} and \code{y} of the points that were drawn is
returned invisibly.
}
\examples{\dontshow{opar <- graphics::par(mar = c(5, 4.4, 0.4, 0.4) + 0.1)}
x <- as.polynomial(Hermite(5))
plot(x, -3, 6, ylim = c(-10, 20))


## 'plot.polynomial' can be combined with 'reduce' for a different y label
## use 'add = NA' to use [from, to] from previous plot
plot(reduce(x, decreasing = TRUE), add = NA, ylim = c(-10, 20))


x <- as.polynomial(Chebyshev1(4))
plot(x, -1, 1, ylim = c(-1, 1))


## use 'xname' for a different x label (and polynomial variable)
plot(x, add = NA, ylim = c(-1, 1), xname = expression(cos(theta)))


## plot Chebyshev polynomials for 'n = 0, ..., 4'
plot(
    as.polynomial(Chebyshev1(0)),
    from = -1, to = 1, ylim = c(-1, 1),
    xname = expression(cos(theta)),
    ylab = expression(T[n](cos(theta)))
  )
lines (as.polynomial(Chebyshev1(1)), col = "red")
lines (as.polynomial(Chebyshev1(2)), col = "green3")
points(as.polynomial(Chebyshev1(3)), col = "blue")
points(as.polynomial(Chebyshev1(4)), col = "deeppink2")
\dontshow{graphics::par(opar)}}
