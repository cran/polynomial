\name{CommonPolynomials}
\alias{CommonPolynomials}
\alias{Abel}
\alias{Bessel}
\alias{rBessel}
\alias{Chebyshev}
\alias{Chebyshev1}
\alias{Chebyshev2}
\alias{ChebyshevT}
\alias{ChebyshevU}
\alias{Cyclotomic}
\alias{Hermite}
\title{
Common Polynomials
}
\description{
Determines the polynomial coefficients of some common \code{n}-th order
polynomials.
}
\usage{
Abel(n, a = 1)

Bessel(n)
rBessel(n)

Chebyshev1(n)
Chebyshev2(n)
ChebyshevT(n)
ChebyshevU(n)

Cyclotomic(n)

Hermite(n)
}
\arguments{
  \item{n}{
a non-negative integer; the order of the polynomial. For \code{Cyclotomic},
\code{n} must be positive.
}
  \item{a}{
a real or complex number; most commonly a positive integer.
}
}
\details{
The \code{n}-th order Abel polynomial is defined as

\deqn{p_n(x) = x (x - a n)^{n - 1}}{
      p[n](x) = x (x - a n)^{n - 1}}

The \code{n}-th order Bessel polynomial (\eqn{y_n}{y[n]}) and reverse Bessel
polynomial (\eqn{\theta_n}{\theta[n]}) are defined as

\deqn{y_n(x) = \sum_{k=0}^n\frac{(n + k)!}{(n - k)!\,k!} (\frac{x}{2})^k}{
      y[n](x) = \sum (n + k)!/((n - k)! k!) (x/2)^k for k = 0, \ldots, n}

\deqn{\theta_n(x) = x^n y_n(1/x) = \sum_{k=0}^n\frac{(n + k)!}{(n - k)!\,k!}
\frac{x^{n - k}}{2^k}}{
      \theta[n](x) = x^n y[n](1/x) = \sum (n + k)!/((n - k)! k!) x^(n - k)/2^k
for k = 0, \ldots, n}

The \code{n}-th order Chebyshev polynomials of the first kind (\eqn{T_n}{T[n]})
and the second kind (\eqn{U_n}{U[n]}) are defined as

\deqn{T_n(cos(\theta)) = cos(n \theta)}{
      T[n](cos(\theta)) = cos(n \theta)}

\deqn{U_n(cos(\theta)) sin(\theta) = sin((n + 1) \theta)}{
      U[n](cos(\theta)) sin(\theta) = sin((n + 1) \theta)}

where upon defining \eqn{x = cos(\theta)} and using trig identities to replace
\eqn{cos(n \theta)} and \eqn{sin((n + 1) \theta)}, we find two polynomials in
terms of \code{x}.

The \code{n}-th order Cyclotomic polynomial is defined as

\deqn{\Phi_n(x) = \prod_{\stackrel{1 \le k \le n}{\gcd(k, n) = 1}} \left(x -
e^{2 i \pi \frac{k}{n}}\right)}{
      \Phi[n](x) = \prod (x - e^(2 i \pi k/n)) for k = 1, \ldots, n and k is
coprime with n}

The \code{n}-th order Hermite polynomial is defined as

\deqn{He(x) = (-1)^n exp(\frac{x^2}{2}) \frac{d^n}{dx^n} exp(-\frac{x^2}{2}}{
      He(x) = (-1)^n e^(x^2/2) d^n/dx^n e^(-x^2/2)}
}
\value{
For \code{Cyclotomic}, a numeric vector of length \code{sum(\link{coprime}(1:n,
n)) + 1}; the coefficients of the \code{n}-th order Cyclotomic polynomial.

For all others, a numeric vector of length \code{n + 1}; the coefficients of the
\code{n}-th order polynomial.
}
\examples{\dontshow{opar <- graphics::par(mar = c(5.1, 4.5, 2.1, 0.4))}
colours <- grDevices::palette.colors()[c("vermillion", "bluishgreen",
    "blue", "skyblue", "reddishpurple", "orange")]


## Abel polynomials
n <- 0:4
plot(
    xlim = c(-3, 6), ylim = c(-20, 30),
    panel.first = graphics::grid(col = "gray69"),
    x = NA_real_, y = NA_real_,
    xlab = "x", ylab = as.call(list(call("[", as.symbol("p"),
        as.symbol("n")), as.symbol("x"))),
    main = "Abel Polynomials"
)
for (i in seq_along(n)) {
    graphics::lines(as.polynomial(Abel(n[[i]])),
        col = colours[[i]], n = 1001, lwd = 2,
        xlim = graphics::par("usr")[1:2])
}
graphics::box()
graphics::legend(
    x        = "bottomright",
    legend   = as.expression(lapply(X = n, FUN = function(n1) {
        call("==", as.symbol("n"), n1)
    })),
    col      = colours[seq_along(n)],
    lwd      = 2,
    bty      = "n"
)
graphics::title(sub = ~"All with" ~ list(a == 1), adj = 1)


## Bessel polynomials
n <- 0:5
plot(
    xlim = c(-3, 6), ylim = c(-10, 20),
    panel.first = graphics::grid(col = "gray69"),
    x = NA_real_, y = NA_real_,
    xlab = "x", ylab = as.call(list(call("[", as.symbol("y"),
        as.symbol("n")), as.symbol("x"))),
    main = "Bessel Polynomials"
)
for (i in seq_along(n)) {
    graphics::lines(as.polynomial(Bessel(n[[i]])),
        col = colours[[i]], n = 1001, lwd = 2,
        xlim = graphics::par("usr")[1:2])
}
graphics::box()
graphics::legend(
    x        = "bottomright",
    legend   = as.expression(lapply(X = n, FUN = function(n1) {
        call("==", as.symbol("n"), n1)
    })),
    col      = colours[seq_along(n)],
    lwd      = 2,
    bty      = "n"
)


## Reverse Bessel polynomials
n <- 0:5
plot(
    xlim = c(-10, 10), ylim = c(-10, 20),
    panel.first = graphics::grid(col = "gray69"),
    x = NA_real_, y = NA_real_,
    xlab = "x", ylab = as.call(list(call("[", as.symbol("theta"),
        as.symbol("n")), as.symbol("x"))),
    main = "Reverse Bessel Polynomials"
)
for (i in seq_along(n)) {
    graphics::lines(as.polynomial(rBessel(n[[i]])),
        col = colours[[i]], n = 1001, lwd = 2,
        xlim = graphics::par("usr")[1:2])
}
graphics::box()
graphics::legend(
    x        = "bottomright",
    legend   = as.expression(lapply(X = n, FUN = function(n1) {
        call("==", as.symbol("n"), n1)
    })),
    col      = colours[seq_along(n)],
    lwd      = 2,
    bty      = "n"
)


## Chebyshev polynomials of the first kind
n <- 0:4
plot(
    xlim = c(-1, 1), ylim = c(-1, 1),
    panel.first = graphics::grid(col = "gray69"),
    x = NA_real_, y = NA_real_,
    xlab = "x", ylab = as.call(list(call("[", as.symbol("T"),
        as.symbol("n")), as.symbol("x"))),
    main = "Chebyshev polynomials of the first kind"
)
for (i in seq_along(n)) {
    graphics::lines(as.polynomial(Chebyshev1(n[[i]])),
        col = colours[[i]], n = 1001, lwd = 2,
        xlim = graphics::par("usr")[1:2]
    )
}
graphics::box()
graphics::legend(
    x        = "bottomright",
    legend   = as.expression(lapply(X = n, FUN = function(n1) {
        call("==", as.symbol("n"), n1)
    })),
    col      = colours[seq_along(n)],
    lwd      = 2,
    bty      = "n"
)


## Chebyshev polynomials of the second kind
n <- 0:4
plot(
    xlim = c(-1, 1), ylim = c(-4, 5),
    panel.first = graphics::grid(col = "gray69"),
    x = NA_real_, y = NA_real_,
    xlab = "x", ylab = as.call(list(call("[", as.symbol("U"),
        as.symbol("n")), as.symbol("x"))),
    main = "Chebyshev polynomials of the second kind"
)
for (i in seq_along(n)) {
    graphics::lines(as.polynomial(Chebyshev2(n[[i]])),
        col = colours[[i]], n = 1001, lwd = 2,
        xlim = graphics::par("usr")[1:2]
    )
}
graphics::box()
graphics::legend(
    x        = "bottomright",
    legend   = as.expression(lapply(X = n, FUN = function(n1) {
        call("==", as.symbol("n"), n1)
    })),
    col      = colours[seq_along(n)],
    lwd      = 2,
    bty      = "n"
)


## Cyclotomic polynomials
n <- 1:5
plot(
    xlim = c(-3, 3), ylim = c(-10, 10),
    panel.first = graphics::grid(col = "gray69"),
    x = NA_real_, y = NA_real_,
    xlab = "x", ylab = as.call(list(call("[", as.symbol("Phi"),
        as.symbol("n")), as.symbol("x"))),
    main = "Cyclotomic Polynomials"
)
for (i in seq_along(n)) {
    graphics::lines(as.polynomial(Cyclotomic(n[[i]])),
        col = colours[[i]], n = 1001, lwd = 2,
        xlim = graphics::par("usr")[1:2]
    )
}
graphics::box()
graphics::legend(
    x        = "bottomright",
    legend   = as.expression(lapply(X = n, FUN = function(n1) {
        call("==", as.symbol("n"), n1)
    })),
    col      = colours[seq_along(n)],
    lwd      = 2,
    bty      = "n"
)


## Hermite polynomials
n <- 0:5
plot(
    xlim = c(-3, 6), ylim = c(-10, 20),
    panel.first = graphics::grid(col = "gray69"),
    x = NA_real_, y = NA_real_,
    xlab = "x", ylab = as.call(list(call("[", as.symbol("H"),
        as.symbol("n")), as.symbol("x"))),
    main = "Hermite Polynomials"
)
for (i in seq_along(n)) {
    graphics::lines(as.polynomial(Hermite(n[[i]])),
        col = colours[[i]], n = 1001, lwd = 2,
        xlim = graphics::par("usr")[1:2]
    )
}
graphics::box()
graphics::legend(
    x        = "bottomright",
    legend   = as.expression(lapply(X = n, FUN = function(n1) {
        call("==", as.symbol("n"), n1)
    })),
    col      = colours[seq_along(n)],
    lwd      = 2,
    bty      = "n"
)
\dontshow{graphics::par(opar)}}
