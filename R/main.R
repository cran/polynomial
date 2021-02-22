methods::setClass(
    Class     = "polynomial",
    contains  = c("numbers", "oldClass"),
    slots     = c(exponents = "integer"),
    prototype = methods::prototype(.S3Class = "polynomial")
)


methods::setMethod(
    f = "coerce",
    signature = c(from = "ANY", to = "polynomial"),
    definition = function (from, to, strict = TRUE)
{
    value <- as.polynomial(from)
    if (strict)
        attributes(value) <- attributes(value)[c("names", "class", "exponents")]
    value
})





polynomial <- function (length = 0L)
.polynomial(numeric(length))


as.polynomial <- function (x, exponents = NULL, ...)
{
    if (missing(x))
        return(polynomial(0))
    if (inherits(x, "polynomial"))
        return(x)
    x <- structure(as.numbers(x = x, ...), dim = dim(x), dimnames = dimnames(x),
        names = names(x), class = "polynomial")
    exponents(x) <- exponents
    x
}


is.polynomial <- function (x)
inherits(x, "polynomial")


.polynomial <- function (xx, exponents = integer(0), cl = "polynomial")
{
    class(xx) <- cl
    attr(xx, "exponents") <- exponents
    xx
}





is.numeric.polynomial <- function (x)
FALSE


is.complex.polynomial <- function (x)
FALSE





Abel <- function (n, a = 1)
{
    n <- as.scalar.integer(n)
    if (is.na(n) || n < 0L)
        stop("'n' must be a non-negative integer")
    if (n == 0L) {
        1
    }
    else if (n == 1L) {
        c(0, 1)
    }
    else {
        a <- as.scalar.number(a)
        k <- (n - 1L):0
        c(0, choose(n - 1, k) * (-a * n)^k)
    }
}


Bessel <- function (n)
{
    n <- as.scalar.integer(n)
    if (is.na(n) || n < 0L)
        stop("'n' must be a non-negative integer")
    vapply(X = 0:n, FUN = function(k) {
        prod(seq.int(to = n + k, length.out = 2L * k)/c(seq_len(k), rep(2L, k)))
    }, FUN.VALUE = 0)
}


rBessel <- function (n)
rev(Bessel(n))


Hermite <- function (n)
.Call(C_hermite, n)





.Chebyshev <- function (x)
{
    value <- function(n) NULL
    body(value) <- substitute({
        n <- as.scalar.integer(n)
        if (is.na(n) || n < 0L)
            stop("'n' must be a non-negative integer")
        Ts <- vector("list", length = n + 1L)
        Ts[1:2] <- list(1, c(0, x))
        fun <- function(n) {
            if (!is.null(t <- Ts[[n]]))
                t
            else (Ts[[n]] <<- c(0, 2 * fun(n - 1L)) - c(fun(n - 2L), 0, 0))
        }
        fun(n + 1L)
    }, list(x = x))
    environment(value) <- parent.frame()
    return(value)
}


Chebyshev1 <- ChebyshevT <- .Chebyshev(1)


Chebyshev2 <- ChebyshevU <- .Chebyshev(2)


remove(.Chebyshev)





gcd <- function (..., na.rm = FALSE)
.External(C_gcd, na.rm, ...)


pgcd <- function (..., na.rm = FALSE)
.External(C_pgcd, na.rm, ...)


lcm <- function (..., na.rm = FALSE)
.External(C_lcm, na.rm, ...)


plcm <- function (..., na.rm = FALSE)
.External(C_plcm, na.rm, ...)


coprime <- function (...)
gcd(...) == 1L


pcoprime <- function (...)
pgcd(...) == 1L


Cyclotomic <- function (n)
{
    n <- as.scalar.integer(n)
    if (is.na(n) || n <= 0L)
        stop("'n' must be a positive integer")
    k <- seq_len(n)
    k <- lapply(X = -exp(2i * pi * k[pcoprime(k, n)]/n), FUN = "c", 1)
    value <- k[[1L]]
    fun <- function(x, y) {
        value <- numeric(length(x) + 1L)
        value[-length(value)] <- value[-length(value)] + x * y[1L]
        value[-1L] <- value[-1L] + x * y[2L]
        value
    }
    for (y in k[-1L]) value <- fun(value, y)
    round(Re(value))
}





degree <- function (x)
{
    if (!length(k <- attr(x, "exponents")))
        max(0L, length(x) - 1L)
    else max(0L, exponents(x), na.rm = TRUE)
}


exponents <- function (x)
{
    if (length(value <- attr(x, "exponents")))
        value
    else seq.int(0L, along.with = x)
}


"exponents<-" <- function (x, value)
{
    if (!is.polynomial(x))
        x <- as.polynomial(x)
    if (is.null(value)) {
        attr(x, "exponents") <- integer(0)
        return(x)
    }
    value <- as.integer(value)
    if (length(value) != length(x))
        stop("invalid 'exponents' length")
    else if (anyNA(value))
        stop("missing values in 'exponents' are not allowed")
    else if (any(value < 0L))
        stop("negative values in 'exponents' are not allowed")
    attr(x, "exponents") <- value
    x
}





as.body <- function (x, var = "x", ...)
{
    x <- reduce(x = x, ...)
    if (is.expression(var))
        var <- aslength1(var)[[1L]]
    if (!is.call(var))
        var <- as.symbol(var)
    if (!length(x))
        return(substitute(fun(length(var)), list(fun = if (mode(x) ==
            "complex") as.symbol("complex") else as.symbol("numeric"),
            var = var)))
    dim(x) <- dimnames(x) <- names(x) <- NULL
    a <- coef(x)
    k <- exponents(x)
    if (length(x) == 1L) {
        if (k != 1L)
            var <- call("^", var, Re(k))
        return(if (is.na(a) || is.complex(x))
            call("*", a, var)
        else if (a == 1)
            var
        else if (a == -1)
            call("-", var)
        else call("*", a, var))
    }
    vars <- lapply(X = Re(k), FUN = function(y) {
        if (y == 1)
            var
        else if (!y)
            NULL
        else call("^", var, y)
    })
    if (all(k == 0L))
        vars[[length(vars)]] <- call("^", var, 0)
    if (is.numeric(a)) {
        signs <- sign(a)[-1L]
        signs[is.na(signs)] <- 1
        a <- c(1, signs) * a
    }
    else {
        signs <- sign(Re(a))[-1L]
        i <- is.na(signs) | signs == 0
        signs[i] <- sign(Im(a))[-1L][i]
        signs[is.na(signs)] <- 1
        a <- complex(real = c(1, signs) * Re(a),
            imaginary = c(1, signs) * Im(a))
        a <- lapply(a, function(x) {
            if (!is.na(x) && !Im(x))
                Re(x)
            else x
        })
        if (!any(vapply(a, "is.complex", NA)))
            a[[length(a)]] <- a[[length(a)]] + 0i
    }
    vars <- .mapply(function(e1, e2) {
        if (is.null(e2))
            e1
        else if (!is.na(e1) && e1 == 1)
            e2
        else call("*", e1, e2)
    }, list(a, vars), NULL)
    i <- signs == -1
    signs[i] <- "-"
    signs[!i] <- "+"
    value <- vars[[1L]]
    for (i in seq_along(signs)) {
        value <- call(signs[[i]], value, vars[[i + 1L]])
    }
    return(value)
}


integral <- function (x, lower, upper, C, ...)
{
    if (!is.polynomial(x))
        x <- as.polynomial(x)
    a <- coef(x)
    k <- exponents(x) + 1L
    if (missing(C))
        return(diff(predict(object = .polynomial(a/k, k), newdata = c(lower, upper), ...)))
    C <- as.scalar.number(C)
    if (!length(attr(x, "exponents")))
        .polynomial(c(C, a/k))
    else .polynomial(c(C, a/k), c(0L, k))
}





removeInvalidExponents <- function (x)
{
    if (!length(k <- attr(x, "exponents")))
        x
    else if (any(i <- is.na(k) | k < 0L))
        x[!i]
    else x
}


sortExponents <- function (x, decreasing = FALSE)
{
    x[order(exponents(x), decreasing = decreasing)]
}


fixDuplicateExponents <- function (x)
{
    if (!length(k <- attr(x, "exponents"))) {
        x
    }
    else if (anyDuplicated(k)) {
        dim(x) <- dimnames(x) <- names(x) <- NULL
        a <- coef(x)
        dk <- duplicated(k)
        u <- unique(k[dk])
        i <- k %in% u
        a[match(u, k)] <- vapply(X = split(a[i], factor(k[i], u)),
            FUN = "sum", FUN.VALUE = vector(mode = mode(a), length = 1L))
        .polynomial(a[!dk], k[!dk])
    }
    else x
}


withDefaultExponents <- function (x)
{
    if (!length(k <- attr(x, "exponents")))
        return(x)
    value <- polynomial(max(k, na.rm = TRUE) + 1L)
    value[k + 1L] <- coef(x)
    return(value)
}


removeUnwantedCoefficients <- function (x, zero.rm = TRUE, na.rm = FALSE, finite = FALSE)
{
    i <- if (zero.rm) {
        if (finite)
            is.finite(x) & x
        else if (na.rm)
            !is.na(x) & x
        else is.na(x) | x
    }
    else if (finite) {
        is.finite(x)
    }
    else if (na.rm) {
        !is.na(x)
    }
    else return(x)
    if (all(i))
        x
    else x[i]
}


maybeDefaultExponents <- function (x)
{
    if (length(x) == 0L)
        x
    else if ((degree(x) + 1L)/2 <= length(x))
        withDefaultExponents(x)
    else x
}


reduce <- function (x, zero.rm = TRUE, na.rm = FALSE, finite = FALSE, dflt.exp = FALSE,
    fix.dup.exp = is.na(dflt.exp) || dflt.exp, decreasing = NA,
    ...)
{
    if (!is.polynomial(x))
        x <- as.polynomial(x = x, ...)
    x <- removeInvalidExponents(x)
    x <- removeUnwantedCoefficients(x, zero.rm = zero.rm, na.rm = na.rm,
        finite = finite)
    dflt.exp <- as.scalar.logical(dflt.exp)
    if (fix.dup.exp) {
        x <- fixDuplicateExponents(x)
        x <- removeUnwantedCoefficients(x, zero.rm = zero.rm,
            na.rm = na.rm, finite = finite)
        if (is.na(dflt.exp))
            x <- maybeDefaultExponents(x)
        else if (dflt.exp)
            x <- withDefaultExponents(x)
    }
    else if (!anyDuplicated(exponents(x))) {
        if (is.na(dflt.exp))
            x <- maybeDefaultExponents(x)
        else if (dflt.exp)
            x <- withDefaultExponents(x)
    }
    decreasing <- as.scalar.logical(decreasing)
    if (is.na(decreasing))
        x
    else sortExponents(x, decreasing = decreasing)
}





as.function.polynomial <- function (x, envir = parent.frame(), xname = "x", var = xname,
    ...)
{
    xname <- as.symbol(xname)
    value <- function(x) NULL
    names(formals(value)) <- as.character(xname)
    body(value) <- as.body(x = x, var = var, ...)
    environment(value) <- envir
    return(value)
}


as.list.polynomial <- function (x, ...)
{
    value <- .mapply(.polynomial, list(coef(x), exponents(x)), NULL)
    names(value) <- names(x)
    return(value)
}


print.polynomial <- function (x, ...)
{
    if (length(x) == 0L) {
        if (length(d <- dim(x)) > 1L)
            cat(sprintf("<%s polynomial>\n", paste0(d, collapse = " x ")))
        else if (!is.null(names(x)))
            cat("named polynomial(0)\n")
        else cat("polynomial(0)\n")
    }
    else {
        xx <- coef(x)
        keepAttrs <- setdiff(names(attributes(x)), c("exponents",
            "class"))
        attributes(xx)[keepAttrs] <- attributes(x)[keepAttrs]
        print(xx)
        cat("Exponents:\n")
        print(exponents(x))
    }
    invisible(x)
}


"-.polynomial" <- function (e1, e2)
{
    if (!is.polynomial(e1))
        e1 <- as.polynomial(e1)
    if (nargs() == 1L)
        return(.polynomial(-coef(e1), attr(e1, "exponents")))
    if (!is.polynomial(e2))
        e2 <- as.polynomial(e2)
    if (!length(e2)) {
        e1
    }
    else if (!length(e1)) {
        .polynomial(-coef(e2), attr(e2, "exponents"))
    }
    else if (!length(attr(e1, "exponents")) && !length(attr(e2, "exponents"))) {
        e1 <- coef(e1)
        e2 <- coef(e2)
        if (length(e1) > length(e2)) {
            names(e2) <- NULL
            e2 <- c(e2, rep(0, length(e1) - length(e2)))
        }
        else if (length(e1) < length(e2)) {
            names(e1) <- NULL
            e1 <- c(e1, rep(0, length(e2) - length(e1)))
        }
        .polynomial(e1 - e2)
    }
    else reduce(.polynomial(c(coef(e1), -coef(e2)), c(exponents(e1), exponents(e2))), zero.rm = FALSE, fix.dup.exp = TRUE)
}


"*.polynomial" <- function (e1, e2)
{
    if (nargs() == 1L)
        stop("invalid unary operator")
    if (!is.polynomial(e1))
        e1 <- as.polynomial(e1)
    if (!is.polynomial(e2))
        e2 <- as.polynomial(e2)
    if (!length(e1) || !length(e2)) {
        polynomial(0)
    }
    else if (!length(attr(e1, "exponents")) && !length(attr(e2, "exponents"))) {
        e1 <- coef(e1)
        e2 <- coef(e2)
        if (length(e1) != 1L && length(e2) != 1L) {
            value <- numeric(length(e1) + length(e2) - 1L)
            n <- seq.int(0L, along.with = e1)
            for (i in seq_along(e2)) {
                value[n + i] <- value[n + i] + e1 * e2[[i]]
            }
            .polynomial(value)
        }
        else .polynomial(e1 * e2)
    }
    else reduce(.polynomial(do.call("c", c(lapply(coef(e1), "*", coef(e2)), list(use.names = FALSE))),
        do.call("c", c(lapply(exponents(e1), "+", exponents(e2)), list(use.names = FALSE)))),
        zero.rm = FALSE, fix.dup.exp = TRUE)
}


"/.polynomial" <- function (e1, e2)
{
    if (nargs() == 1L)
        stop("invalid unary operator")
    if (!is.polynomial(e1))
        e1 <- as.polynomial(e1)
    if (!is.polynomial(e2))
        e2 <- as.polynomial(e2)
    if (!length(e1) || !length(e2))
        return(polynomial(0))
    if (!length(attr(e2, "exponents")) && length(e2) == 1L)
        return(.polynomial(coef(e1)/coef(e2), attr(e1, "exponents")))
    e1 <- rev(coef(reduce(e1, dflt.exp = TRUE)))
    e2 <- rev(coef(reduce(e2, dflt.exp = TRUE)))
    if (!length(e1) || !length(e2))
        return(polynomial(0))
    len <- max(0L, length(e1) - length(e2) + 1L)
    value <- vector(mode(e1[0L] + e2[0L]), len)
    for (i in seq.int(to = 1L, by = -1L, length.out = len)) {
        value[[i]] <- e1[[1L]]/e2[[1L]]
        e1 <- e1[-1L] - c(value[[i]] * e2[-1L], integer(i - 1L))
    }
    if (anyNA(e1) || any(e1 != 0))
        warning("non-zero remainder in polynomial division")
    .polynomial(value)
}


"[.polynomial" <- function (x, ...)
{
    k <- exponents(x)
    a <- x <- coef(x)
    suppressWarnings(storage.mode(a) <- "integer")
    a[] <- seq_along(a)
    .polynomial(x[...], k[as.vector(a[...])])
}


"[[.polynomial" <- function (x, ...)
{
    k <- exponents(x)
    a <- x <- coef(x)
    a[] <- seq_along(a)
    .polynomial(x[[...]], k[[as.vector(a[[...]])]])
}


"+.polynomial" <- function (e1, e2)
{
    if (!is.polynomial(e1))
        e1 <- as.polynomial(e1)
    if (nargs() == 1L)
        return(.polynomial(+coef(e1), attr(e1, "exponents")))
    if (!is.polynomial(e2))
        e2 <- as.polynomial(e2)
    if (!length(e2)) {
        e1
    }
    else if (!length(e1)) {
        .polynomial(+coef(e2), attr(e2, "exponents"))
    }
    else if (!length(attr(e1, "exponents")) && !length(attr(e2, "exponents"))) {
        e1 <- coef(e1)
        e2 <- coef(e2)
        if (length(e1) > length(e2)) {
            names(e2) <- NULL
            e2 <- c(e2, rep(0, length(e1) - length(e2)))
        }
        else if (length(e1) < length(e2)) {
            names(e1) <- NULL
            e1 <- c(e1, rep(0, length(e2) - length(e1)))
        }
        .polynomial(e1 + e2)
    }
    else reduce(.polynomial(c(coef(e1), coef(e2)), c(exponents(e1), exponents(e2))), zero.rm = FALSE, fix.dup.exp = TRUE)
}


"^.polynomial" <- function (e1, e2)
{
    if (nargs() == 1L)
        stop("invalid unary operator")
    if (!is.polynomial(e1))
        e1 <- as.polynomial(e1)
    if (is.polynomial(e2)) {
        if (length(attr(e2, "exponents"))) {
            tmp <- coef(e2)[!exponents(e2) %in% 0L]
            if (anyNA(tmp) || any(tmp != 0))
                warning("non-constant parts of 'e2' discarded in coercion")
            e2 <- as.integer(sum(coef(e2)[exponents(e2) %in% 0L]))
        }
        else if (length(e2)) {
            e2 <- as.integer(coef(e2)[[1L]])
        }
        else return(polynomial(0))
    }
    else if (length(e2)) {
        e2 <- as.scalar.integer(e2)
    }
    else return(polynomial(0))
    if (is.na(e2)) {
        e1[] <- NA
        return(e1)
    }
    if (e2 < 0L) {
        warning("'e2' should be a non-negative integer")
        return(polynomial(0))
    }
    if (e2 == 0L)
        return(as.polynomial(1))
    if (e2 == 1L)
        return(e1)
    value <- e1
    for (i in logical(e2 - 1L)) value <- value * e1
    return(value)
}





plot.polynomial <- function (x, y = NULL, to = NULL, from = y, n = 101, add = FALSE,
    type = "l", xname = "x", xlab = xname, ylab = NULL, log = NULL,
    xlim = NULL, ...)
{
    if (dev.cur() == 1L && !isFALSE(add)) {
        warning("'add' will be ignored as there is no existing plot")
        add <- FALSE
    }
    addF <- isFALSE(add)
    if (is.null(ylab) && !isTRUE(add))
        try(ylab <- as.body(signif(x, 7L), var = xname, finite = TRUE))
    if (is.null(from) || is.null(to)) {
        xl <- if (!is.null(xlim))
            xlim
        else if (!addF) {
            pu <- par("usr")[1:2]
            if (par("xaxs") == "r")
                pu <- extendrange(pu, f = -1/27)
            if (par("xlog"))
                10^pu
            else pu
        }
        else c(0, 1)
        if (is.null(from))
            from <- xl[1L]
        if (is.null(to))
            to <- xl[2L]
    }
    lg <- if (length(log))
        log
    else if (!addF && par("xlog"))
        "x"
    else ""
    y <- as.polynomial(x)
    if (grepl("x", lg, fixed = TRUE)) {
        if (from <= 0 || to <= 0)
            stop("'from' and 'to' must be > 0 with log=\"x\"")
        x <- exp(seq.int(log(from), log(to), length.out = n))
    }
    else x <- seq.int(from, to, length.out = n)
    y <- predict(y, x, finite = TRUE)
    if (isTRUE(add))
        lines(x = x, y = y, type = type, ...)
    else plot(x = x, y = y, type = type, xlab = xlab, ylab = ylab,
        xlim = xlim, log = lg, ...)
    invisible(list(x = x, y = y))
}


points.polynomial <- function (x, y = NULL, to = NULL, from = y, n = 101, add = FALSE,
    type = "p", ...)
{
    if (dev.cur() == 1L)
        stop("plot.new has not been called yet")
    plot.polynomial(x = x, y = y, to = to, from = from, n = n,
        add = TRUE, type = type, ...)
}


lines.polynomial <- function (x, y = NULL, to = NULL, from = y, n = 101, add = FALSE,
    type = "l", ...)
{
    if (dev.cur() == 1L)
        stop("plot.new has not been called yet")
    plot.polynomial(x = x, y = y, to = to, from = from, n = n,
        add = TRUE, type = type, ...)
}





coef.polynomial <- function (object, ...)
object@.Data


"coef<-" <- "coefficients<-" <- function (object, ..., value)
UseMethod("coef<-")


"coef<-.default" <- function (object, ..., value)
{
    if (length(object) != length(value))
        stop(sprintf("invalid 'value', length of value [%d] should be equal to length of object [%d]",
             length(value), length(object)))
    if (!is.polynomial(object))
        object <- as.polynomial(object)
    value <- as.numbers(x = value, ...)
    attributes(value) <- c(attributes(object), class = "polynomial")
    return(value)
}


deriv.polynomial <- function (expr, ...)
{
    k <- exponents(expr)
    reduce(.polynomial(k * coef(expr), k - 1L), zero.rm = FALSE)
}


predict.polynomial <- function (object, newdata = seq.int(0, 1, 0.01), ...)
{
    object <- reduce(x = object, ...)
    value <- vector(mode = mode(object), length = length(newdata))
    for (e2 in .mapply(function(a, k) {
        a * newdata^k
    }, list(coef(object), exponents(object)), NULL)) {
        value <- value + e2
    }
    return(value)
}
