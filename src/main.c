#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>





// replace an NA with 0 or 1
#define filter_na0(X) ((X == NA_INTEGER) ? 0 : X)
#define filter_na1(X) ((X == NA_INTEGER) ? 1 : X)





SEXP hermite(SEXP n)
{
    int degree = asInteger(n);
    SEXP value;

    if (xlength(n) != 1)
        warning("first element used of '%s' argument", "n");

    if (degree == NA_INTEGER || degree < 0)
        error("'%s' must be a non-negative integer", "n");
    if (degree > 65535)
        error("'%s' too large", "n");

    value = PROTECT(allocVector(REALSXP, degree + 1));
    double tmp[degree], *rvalue = REAL(value);

    // first we populate rvalue with the appropriate coefficients:
    // 1 at the first position and 0 in all the rest
    memset(rvalue, 0.0, (degree + 1) * sizeof(double));
    rvalue[0] = 1.0;

    // each time we loop through, we calculate the next set of hermite
    // polynomial coefficients
    for (int i = 0; i < degree; i++) {

        // the first coefficient has no lower neighbour so its value is the
        // negative of its upper neighbour
        tmp[0] = -rvalue[1];

        // the next coefficient is its lower neighbour minus the product of its
        // upper neighbour and degree
        for (int j = 1; j <= i; j++)
            tmp[j] = rvalue[j - 1] - rvalue[j + 1] * (j + 1);

        // the last coefficient has no upper neighbour so its value is its
        // lower neighbour
        tmp[i + 1] = rvalue[i];

        // copy the coefficients from 'tmp' into 'rvalue'
        for (int j = 0; j <= i + 1; j++) rvalue[j] = tmp[j];
    }

    // a warning message for non-finite coefficients
    for (int i = 0; i < degree + 1; i++) {
        if (!R_FINITE(rvalue[i])) {
            warning("non-finite coefficients produced");
            break;
        }
    }

    UNPROTECT(1);
    return value;
}





int _gcd(int x, int y)
{
    if (x == y) {
        return x;
    }
    else if (x % 2 == 0) {
        if (y % 2 == 0)
            return 2 * _gcd(x/2, y/2);
        else return _gcd(x/2, y);
    }
    else if (y % 2 == 0) {
        return _gcd(x, y/2);
    }
    else if (x > y) {
        return _gcd((x - y)/2, y);
    }
    else return _gcd((y - x)/2, x);
}


int gcd(int x, int y)
{
    if (x == NA_INTEGER || y == NA_INTEGER) return NA_INTEGER;
    x = abs(x), y = abs(y);
    if (x == 0) return y;
    if (y == 0) return x;
    return _gcd(x, y);
}


int lcm(int x, int y)
{
    if (x == NA_INTEGER || y == NA_INTEGER) return NA_INTEGER;
    if (x == 0 || y == 0) return 0;
    x = abs(x), y = abs(y);
    return x / gcd(x, y) * y;
}





/*
SEXP gcd_vectorized(SEXP x, SEXP y)
{
    R_xlen_t len, i, lenx = xlength(x), leny = xlength(y);
    if (lenx == 0 || leny == 0) return allocVector(INTSXP, 0);
    x = PROTECT(coerceVector(x, INTSXP));
    y = PROTECT(coerceVector(y, INTSXP));
    len = fmax(lenx, leny);

    SEXP value = PROTECT(allocVector(INTSXP, len));
    int *ix = INTEGER(x), *iy = INTEGER(y), *ivalue = INTEGER(value);
    for (i = 0; i < len; i++) {
        ivalue[i] = gcd(ix[i % lenx], iy[i % leny]);
    }

    UNPROTECT(3);
    return value;
}
 */





/* old version of 'pgcd'

void pgcd4(int *value, SEXP x, R_xlen_t len, int na_rm)
{
    R_xlen_t i, lenx = xlength(x);
    if (lenx == 0) {
        error("internal error, 'x' should not be length 0");
        return;
    }
    switch(TYPEOF(x)) {
    case INTSXP:
    case LGLSXP:
        break;
    case REALSXP:
    case CPLXSXP:
        x = coerceVector(x, INTSXP);
        break;
    default:
        error("invalid 'type' (%s) of argument", type2char(TYPEOF(x)));
        return;
    }
    int *ix = INTEGER(x);
    if (lenx == 1) {
        int x1 = ix[0];

        if (na_rm && x1 == NA_INTEGER)
            x1 = 0;

        for (i = 0; i < len; i++) {
            value[i] = gcd(value[i], x1);
        }
    }
    else if (lenx == len) {
        if (na_rm) {
            for (i = 0; i < len; i++) {
                value[i] = gcd(value[i], filter_na0(ix[i]));
            }
        }
        else {
            for (i = 0; i < len; i++) {
                value[i] = gcd(value[i], ix[i]);
            }
        }
    }
    else {
        if (na_rm) {
            for (i = 0; i < len; i++) {
                value[i] = gcd(value[i], filter_na0(ix[i % len]));
            }
        }
        else {
            for (i = 0; i < len; i++) {
                value[i] = gcd(value[i], ix[i % len]);
            }
        }
    }
}


SEXP pgcd(SEXP args)
{
    SEXP a, x;
    R_xlen_t n, len = 1;
    int na_rm;

    args = CDR(args);

    na_rm = asLogical(CAR(args));
    args = CDR(args);

    if (args == R_NilValue) return allocVector(INTSXP, 0);

    for (a = args; a != R_NilValue; a = CDR(a)) {
        x = CAR(a);
        switch(TYPEOF(x)) {
        case INTSXP:
        case LGLSXP:
        case NILSXP:
        case REALSXP:
        case CPLXSXP:
            break;
        default:
            error("invalid 'type' (%s) of argument", type2char(TYPEOF((x))));
            return R_NilValue;
        }
        if (len) {
            n = xlength(x);
            len = (n == 0) ? n : ((len < n) ? n : len);
        }
    }

    if (len == 0) return allocVector(REALSXP, 0);

    for (a = args; a != R_NilValue; a = CDR(a)) {
        n = xlength(CAR(a));
        if (len % n) {
            warning("an argument will be fractionally recycled");
            break;
        }
    }

    SEXP value = PROTECT(allocVector(INTSXP, len));
    int *ivalue = INTEGER(value);
    memset(ivalue, 0, len * sizeof(int));

    for (a = args; a != R_NilValue; a = CDR(a)) {
        pgcd4(ivalue, CAR(a), len, na_rm);
    }

    for (a = args; a != R_NilValue; a = CDR(a)) {

        x = CAR(a);
        // the argument must be the same length as 'value' in
        // order to copy 'dim', 'dimnames' and 'names' attributes
        if (xlength(x) == len) {

            // if 'value' does not have a 'dim' attribute
            if (getAttrib(value, R_DimSymbol) == R_NilValue) {

                // if 'value' does not have a 'dim' attribute
                // and 'x' has a 'dim' attribute
                if (getAttrib(x, R_DimSymbol) != R_NilValue) {
                    setAttrib(value, R_DimSymbol, getAttrib(x, R_DimSymbol));
                    setAttrib(value, R_NamesSymbol, R_NilValue);

                    // if 'a' also has a 'dimnames' attribute,
                    // copy them and immediately return 'value'
                    if (getAttrib(x, R_DimNamesSymbol) != R_NilValue) {
                        setAttrib(value, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));
                        UNPROTECT(1);
                        return value;
                    }
                }
                else if (getAttrib(value, R_NamesSymbol) == R_NilValue &&
                    getAttrib(x, R_NamesSymbol) != R_NilValue) {
                    setAttrib(value, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
                }
            }
            else if (conformable(value, x)) {
                if (getAttrib(x, R_DimNamesSymbol) != R_NilValue) {
                    setAttrib(value, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));
                    UNPROTECT(1);
                    return value;
                }
            }
        }
    }

    UNPROTECT(1);
    return value;
}
 */





/* old version of 'ngcd'

void do_gcd3(int *value, SEXP x, int na_rm)
{
    R_xlen_t i, len = xlength(x);

    switch(TYPEOF(x)) {
    case INTSXP:
    case LGLSXP:
        break;
    case REALSXP:
    case CPLXSXP:
        x = coerceVector(x, INTSXP);
        break;
    default:
        error("invalid 'type' (%s) of argument", type2char(TYPEOF(x)));
        return;
    }
    int *ix = INTEGER(x);
    if (na_rm) {
        for (i = 0; i < len; i++) {
            *value = gcd(*value, filter_na0(ix[i]));
        }
    }
    else {
        for (i = 0; i < len; i++) {
            *value = gcd(*value, ix[i]);
        }
    }
}


// n-dimensional gcd
SEXP ngcd(SEXP args)
{
    SEXP a, x;
    int na_rm;

    args = CDR(args);
    na_rm = asLogical(CAR(args)), args = CDR(args);

    if (args == R_NilValue) return ScalarInteger(0);

    for (a = args; a != R_NilValue; a = CDR(a)) {
        x = CAR(a);
        switch(TYPEOF(x)) {
        case INTSXP:
        case LGLSXP:
        case NILSXP:
        case REALSXP:
        case CPLXSXP:
            break;
        default:
            error("invalid 'type' (%s) of argument", type2char(TYPEOF(x)));
            return R_NilValue;
        }
    }

    int value = 0;
    for (a = args; a != R_NilValue; a = CDR(a)) {
        do_gcd3(&value, CAR(a), na_rm);
    }
    return ScalarInteger(value);
}
 */





SEXP pgcd(SEXP args)
{
    SEXP a, x;
    R_xlen_t i, n, len = 1;
    int na_rm;

    args = CDR(args);
    na_rm = asLogical(CAR(args));
    args = CDR(args);

    if (args == R_NilValue) return allocVector(INTSXP, 0);

    for (a = args; a != R_NilValue; a = CDR(a)) {
        x = CAR(a);
        switch(TYPEOF(x)) {
        case INTSXP:
        case LGLSXP:
        case NILSXP:
        case REALSXP:
        case CPLXSXP:
            break;
        default:
            error("invalid 'type' (%s) of argument", type2char(TYPEOF((x))));
            return R_NilValue;
        }
        if (len) {
            n = xlength(x);
            if (n == 0 || n > len)
                len = n;
        }
    }

    if (len == 0) return allocVector(INTSXP, 0);

    for (a = args; a != R_NilValue; a = CDR(a)) {
        if (len % xlength(CAR(a))) {
            warning("an argument will be fractionally recycled");
            break;
        }
    }

    SEXP value = PROTECT(allocVector(INTSXP, len));
    int *ivalue = INTEGER(value);
    memset(ivalue, 0, len * sizeof(int));

    for (a = args; a != R_NilValue; a = CDR(a)) {
        x = CAR(a);
        n = xlength(x);
        switch(TYPEOF(x)) {
        case INTSXP:
        case LGLSXP:
            break;
        case REALSXP:
        case CPLXSXP:
            x = coerceVector(x, INTSXP);
            break;
        default:
            error("invalid 'type' (%s) of argument", type2char(TYPEOF(x)));
            return R_NilValue;
        }
        int *ix = INTEGER(x);
        if (n == 1) {
            int x1 = ix[0];

            if (na_rm && x1 == NA_INTEGER)
                x1 = 0;

            for (i = 0; i < len; i++) {
                ivalue[i] = gcd(ivalue[i], x1);
            }
        }
        else if (n == len) {
            if (na_rm) {
                for (i = 0; i < len; i++) {
                    ivalue[i] = gcd(ivalue[i], filter_na0(ix[i]));
                }
            }
            else {
                for (i = 0; i < len; i++) {
                    ivalue[i] = gcd(ivalue[i], ix[i]);
                }
            }
        }
        else {
            if (na_rm) {
                for (i = 0; i < len; i++) {
                    ivalue[i] = gcd(ivalue[i], filter_na0(ix[i % len]));
                }
            }
            else {
                for (i = 0; i < len; i++) {
                    ivalue[i] = gcd(ivalue[i], ix[i % len]);
                }
            }
        }
    }

    for (a = args; a != R_NilValue; a = CDR(a)) {

        x = CAR(a);
        // the argument must be the same length as 'value' in
        // order to copy 'dim', 'dimnames' and 'names' attributes
        if (xlength(x) == len) {

            // if 'value' does not have a 'dim' attribute
            if (getAttrib(value, R_DimSymbol) == R_NilValue) {

                // if 'value' does not have a 'dim' attribute
                // and 'x' has a 'dim' attribute
                if (getAttrib(x, R_DimSymbol) != R_NilValue) {
                    setAttrib(value, R_DimSymbol, getAttrib(x, R_DimSymbol));
                    setAttrib(value, R_NamesSymbol, R_NilValue);

                    // if 'a' also has a 'dimnames' attribute,
                    // copy them and immediately return 'value'
                    if (getAttrib(x, R_DimNamesSymbol) != R_NilValue) {
                        setAttrib(value, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));
                        UNPROTECT(1);
                        return value;
                    }
                }
                else if (getAttrib(value, R_NamesSymbol) == R_NilValue &&
                    getAttrib(x, R_NamesSymbol) != R_NilValue) {
                    setAttrib(value, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
                }
            }
            else if (conformable(value, x)) {
                if (getAttrib(x, R_DimNamesSymbol) != R_NilValue) {
                    setAttrib(value, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));
                    UNPROTECT(1);
                    return value;
                }
            }
        }
    }

    UNPROTECT(1);
    return value;
}


// n-dimensional gcd
SEXP ngcd(SEXP args)
{
    SEXP a, x;
    R_xlen_t i, n;
    int na_rm;

    args = CDR(args);
    na_rm = asLogical(CAR(args)), args = CDR(args);

    if (args == R_NilValue) return ScalarInteger(0);

    for (a = args; a != R_NilValue; a = CDR(a)) {
        x = CAR(a);
        switch(TYPEOF(x)) {
        case INTSXP:
        case LGLSXP:
        case NILSXP:
        case REALSXP:
        case CPLXSXP:
            break;
        default:
            error("invalid 'type' (%s) of argument", type2char(TYPEOF(x)));
            return R_NilValue;
        }
    }

    int value = 0;
    for (a = args; a != R_NilValue; a = CDR(a)) {
        x = CAR(a);
        n = xlength(x);

        switch(TYPEOF(x)) {
        case INTSXP:
        case LGLSXP:
            break;
        case NILSXP:
            continue;
        case REALSXP:
        case CPLXSXP:
            x = coerceVector(x, INTSXP);
            break;
        default:
            error("invalid 'type' (%s) of argument", type2char(TYPEOF(x)));
            return ScalarInteger(NA_INTEGER);
        }
        int *ix = INTEGER(x);
        if (na_rm) {
            for (i = 0; i < n; i++) {
                value = gcd(value, filter_na0(ix[i]));
            }
        }
        else {
            for (i = 0; i < n; i++) {
                value = gcd(value, ix[i]);
            }
        }
    }
    return ScalarInteger(value);
}





SEXP plcm(SEXP args)
{
    SEXP a, x;
    R_xlen_t i, n, len = 1;
    int na_rm;

    args = CDR(args);
    na_rm = asLogical(CAR(args));
    args = CDR(args);

    if (args == R_NilValue) return allocVector(INTSXP, 0);

    for (a = args; a != R_NilValue; a = CDR(a)) {
        x = CAR(a);
        switch(TYPEOF(x)) {
        case INTSXP:
        case LGLSXP:
        case NILSXP:
        case REALSXP:
        case CPLXSXP:
            break;
        default:
            error("invalid 'type' (%s) of argument", type2char(TYPEOF((x))));
            return R_NilValue;
        }
        if (len) {
            n = xlength(x);
            if (n == 0 || n > len)
                len = n;
        }
    }

    if (len == 0) return allocVector(INTSXP, 0);

    for (a = args; a != R_NilValue; a = CDR(a)) {
        if (len % xlength(CAR(a))) {
            warning("an argument will be fractionally recycled");
            break;
        }
    }

    SEXP value = PROTECT(allocVector(INTSXP, len));
    int *ivalue = INTEGER(value);
    for (i = 0; i < len; i++) ivalue[i] = 1;

    for (a = args; a != R_NilValue; a = CDR(a)) {
        x = CAR(a);
        n = xlength(x);
        switch(TYPEOF(x)) {
        case INTSXP:
        case LGLSXP:
            break;
        case REALSXP:
        case CPLXSXP:
            x = coerceVector(x, INTSXP);
            break;
        default:
            error("invalid 'type' (%s) of argument", type2char(TYPEOF(x)));
            return R_NilValue;
        }
        int *ix = INTEGER(x);
        if (n == 1) {
            int x1 = ix[0];

            if (na_rm && x1 == NA_INTEGER)
                x1 = 1;

            for (i = 0; i < len; i++) {
                ivalue[i] = lcm(ivalue[i], x1);
            }
        }
        else if (n == len) {
            if (na_rm) {
                for (i = 0; i < len; i++) {
                    ivalue[i] = lcm(ivalue[i], filter_na1(ix[i]));
                }
            }
            else {
                for (i = 0; i < len; i++) {
                    ivalue[i] = lcm(ivalue[i], ix[i]);
                }
            }
        }
        else {
            if (na_rm) {
                for (i = 0; i < len; i++) {
                    ivalue[i] = lcm(ivalue[i], filter_na1(ix[i % len]));
                }
            }
            else {
                for (i = 0; i < len; i++) {
                    ivalue[i] = lcm(ivalue[i], ix[i % len]);
                }
            }
        }
    }

    for (a = args; a != R_NilValue; a = CDR(a)) {

        x = CAR(a);
        // the argument must be the same length as 'value' in
        // order to copy 'dim', 'dimnames' and 'names' attributes
        if (xlength(x) == len) {

            // if 'value' does not have a 'dim' attribute
            if (getAttrib(value, R_DimSymbol) == R_NilValue) {

                // if 'value' does not have a 'dim' attribute
                // and 'x' has a 'dim' attribute
                if (getAttrib(x, R_DimSymbol) != R_NilValue) {
                    setAttrib(value, R_DimSymbol, getAttrib(x, R_DimSymbol));
                    setAttrib(value, R_NamesSymbol, R_NilValue);

                    // if 'a' also has a 'dimnames' attribute,
                    // copy them and immediately return 'value'
                    if (getAttrib(x, R_DimNamesSymbol) != R_NilValue) {
                        setAttrib(value, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));
                        UNPROTECT(1);
                        return value;
                    }
                }
                else if (getAttrib(value, R_NamesSymbol) == R_NilValue &&
                    getAttrib(x, R_NamesSymbol) != R_NilValue) {
                    setAttrib(value, R_NamesSymbol, getAttrib(x, R_NamesSymbol));
                }
            }
            else if (conformable(value, x)) {
                if (getAttrib(x, R_DimNamesSymbol) != R_NilValue) {
                    setAttrib(value, R_DimNamesSymbol, getAttrib(x, R_DimNamesSymbol));
                    UNPROTECT(1);
                    return value;
                }
            }
        }
    }

    UNPROTECT(1);
    return value;
}


// n-dimensional lcm
SEXP nlcm(SEXP args)
{
    SEXP a, x;
    R_xlen_t i, n;
    int na_rm;

    args = CDR(args);
    na_rm = asLogical(CAR(args)), args = CDR(args);

    if (args == R_NilValue) return ScalarInteger(1);

    for (a = args; a != R_NilValue; a = CDR(a)) {
        x = CAR(a);
        switch(TYPEOF(x)) {
        case INTSXP:
        case LGLSXP:
        case NILSXP:
        case REALSXP:
        case CPLXSXP:
            break;
        default:
            error("invalid 'type' (%s) of argument", type2char(TYPEOF(x)));
            return R_NilValue;
        }
    }

    int value = 1;
    for (a = args; a != R_NilValue; a = CDR(a)) {
        x = CAR(a);
        n = xlength(x);

        switch(TYPEOF(x)) {
        case INTSXP:
        case LGLSXP:
            break;
        case NILSXP:
            continue;
        case REALSXP:
        case CPLXSXP:
            x = coerceVector(x, INTSXP);
            break;
        default:
            error("invalid 'type' (%s) of argument", type2char(TYPEOF(x)));
            return ScalarInteger(NA_INTEGER);
        }
        int *ix = INTEGER(x);
        if (na_rm) {
            for (i = 0; i < n; i++) {
                value = lcm(value, filter_na1(ix[i]));
            }
        }
        else {
            for (i = 0; i < n; i++) {
                value = lcm(value, ix[i]);
            }
        }
    }
    return ScalarInteger(value);
}





static const R_CallMethodDef callRoutines[] = {
    {"hermite", (DL_FUNC) &hermite       , 1},

    // {"gcd"    , (DL_FUNC) &gcd_vectorized, 2},

    {NULL, NULL, 0}
};


static const R_ExternalMethodDef externalRoutines[] = {
    {"gcd" , (DL_FUNC) &ngcd, -1},
    {"pgcd", (DL_FUNC) &pgcd, -1},
    {"lcm" , (DL_FUNC) &nlcm, -1},
    {"plcm", (DL_FUNC) &plcm, -1},
    {NULL, NULL, 0}
};


void R_init_polynomial(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callRoutines, NULL, externalRoutines);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
