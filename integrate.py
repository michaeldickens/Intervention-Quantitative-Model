"""

integrate.py
------------

Author: Michael Dickens <mdickens93@gmail.com>
Created: 2016-03-17

Accurately compute numerical integrals over big ranges where function drops off quickly.

""" 

import numpy as np

def trapezoid(f, x_lo, x_hi):
    y_lo = f(x_lo)
    y_hi = f(x_hi)
    avg = (y_lo + y_hi) / 2
    delta = x_hi - x_lo
    total = avg * delta
    error = max(abs(y_hi - avg) * delta, abs(avg - y_lo) * delta)
    return (total, error)

"""
Numeerically computes the integral of a function with a very large
upper bound `x_f` where the function drops off rapidly.

`error` gives an upper bound on how far off the result could be from
the true value.
"""
def integral(f, x_0, x_f, limit=2**10):
    total = 0
    error = 0

    if x_f == np.inf:
        x_f = 1e99
    if x_0 == 0:
        t, e = trapezoid(f, 0, 1)
        total += t
        error += e
        x_0 = 1

    step = (x_f / x_0)**(1.0/(limit + 1))
    x_lo = x_0
    x_hi = x_lo * step

    for _ in range(limit):
        t, e = trapezoid(f, x_lo, x_hi)
        total += t
        error += e
        x_lo = x_hi
        x_hi *= step

    return (total, error)
    


