"""

EV.py
----------------

Author: Michael Dickens <mdickens93@gmail.com>
Created: 2016-03-07

Probability distribution of intervention cost-effectiveness.

""" 

import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import scipy.stats as stats

robust_sigma = 1.0 # 95% chance true value is within factor of  2e (~5)
mid_sigma = 2.0    # 95% chance true value is within factor of  6e (~16)
weak_sigma = 6.0   # 95% chance true value is within factor of 12e (~32)

x0 = 100.0
a = 1.0/3
sigma = mid_sigma

def pareto_pdf(x):
    global x0, a
    return stats.pareto.pdf(x, a, scale=x0)

def lognorm_pdf(x):
    global sigma
    return stats.lognorm.pdf(x, sigma)

def pareto_cdf(x):
    global x0, a
    return stats.pareto.cdf(x, a, scale=x0)

def lognorm_cdf(x):
    global sigma
    return stats.lognorm.cdf(x, sigma)

"""
Prior probability that `Measurement = m`, possibly with an upper
bound for utility `u`.
"""
def measurement_prior(m, u_max=np.inf):
    # TODO: Python has a lot of trouble with this integral, it
    # converges too slowly
    return integrate.quad(lambda u: lognorm_pdf(m / u) * pareto_pdf(u),
                          0, u_max, epsabs=0, limit=100)[0]

def log_measurement_prior(m, u_max=np.inf):
    return integrate.quad(lambda u: lognorm_pdf(m - u) * pareto_pdf(np.exp(u)))

"""
Returns a probability function `P(U|M = m) that, when called,
gives P(U = u|M = m)`.
"""
def posterior(m):
    const = measurement_prior(m)
    def f(u):
        prior = pareto_pdf(u)
        update = lognorm_pdf(m / u)
        return prior * update / const
    return f

"""
Returns a probability function that, when called, gives
`P(U < u|M = m)`.
"""
def posterior_cumulative_1(m):
    pdf = posterior(m)
    def f(u):
        return integrate.quad(pdf, 0, u)[0]
    return f

"""
Returns a probability function that, when called, gives
`P(U < u|M = m)`.
"""
def posterior_cumulative_2(m):
    const = measurement_prior(m, u_max=1e10)
    def f(u):
        prior = pareto_cdf(u)
        update = measurement_prior(m, u_max=u)
        return prior * update / const
    return f

m = 10000
cdf = posterior_cumulative_2(m)
def test(u):
    global pdf, cdf
    print "%.3f %.3f" % (pareto_cdf(u), cdf(u))

test(1000)
test(10000)
test(100000)
test(1000000)
