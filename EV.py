"""

EV.py
----------------

Author: Michael Dickens <mdickens93@gmail.com>
Created: 2016-03-07

Probability distribution of intervention cost-effectiveness.

""" 

import numpy as np
import scipy.integrate as integrate
import scipy.stats as stats

m = 1.0/3
a = 100.0
sigma = 1.5

# TODO: some (all?) of these PDFs should be CDFs. We can't get posterior probability that x = k but we can get probability that x < k.

def pareto_pdf(x):
    global m, a
    return stats.pareto.pdf(x, a, scale=m)

def lognorm_pdf(mu, x):
    global sigma
    return stats.lognorm.pdf(x, sigma, scale=mu)

def pareto_cdf(x):
    global m, a
    return stats.pareto.cdf(x, a, scale=m)

def lognorm_cdf(x):
    global sigma
    return stats.lognorm.cdf(x, sigma)

"""
Probability of getting `measurement` if true utility is `utility`.
"""
def p_measurement(measurement, utility):
    return lognorm_pdf(1, measurement / float(utility))
    # return lognorm_pdf(utility, measurement)

def get_const(measurement):
    return integrate.quad(lambda u: posterior_density(u, measurement),
                          0, np.inf)[0]

def posterior_density(utility, measurement, const=1):
    prior = pareto_pdf(utility)
    update = integrate.quad(lambda x: p_measurement(measurement, x), 0, utility)[0]
    return prior * update / const

"""
Probability that `Utility < utility` given `measurement`.
"""
def posterior(utility, measurement, const=1):
    prior = pareto_cdf(utility)
    update = integrate.quad(lambda x: p_measurement(measurement, x), 0, utility)[0]
    return prior * update / const

const = get_const(1000)
const = 1
print const
print posterior(1000, 1000, const)
print posterior(10000, 1000, const)
print posterior(100000, 1000, const)
