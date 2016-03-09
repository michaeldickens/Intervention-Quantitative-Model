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

x0 = 100.0
a = 1.0/3
sigma = 1.5

# TODO: some (all?) of these PDFs should be CDFs. We can't get posterior probability that x = k but we can get probability that x < k.

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

def prior_measurement(m):
    return integrate.quad(lambda u: lognorm_pdf(m / u) * pareto_pdf(u),
                          0, np.inf)

print prior_measurement(100)
print prior_measurement(1000)
print prior_measurement(10000)
