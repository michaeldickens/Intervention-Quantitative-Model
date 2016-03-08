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

def pareto(x):
    global m, a
    return stats.pareto.pdf(x, a, scale=m)

def lognorm(mu, x):
    global sigma
    return stats.lognorm.pdf(x, sigma, scale=np.exp(mu))

def prior_EV_calculation(mu):
    return integrate.quad(lambda x: pareto(lognorm(mu, x)),
                   0, np.inf)

def posterior(val):
    mu = mu_posterior(val)
    prior = pareto(val)
    conditional = lognorm(mu, val)
    denom = prior_EV_calculation(mu)
    return prior * conditional / denom

def mu_posterior(mu_p, val):
    global m, a, sigma
    prior = pareto(mu_p)
    conditional = lognorm(mu_p, val)
    denom = pareto(val)
    return prior * conditional / denom

print posterior(1)
