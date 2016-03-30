"""

MH.py
-----

Author: Michael Dickens <mdickens93@gmail.com>
Created: 2016-03-29

Implementation of the Metropolis-Hastings algorithm for approximating Bayesian posteriors. The main purpose of implementing here is so I can easily test and get a handle on the algorithm before implementing it in Visual Basic.

""" 

import numpy as np

def get_candidate(q, given=None):
    pass
    
def metropolis_hastings(q, num_iters):
    """
    `q`: proposal distribution
    """
    x = get_candidate(q)
    for _ in range(num_iters):
        x_cand = get_candidate(q, x)
        p_accept = min(1, (q(x, x_cand) * pi(x_cand)) / \
                       (q(x_cand, x) * pi(x)))
        if np.random.rand() < p_accept:
            x = x_cand
