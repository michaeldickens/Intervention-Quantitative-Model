"""

sum.py
------

Author: Michael Dickens <mdickens93@gmail.com>
Created: 2016-04-01

""" 

import numpy as np
import scipy.stats as stats

# num_buckets = 1000
# exp_offset = 50
# step = 1.25

num_buckets = 200
exp_offset = 30
step = 1.35

def dist_to_buckets(cdf):
    """
    Bucket `i` gives the probability density between `step^i` and `step^(i+1)` (subtracting `bucket_offset`).
    """
    buckets = []
    for i in range(num_buckets):
        e = i - exp_offset
        buckets.append(cdf(step**(e+1)) - cdf(step**e))
    return buckets
 
# TODO: If a[i] is much larger than b[j] or vice versa, the sum is
# approximately equal to the max of the two. You can exploit this to
# improve runtime but it's a bit tricky.
def sum_buckets(a, b):
    res = [ 0 for _ in range(num_buckets) ]
    for i in range(num_buckets):
        for j in range(num_buckets):
            p = a[i] * b[j]
            index = int(np.floor(np.log(step**(i - exp_offset) + step**(j - exp_offset)) / np.log(step))) + exp_offset
            if index >= num_buckets:
                # this makes a better approximation than just ignoring
                # values that go above step**num_buckets
                index = num_buckets - 1
            res[index] += p
    return res

a = dist_to_buckets(lambda x: stats.lognorm.cdf(x, 1))
b = dist_to_buckets(lambda x: stats.lognorm.cdf(x, 2))
s = sum_buckets(a, b)
print sum(s)
