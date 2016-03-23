"""

universe.py
-----------

Author: Michael Dickens <mdickens93@gmail.com>
Created: 2016-03-09

Simulates the universe to see how you'd expect the distribution of charities to work if they're completely random.

""" 

import numpy as np
import scipy.sparse as sparse

def utility(mat, vec, num_iters):
    for _ in range(num_iters):
        vec = np.dot(mat, vec) / len(vec)
    return sum(vec)

def rand_matrix(size):
    data = np.zeros((size, size))
    for _ in range(int(size*np.sqrt(size))):
        i = np.random.randint(size)
        j = np.random.randint(size)
        data[i][j] = np.random.rand()
    mat = data
    return mat

def score(mat, vec, i):
    num_iters = 1
    vec1 = np.array([ x for x in vec ])
    vec1[i] += 1
    return abs(utility(mat, vec, num_iters) - utility(mat, vec1, num_iters))

size = 1000 # 1000 is quick, 5000 takes a while
mat = rand_matrix(size)
vec = np.random.rand(size)
scores = [ score(mat, vec, i) for i in range(size) ]
print "Value"
print "\n".join([ "%f" % (x) for (i, x) in enumerate(scores) ])
