import numpy as np
import dendropy
import matplotlib.pyplot as plt
import warnings
import EMDUnifrac as EMD

# Very simple example first off

# Tree of depth 2
b = 2

# all branch lengths 1
l = dict()
for i in range(2, 2**(b+1)):
	l[i] = 1

# number of vertices
d = 7

# populate the w guy (iterations line 1)
w = np.zeros((d, d))  # note, these are vectors of length d, not d-1 as in the manuscript
for i in np.array(range(2**b, 2**(b+1)))-1:
	w[i, i] = l[i]

# iteration line 4
for i in np.array(range(2**b-1, 2, -1))-1:  # negative stride
	w[i] = w[2*i] + w[2*i+1]  # add up the entire vectors ("all of the subtrees for the daughters of i)
	w[i, i] = l[i]
