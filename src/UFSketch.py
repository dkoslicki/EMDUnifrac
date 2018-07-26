import numpy as np
import dendropy
import matplotlib.pyplot as plt
import warnings
import EMDUnifrac as EMD

# Very simple example first off

# Tree of depth 2
b = 2

# The sample distributions
p1 = np.array([0, 0, 0, 0, 0, 1/2., 1/2.])
p2 = np.array([0, 0, 1/3., 0, 0, 1/3., 1/3.])
p3 = np.array([0, 0, 0, 1/2., 1/2., 0, 0])
p4 = np.array([0, 2/3., 0, 1/6., 1/6., 0, 0])
p5 = np.array([0, 1/2., 1/2., 0, 0, 0, 0])

# The partitions
S1 = [p1, p2]
S2 = [p3, p4]
S3 = [p5]

# all samples
S = S1 + S2 + S3

# all branch lengths 1
l = dict()
for i in np.array(range(2, 2**(b+1)))-1:
	l[i] = 1

# number of vertices
d = 2**(b+1)-1

# populate the w guy (iterations line 1)
w = np.zeros((d, d))  # note, these are vectors of length d, not d-1 as in the manuscript
for i in np.array(range(int((d-1)/2.), d+1))-1:
	w[i, i] = l[i]

# iteration line 4
for i in np.array(range(int((d-1)/2.), 1, -1))-1:  # negative stride
	w[i] = w[2*(i+1)-1] + w[2*(i+1)+1-1]  # add up the entire vectors ("all of the subtrees for the daughters of i)
	w[i, i] = l[i]
