import numpy as np
import dendropy
import matplotlib.pyplot as plt
import warnings
import EMDUnifrac as EMD
import random
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
S_unflat = [S1, S2, S3]
t = len(S_unflat)
U = []
k = len(S)

# Sample size variance estimate
K = 3

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

# line 8: generate a random sample without repetition, forcing distribution to respect population subsets

for i in range(t):
	num_sample = int(np.floor(len(S_unflat[i]) * K / float(len(S))))
	U.append(random.sample(S_unflat[i], k=max(1, num_sample)))  # TODO: check with jason that taking max with 1 is ok

# line 13, covariance matrix
BigSigmaU = np.cov(np.array([item for sublist in U for item in sublist]).transpose())  # TODO: check with Jason this works (note the transpose) since it differs from his calculations with U = {p2, p3, p5}

littlesigma = np.zeros(d)

# line 15
for i in range(d):
	littlesigma[i] = np.dot(w[i].transpose(), np.matmul(BigSigmaU, w[i]))

O = np.argsort(littlesigma)

# Algorithm 3
tau = np.median(littlesigma)
nt = np.where(littlesigma[O] <= tau)[0][-1]  # largest index <= threshold
alpha = .5
D = dict()

# Set difference between {2,...,(d-1)/2} - {O(1),...,O(nt)}
J = [item for item in range(1, int((d-1)/2.)) if item not in [O[i] for i in range(O[nt])]]
J.reverse()

# Iterations
for i in range(k):
	si = S[i]
	for j in J:
		si[j] = si[j] + l[2*(j+1)-1]*si[2*(j+1)-1] + l[2*(j+1)+1-1]*si[2*(j+1)+1-1]

# line 8
F = set()
for i in range(k-1):
	for j in range(i+1, k-1):
		F.add((i, j))

for i in range(d-1, nt, -1):
	if not F:
		break
	for tup in list(F):
		a = tup[0]
		b = tup[1]
		if tup in D:
			D[tup] = D[tup] + abs(S[a][O[i]] - S[b][O[i]])
		else:
			D[tup] = abs(S[a][O[i]] - S[b][O[i]])
		if D[tup] > alpha:
			F.remove(tup)

# TODO: looks like the distances that I'm getting are exactly 1/2 of what they should be...

# Larger example
b = 10
S1 = np.random.dirichlet(np.random.beta(1, 2, size=2**(b+1)-1), size=1000)
S2 = np.random.dirichlet(np.random.beta(2, 1, size=2**(b+1)-1), size=1000)
S3 = np.random.dirichlet(np.random.beta(5, .01, size=2**(b+1)-1), size=1000)
