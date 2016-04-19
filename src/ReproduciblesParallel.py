# This is the same as Reproducibles.py, but since the timing is embarrassingly parallel, I'll parallelize over the tree size
# simply to decrease amount of time to wait for results.
import EMDUnifracPackage as EMDU
from cogent.parse.tree import DndParser
from cogent.maths.unifrac.fast_unifrac import fast_unifrac
from cogent.maths.unifrac.fast_tree import UniFracTreeNode
import numpy as np
import timeit
from ete2 import Tree
import timeit
from itertools import *
from multiprocessing import Pool, freeze_support

num_threads = 47
tree_sizes = range(10,110010,5000)  # Number of tree leaves to iterate over
num_trees = 10
num_samples = 10


def timing_star(arg):
	return timing(*arg)


def timing(tree_size, num_trees, num_samples):
	FastUnifrac_times = list()
	EMDUnifrac_times = list()
	EMDUnifrac_flow_times = list()
	for tree_it in range(num_trees):
		t = Tree()
		t.populate(tree_size, random_branches = True)
		tree_str = t.write(format=1)
		tr = DndParser(tree_str, UniFracTreeNode)
		(T,l,nodes_in_order) = EMDU.parse_tree(tree_str)
		for it in range(num_samples):
			envs = EMDU.simulate_data(t.get_leaf_names())  # FastUnifrac can only take weight on leaf nodes
			(envs_prob_dict, samples) = EMDU.parse_envs(envs, nodes_in_order)
			P = envs_prob_dict[samples[0]]
			Q = envs_prob_dict[samples[1]]
			#EMDUnifrac with flow
			t0 = timeit.default_timer()
			(Z, Flow) = EMDU.EMDUnifrac_weighted_flow(T, l, nodes_in_order, P, Q)
			t1 = timeit.default_timer()
			EMDUnifrac_flow_times.append(t1-t0)
			#EMDUnifrac no flow
			t0 = timeit.default_timer()
			Z = EMDU.EMDUnifrac_weighted(T, l, nodes_in_order, P, Q)
			t1 = timeit.default_timer()
			EMDUnifrac_times.append(t1-t0)
			#FastUnifrac weighted
			t0 = timeit.default_timer()
			res = fast_unifrac(tr, envs, weighted=True, modes=set(['distance_matrix']))
			t1 = timeit.default_timer()
			FastUnifrac_times.append(t1-t0)
	return  (np.array(EMDUnifrac_times).mean(), np.array(EMDUnifrac_flow_times).mean(), np.array(FastUnifrac_times).mean())

#  Parallelize over the tree sizes
pool = Pool(processes = num_threads)
results = pool.map(timing_star, izip(tree_sizes, repeat(num_trees), repeat(num_samples)))
pool.close()
pool.join()
FastUnifrac_times = list()
EMDUnifrac_times = list()
EMDUnifrac_flow_times = list()
# Gather results
for result in results:
	EMDUnifrac_times.append(result[0])
	EMDUnifrac_flow_times.append(result[1])
	FastUnifrac_times.append(result[2])


EMDUnifrac_times = np.array(EMDUnifrac_times)
EMDUnifrac_flow_times = np.array(EMDUnifrac_flow_times)
FastUnifrac_times = np.array(FastUnifrac_times)

# Export results
np.savetxt('EMDU_mean_times.txt', EMDUnifrac_times, delimiter=',')
np.savetxt('EMDU_flow_mean_times.txt', EMDUnifrac_flow_times, delimiter=',')
np.savetxt('FastUnifrac_mean_times.txt', FastUnifrac_times, delimiter=',')



