import EMDUnifrac as EMDU
from cogent.parse.tree import DndParser  # For comparison to FastUnifrac
from cogent.maths.unifrac.fast_unifrac import fast_unifrac  # For comparison to FastUnifrac
from cogent.maths.unifrac.fast_tree import UniFracTreeNode  # For comparison to FastUnifrac
import numpy as np
import timeit
from ete2 import Tree  # To generate random trees

#Timing for EMDUnifrac and FastUnifrac. Figure XXXX

# Now for the full timing
num_samples = 5
num_trees = 5
tree_sizes = range(10,100000,5000)  # Number of leaf nodes in the tree
mean_EMDUnifrac_times = np.zeros(len(tree_sizes), dtype = np.float64)
mean_EMDUnifrac_flow_times = np.zeros(len(tree_sizes), dtype = np.float64)
mean_FastUnifrac_times = np.zeros(len(tree_sizes), dtype = np.float64)
#Iterate over tree sizes
for tree_size in tree_sizes:
	FastUnifrac_times = list()
	EMDUnifrac_times = list()
	EMDUnifrac_flow_times = list()
	for tree_it in range(num_trees):
		t = Tree()
		t.populate(tree_size, random_branches = True)  # Randomly generate tree with random branch lengths
		tree_str = t.write(format=1)  # Convert to newick and preserve branch lengths
		tr = DndParser(tree_str, UniFracTreeNode)  # FastUnifrac data parsing
		(T,l,nodes_in_order) = EMDU.parse_tree(tree_str)  # EMDUnifrac data parsing
		i=0
		for it in range(num_samples):
			envs = EMDU.simulate_data(t.get_leaf_names())  # FastUnifrac can only take weight on leaf nodes
			(envs_prob_dict, samples) = EMDU.parse_envs(envs, nodes_in_order)
			P = envs_prob_dict[samples[0]]
			Q = envs_prob_dict[samples[1]]
			#EMDUnifrac with flow
			t0 = timeit.default_timer()
			(Z, Flow, diffab) = EMDU.EMDUnifrac_weighted_flow(T, l, nodes_in_order, P, Q)
			t1 = timeit.default_timer()
			EMDUnifrac_flow_times.append(t1-t0)
			#EMDUnifrac no flow
			t0 = timeit.default_timer()
			(Z, diffab) = EMDU.EMDUnifrac_weighted(T, l, nodes_in_order, P, Q)
			t1 = timeit.default_timer()
			EMDUnifrac_times.append(t1-t0)
			#FastUnifrac_ weighted
			t0 = timeit.default_timer()
			res = fast_unifrac(tr, envs, weighted=True, modes=set(['distance_matrix']))
			t1 = timeit.default_timer()
			FastUnifrac_times.append(t1-t0)
			i = i+1
	#Save means
	mean_EMDUnifrac_times[tree_sizes.index(tree_size)] = np.array(EMDUnifrac_times).mean()
	mean_EMDUnifrac_flow_times[tree_sizes.index(tree_size)] = np.array(EMDUnifrac_flow_times).mean()
	mean_FastUnifrac_times[tree_sizes.index(tree_size)] = np.array(FastUnifrac_times).mean()

#  Export all mean times
np.savetxt('EMDU_mean_times.txt', mean_EMDUnifrac_times, delimiter=',')
np.savetxt('EMDU_flow_mean_times.txt', mean_EMDUnifrac_flow_times, delimiter=',')
np.savetxt('FastUnifrac__mean_times.txt', mean_FastUnifrac_times, delimiter=',')


