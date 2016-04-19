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
#tree_sizes = [10,100,1000]
tree_sizes = range(10,100000,5000)
mean_EMDUnifrac_times = np.zeros(len(tree_sizes), dtype = np.float64)
mean_EMDUnifrac_flow_times = np.zeros(len(tree_sizes), dtype = np.float64)
mean_FU_times = np.zeros(len(tree_sizes), dtype = np.float64)
#Iterate over tree sizes
for tree_size in tree_sizes:
	FU_times = np.zeros(num_samples*num_trees, dtype = np.float64)
	EMDUnifrac_times = np.zeros(num_samples*num_trees, dtype = np.float64)
	EMDUnifrac_flow_times = np.zeros(num_samples*num_trees, dtype = np.float64)
	for tree_it in range(num_trees):
		t = Tree()
		t.populate(tree_size, random_branches = True)
		tree_str = t.write(format=1)
		tr = DndParser(tree_str, UniFracTreeNode)
		(T,l,nodes_in_order) = EMDU.parse_tree(tree_str)
		i=0
		for it in range(num_samples):
			envs = EMDU.simulate_data(t.get_leaf_names())  # FastUnifrac can only take weight on leaf nodes
			(envs_prob_dict, samples) = EMDU.parse_envs(envs, nodes_in_order)
			P = envs_prob_dict[samples[0]]
			Q = envs_prob_dict[samples[1]]
			#EMDUnifrac with flow
			t0 = timeit.default_timer()
			(Z, Flow) = EMDU.EMDUnifrac_weighted_flow(T, l, nodes_in_order, P, Q)
			t1 = timeit.default_timer()
			EMDUnifrac_flow_times[i] = t1-t0
			#EMDUnifrac no flow
			t0 = timeit.default_timer()
			Z = EMDU.EMDUnifrac_weighted(T, l, nodes_in_order, P, Q)
			t1 = timeit.default_timer()
			EMDUnifrac_times[i] = t1-t0
			#FastUnifrac weighted
			t0 = timeit.default_timer()
			res = fast_unifrac(tr, envs, weighted=True, modes=set(['distance_matrix']))
			t1 = timeit.default_timer()
			FU_times[i] = t1-t0
			i = i+1
	#Print results
	mean_EMDUnifrac_times[tree_sizes.index(tree_size)] = EMDUnifrac_times.mean()
	mean_EMDUnifrac_flow_times[tree_sizes.index(tree_size)] = EMDUnifrac_flow_times.mean()
	mean_FU_times[tree_sizes.index(tree_size)] = FU_times.mean()
	#print("For tree size %s" % str(tree_size))
	#print("Mean EMDUnifrac flow time: %f" % EMDUnifrac_flow_times.mean())
	#print("Mean EMDUnifrac no flow time: %f" % EMDUnifrac_times.mean())
	#print("Mean FastUnifrac time: %f" % FU_times.mean())
	#print("Speed improvement factor: %f" % (FU_times.mean()/EMDUnifrac_times.mean()))

np.savetxt('EMDU_mean_times.txt', mean_EMDUnifrac_times, delimiter=',')
np.savetxt('EMDU_flow_mean_times.txt', mean_EMDUnifrac_flow_times, delimiter=',')
np.savetxt('FastUnifrac_mean_times.txt', mean_FU_times, delimiter=',')





