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
			(Z, Flow) = EMDU.EMDUnifrac_weighted_flow(T, l, nodes_in_order, P, Q)
			t1 = timeit.default_timer()
			EMDUnifrac_flow_times.append(t1-t0)
			#EMDUnifrac no flow
			t0 = timeit.default_timer()
			Z = EMDU.EMDUnifrac_weighted(T, l, nodes_in_order, P, Q)
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


##########
# Here is the real data that is analyzed in the paper:
# tree_str = "(Acidobacteria:0.03031,Actinobacteria:0.01878,Bacteroidetes:0.00530,Chlorobi:0.01402,Fusobacteria:0.10722,Lentisphaerae:0.04241,Proteobacteria:0.01667,Spirochaetes:0.03298,Synergistetes:0.03566,Tenericutes:0.00230,Verrucomicrobia:0.03222)Bacteria;"
# phyla=["Acidobacteria","Actinobacteria","Bacteroidetes","Chlorobi","Firmicutes","Fusobacteria","Lentisphaerae","Proteobacteria","SR1","Spirochaetes","Synergistetes","TM7","Tenericutes","Verrucomicrobia"]
# healthy = np.loadtxt("/Users/dkoslicki/Dropbox/EMDUnifrac/RealData/HealthyTotals.txt")
# UC = np.loadtxt("/Users/dkoslicki/Dropbox/EMDUnifrac/RealData/UCTotals.txt")
# healthy = np.array([5.29000000e-04, 1.19146400e+00, 2.46148900e+01, 5.29000000e-04, 1.85075170e+01, 1.49653000e-01, 3.41100000e-03, 4.30437200e+00, 0.00000000e+00, 1.77374000e-01, 8.38000000e-03, 0.00000000e+00, 4.18820000e-02, 0.00000000e+00])  # obtained via np.loadtxt("/Users/dkoslicki/Dropbox/EMDUnifrac/RealData/HealthyTotals.txt")
# UC = np.array([ 0.00000000e+00, 5.64607000e-01, 6.00750200e+00, 0.00000000e+00, 7.68366300e+00, 5.19120000e-02, 8.80000000e-04, 1.64073200e+00, 1.57500000e-03, 1.57500000e-03, 0.00000000e+00, 0.00000000e+00, 4.75550000e-02, 0.00000000e+00])   # obtained via np.loadtxt("/Users/dkoslicki/Dropbox/EMDUnifrac/RealData/UCTotals.txt")
# envs = dict()
# for i in range(len(phyla)):
# 	envs[phyla[i]] = dict()
# 	envs[phyla[i]]["healthy"] = healthy[i]
# 	envs[phyla[i]]["UC"] = UC[i]
# (Tint,lint,nodes_in_order) = EMDU.parse_tree(tree_str)
# (nodes_weighted, samples) = EMDU.parse_envs(envs,nodes_in_order)
# (Z, diffab) = EMDU.EMDUnifrac_weighted(Tint,lint,nodes_in_order,nodes_weighted[samples[0]],nodes_weighted[samples[1]])



