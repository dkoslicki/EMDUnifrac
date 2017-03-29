import numpy as np
import sys
import EMDUnifrac as EMDU

# This is to be used as a package. I.e. import CAMI, and then use the included functions as needed. See ParseCAMIParallel.py

# It should be noted that using the profiles misses out on taxa not represented in either of the profiles
# So the trees are constructed *only* from the taxa present in one of the profiles


def CAMI_read_taxonomy_file(file_path):
	'''
	Parse a CAMI profile file
	:param file_path:  location of the profile file
	:return: tax_ids, tax_path, weights
	'''
	tax_path = list()
	tax_ids = list()
	weights = dict()
	with open(file_path, 'r') as read_handler:
		for line in read_handler:
			line = line.rstrip()
			if len(line) == 0:
				continue  # skip blank lines
			if line[0] == '@' and line[1] == '@':
				headers = line.strip().split()
				for header_iter in range(len(headers)):
					header = headers[header_iter]
					header = header.replace('@','')
					if header == 'TAXID':
						tax_id_pos = header_iter
					elif header == 'RANK':
						rank_pos = header_iter
					elif header == 'TAXPATH':
						tax_path_pos = header_iter
					elif header == 'TAXPATHSN':
						tax_path_sn_pos = header_iter
					elif header == 'PERCENTAGE':
						percentage_pos = header_iter
			if line[0] in ['@', '#']:
				continue  # skip comment or header
			if not all([isinstance(x, int) for x in [tax_id_pos, tax_path_pos, percentage_pos]]):
				print("Appears the headers TAXID, TAXPATH, and PERCENTAGE are missing from the header (should start with line @@)")
				sys.exit(2)
			temp_split = line.split('\t')
			tax_path.append(temp_split[tax_path_pos].strip())  # add the whole taxpath
			tax_ids.append(temp_split[tax_id_pos].strip())  # just terminal tax ID
			weights[temp_split[tax_id_pos].strip()] = float(temp_split[percentage_pos].strip())  # the associated weight
			# NEED TO CHECK THAT EACH ONE OF THE TAXIDS IS UNIQUE!!
	return tax_ids, tax_path, weights


def CAMI_form_graph(tax_path_1, tax_path_2):
	'''
	For the graph for the unifrac computation using the tax_paths from a CAMI profile file
	:param tax_path_1: get from CAMI_read_taxonomy_file
	:param tax_path_2: get from CAMI_read_taxonomy_file
	:return: Tint, lint, nodes_in_order, nodes_to_index
	'''
	all_taxpaths = list(set(tax_path_1) | set(tax_path_2))
	# form the graph
	nodes_dict = dict()  # keys are nodes, values are the ancestors tupled with branch lengths
	for tax_path in all_taxpaths:
		tax_ids = tax_path.split('|')
		for i in range(len(tax_ids) - 1, -1, -1):  # iterate backwards through list
			if tax_ids[i] == '':
				continue
			parent_id = '-1'  # in case of unknown path
			branch_len = 1  # This is where we could add the actual branch lengths
			for j in range(i - 1, -1, -1):  # traverse up the taxpath until a non '||' is found, watch out for repeat taxids
				if tax_ids[j] == '' or tax_ids[j] == tax_ids[j+1]:
					branch_len += 1  # you found a '||' so add one to the taxonomy
					continue
				# found the parent, so add the clade
				parent_id = tax_ids[j]
				if tax_ids[i] in nodes_dict:
					if nodes_dict[tax_ids[i]][0] != parent_id:
						print("Warning, taxonomies inconsistent, arbitrarily choosing one of them, but you should be aware.")
				else:
					nodes_dict[tax_ids[i]] = (parent_id, branch_len)
				break
			nodes_dict[tax_ids[i]] = (parent_id, branch_len) ################# If all else fails, connect it to the root
	nodes_in_order = set(nodes_dict.keys())
	nodes_in_order = list(nodes_in_order.union(set([val[0] for val in nodes_dict.values()])))
	# Now create the Tint and lint
	Tint = dict()
	lint = dict()
	nodes_to_index = dict(zip(nodes_in_order, range(len(nodes_in_order))))
	for i in range(len(nodes_in_order)):
		node = nodes_in_order[i]
		if node in nodes_dict:  # It's in there, so pull off the ancestor
			parent = nodes_dict[node][0]
			Tint[i] = nodes_to_index[parent]
			lint[nodes_to_index[node], nodes_to_index[parent]] = nodes_dict[node][1]
		else:  # Node isn't in there, so it's a terminal parent node. Attach it to -1, which will be the last entry of nodes_in_order
			Tint[i] = len(nodes_in_order)
			lint[nodes_to_index[node], len(nodes_in_order)] = 1
	# Append the root node
	nodes_to_index['-1'] = len(nodes_in_order)
	nodes_in_order.append('-1')
	return Tint, lint, nodes_in_order, nodes_to_index


# Make the percentages dict
def CAMI_make_percentages(weights, nodes_to_index):
	'''
	Helper for CAMI_get_probability_distribution
	:param weights: From CAMI_read_taxonomy
	:param nodes_to_index: from CAMI_read_taxonomy
	:return: percentages_dict for use in CAMI_get_probability_distribution
	'''
	percentages_dict = dict()
	weights_keys = weights.keys()
	for i in range(len(weights_keys)):
		node = weights_keys[i]
		percentages_dict[nodes_to_index[node]] = weights[node]
	return percentages_dict

# Interestingly, the percentages_dict is basically the "pushed up" vector. To compute Weighted unifrac, just need to do
# the L1 norm on this vector... but let's get it in the standard non-pushed up version so we can work with flows etc.


def CAMI_get_probability_distribution(nodes_in_order, Tint, percentages_dict):
	'''
	This function gets the probability distributions for the CAMI profile that was used to make percentages_dict.
	:param nodes_in_order: yup
	:param Tint: same ol same ol
	:param percentages_dict: see this function above
	:return: probability_distribution for use in the EMDUnifrac functions
	'''
	# Form the probability distribution for input.
	# The CAMI format automatically sums up over lower taxa, so I need to subtract descendant weights
	npTint_keys = np.array(Tint.keys())
	npTint_values = np.array(Tint.values())
	probability_distribution = np.zeros(len(nodes_in_order), dtype=np.float64)
	for i in range(len(nodes_in_order)):
		if i in percentages_dict:
			probability_distribution[i] = percentages_dict[i]
			# now subtract all of the descendant weights. This is so we end up with a probability distribution.
			desc_nodes = npTint_keys[np.where(npTint_values == i)[0]]  # Get all the descendent nodes to i, one step descendants
			# one step desc is ok to use since we directly connect taxa with intermediate missing taxids: 1|||2 -> 1->2
			for desc_node in desc_nodes:  # just the one-step descendants
				if desc_node in percentages_dict:
					# threshold to zero to get around things like -1.0e-16
					probability_distribution[i] = max([probability_distribution[i] - percentages_dict[desc_node], 0])
	# normalize to make sure it's a prob. distribution
	return probability_distribution / np.sum(probability_distribution)


def test_real_data():
	in_file_1 = '../data/test1.profile'
	in_file_2 = '../data/test2.profile'
	tax_ids_1, tax_path_1, weights_1 = CAMI_read_taxonomy_file(in_file_1)
	tax_ids_2, tax_path_2, weights_2 = CAMI_read_taxonomy_file(in_file_2)
	Tint, lint, nodes_in_order, nodes_to_index = CAMI_form_graph(tax_path_1, tax_path_2)
	perc1 = CAMI_make_percentages(weights_1, nodes_to_index)
	prob1 = CAMI_get_probability_distribution(nodes_in_order, Tint, perc1)
	perc2 = CAMI_make_percentages(weights_2, nodes_to_index)
	prob2 = CAMI_get_probability_distribution(nodes_in_order, Tint, perc2)
	(Z, F, diffab) = EMDU.EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, prob1, prob2)
	assert abs(Z - 0.16659746613180418) < .00001
	assert F == {(68, 68): 0.0035418936986314664, (13, 13): 0.00073486688381858389, (70, 70): 0.048404262374067412, (52, 52): 0.0050164245725866563, (71, 71): 0.014268652871098495, (56, 56): 0.010557928557899702, (78, 64): 0.0090384198880965536, (62, 64): 0.00049539987979265593, (51, 51): 0.0026161596968189331, (115, 115): 0.0013758332620172988, (3, 3): 0.0014041703982190146, (54, 13): 0.00048837561185372967, (79, 52): 0.001183833447680738, (60, 60): 0.00086722562438467851, (48, 64): 0.0091840038427445942, (62, 62): 0.01340134055882653, (63, 63): 0.00030990919703870705, (72, 64): 0.000781957655145967, (3, 64): 0.0013068822485881278, (18, 18): 0.0034242384382712197, (27, 13): 0.00030887922417722737, (71, 42): 0.025362789559096297, (53, 49): 7.634187450514974e-08, (104, 104): 0.031047353518492301, (110, 110): 0.038017493095526496, (101, 101): 0.00043598844519924504, (42, 42): 0.0013937736025219081, (111, 111): 0.0011636981988944345, (100, 100): 0.14932404686451212, (95, 95): 0.098246724876328073, (54, 64): 0.0082366608814708198, (72, 72): 0.00040500580938876624, (82, 82): 0.0088787126866682903, (10, 64): 0.0039956573688116415, (47, 47): 0.0027789969151011938, (77, 64): 0.00096664081476323599}
	assert abs(diffab[(103, 20)] - 0.00037323542437035702) < .00001


def test_artificial_data():
	pass


def run_tests():
	test_real_data()
	test_artificial_data()
