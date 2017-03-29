# This is an example file of how to use the CAMI parser in parallel
import EMDUnifrac as EMDU
import CAMI
import numpy as np
import argparse
import sys
import os

__author__ = 'David Koslicki (dmkoslicki@gmail.com, david.koslicki@math.oregonstate.edu)'
__version__ = '1.0.0'
__date__ = '28 Mar 2017'


def read_params(args):
	parser = argparse.ArgumentParser(description='')
	arg = parser.add_argument
	arg('--input', metavar='files_file', type=str, required=True,
		default=None, help="File of CAMI profile files to use")
	arg('--output', metavar='output_file', required=True, default=None, type=str,
		help="Output file (you should have this end in .csv as it is a matrix)")
	return vars(parser.parse_args())


def make_dist_mat(files_file, output):
	files = []
	fid = open(files_file, 'r')
	for line in fid.readlines():
		files.append(line.strip())

	fid.close()

	D = np.zeros((len(files), len(files)))

	for i in xrange(len(files)):
		for j in xrange(i+1, len(files)):
			in_file_1 = files[i]
			in_file_2 = files[j]
			if os.path.exists(in_file_1) and os.path.exists(in_file_2):
				tax_ids_1, tax_path_1, weights_1 = CAMI.CAMI_read_taxonomy_file(in_file_1)
				tax_ids_2, tax_path_2, weights_2 = CAMI.CAMI_read_taxonomy_file(in_file_2)
				Tint, lint, nodes_in_order, nodes_to_index = CAMI.CAMI_form_graph(tax_path_1, tax_path_2)
				perc1 = CAMI.CAMI_make_percentages(weights_1, nodes_to_index)
				prob1 = CAMI.CAMI_get_probability_distribution(nodes_in_order, Tint, perc1)
				perc2 = CAMI.CAMI_make_percentages(weights_2, nodes_to_index)
				prob2 = CAMI.CAMI_get_probability_distribution(nodes_in_order, Tint, perc2)
				#(Z, F, diffab) = EMDU.EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, prob1, prob2)
				(Z, diffab) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, prob1, prob2)
				D[i, j] = Z
				D[j, i] = Z
	np.savetxt(output, D, delimiter=',', newline='\n')
	#print(D)


if __name__ == '__main__':
	par = read_params(sys.argv)
	make_dist_mat(par['input'], par['output'])

	# ACK!! Looks like EMDUnifrac_weighted is not returning correct results
	# EMDUnifrac_weighted_flow for test files 1 and 2 gives ~0.167 while EMDUnifrac_weighted gives ~2.056

	####################
	#in_file_1 = '/home/dkoslicki/Dropbox/Repositories/EMDUnifrac/data/test0.profile'
	#in_file_2 = '/home/dkoslicki/Dropbox/Repositories/EMDUnifrac/data/test2.profile'
	#tax_ids_1, tax_path_1, weights_1 = CAMI.CAMI_read_taxonomy_file(in_file_1)
	#tax_ids_2, tax_path_2, weights_2 = CAMI.CAMI_read_taxonomy_file(in_file_2)
	#Tint, lint, nodes_in_order, nodes_to_index = CAMI.CAMI_form_graph(tax_path_1, tax_path_2)
	#perc1 = CAMI.CAMI_make_percentages(weights_1, nodes_to_index)
	#prob1 = CAMI.CAMI_get_probability_distribution(nodes_in_order, Tint, perc1)
	#perc2 = CAMI.CAMI_make_percentages(weights_2, nodes_to_index)
	#prob2 = CAMI.CAMI_get_probability_distribution(nodes_in_order, Tint, perc2)
	#(Z, F, diffab) = EMDU.EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, prob1, prob2)
	#(Z2, diffab2) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, prob1, prob2)

#in_file_1 = '/home/dkoslicki/Dropbox/Repositories/EMDUnifrac/data/test3.profile'
#in_file_2 = '/home/dkoslicki/Dropbox/Repositories/EMDUnifrac/data/test4.profile'
#tax_ids_1, tax_path_1, weights_1 = CAMI.CAMI_read_taxonomy_file(in_file_1)
#tax_ids_2, tax_path_2, weights_2 = CAMI.CAMI_read_taxonomy_file(in_file_2)
#Tint, lint, nodes_in_order, nodes_to_index = CAMI.CAMI_form_graph(tax_path_1, tax_path_2)
#perc1 = CAMI.CAMI_make_percentages(weights_1, nodes_to_index)
#prob1 = CAMI.CAMI_get_probability_distribution(nodes_in_order, Tint, perc1)
#perc2 = CAMI.CAMI_make_percentages(weights_2, nodes_to_index)
#prob2 = CAMI.CAMI_get_probability_distribution(nodes_in_order, Tint, perc2)
#(Z, diffab) = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, prob1, prob2)