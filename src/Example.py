import EMDUnifrac as EMDU
import numpy as np
import matplotlib as mpl
mpl.use('TkAgg')

# Run the tests to make sure everything is working
EMDU.run_tests()

# Load in the taxonomic tree
tree_str = "(Acidobacteria:0.03031,Actinobacteria:0.01878,Bacteroidetes:0.00530,Chlorobi:0.01402,Fusobacteria:0.10722," \
			"Lentisphaerae:0.04241,Proteobacteria:0.01667,Spirochaetes:0.03298,Synergistetes:0.03566,Tenericutes:0.00230," \
			"Verrucomicrobia:0.03222)Bacteria;"

# names of phyla (for easy creation of environments)
phyla = ["Acidobacteria", "Actinobacteria", "Bacteroidetes", "Chlorobi", "Firmicutes", "Fusobacteria", "Lentisphaerae",
		"Proteobacteria", "SR1", "Spirochaetes", "Synergistetes", "TM7", "Tenericutes", "Verrucomicrobia"]

# Abundances for healthy samples (via: healthy = np.loadtxt("/Users/dkoslicki/Dropbox/EMDUnifrac/RealData/HealthyTotals.txt"))
healthy = np.array([5.29000000e-04, 1.19146400e+00, 2.46148900e+01,
					5.29000000e-04, 1.85075170e+01, 1.49653000e-01,
					3.41100000e-03, 4.30437200e+00, 0.00000000e+00,
					1.77374000e-01, 8.38000000e-03, 0.00000000e+00,
					4.18820000e-02, 0.00000000e+00])

# Abundances for healthy samples (via: UC = np.loadtxt("/Users/dkoslicki/Dropbox/EMDUnifrac/RealData/UCTotals.txt"))
UC = np.array([0.00000000e+00, 5.64607000e-01, 6.00750200e+00,
				0.00000000e+00, 7.68366300e+00, 5.19120000e-02,
				8.80000000e-04, 1.64073200e+00, 1.57500000e-03,
				1.57500000e-03, 0.00000000e+00, 0.00000000e+00,
				4.75550000e-02, 0.00000000e+00])

# Create the environments
envs = dict()
for i in range(len(phyla)):
	envs[phyla[i]] = dict()
	envs[phyla[i]]["healthy"] = healthy[i]
	envs[phyla[i]]["UC"] = UC[i]

# Run EMDUnifrac (weighted, with flow) and plot the results
(Tint, lint, nodes_in_order) = EMDU.parse_tree(tree_str)
(nodes_weighted, samples) = EMDU.parse_envs(envs, nodes_in_order)
(Z, F, diffab) = EMDU.EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, nodes_weighted['healthy'], nodes_weighted['UC'])
EMDU.plot_diffab(nodes_in_order, diffab, 'healthy', 'Ulcerative Colitis')
#plt.savefig('/Users/dkoslicki/Dropbox/EMDUnifrac/Writeup/FullPaper/DiffAbundNew.png', bbox_inches='tight')


