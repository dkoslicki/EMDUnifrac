# EMDUnifrac
EMDUnifrac is a computational method that computes the Unifrac [1,2,3] distance, but also returns information about which organisms are differentially abundant. This information comes in the form of a *flow* which is a matrix that shows exactly which organism abundances need to be moved where in the computation of Unifrac. A summary of this information is the *differential abundance* vector which shows which organisms are over/under expressed in each sample.

This method utilizes the Earth Mover's distance and details about the algorithm (and proof of correctness) are contained in the [preprint](http://biorxiv.org/content/early/2016/11/11/087171).

## Requirements ##
+ [dendropy](http://www.dendropy.org/)
+ numpy 
+ [ete2](http://etetoolkit.org/download/) (only for the speed comparison to Unifrac)
+ [PyCogent](https://github.com/pycogent/pycogent) (only for the speed comparison to Unifrac)

## Usage ##
Efforts have been made to make the API similar to that of the PyCogent implementation of FastUnifrac.

The basic workflow is:
+ Parse an input taxonomic treen (in Newick format) using `(Tint, lint, nodes_in_order) = parse_tree(tree_str)`.
+ Parse input taxonomic classifications (called "environments") using `(nodes_weighted, samples) = parse_envs(input_environments, nodes_in_order)`.
+ Run a weighted/unweighted version of EMDUnifrac with/without the flow returned. For example: `(unifrac_value, flow, differential_abundance_vector) = EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, input_environments['sample1'], input_environments['sample2'])`.

For the impatient, here is a minimal working example (for more, see the tests contained in the [source](https://github.com/dkoslicki/EMDUnifrac/blob/master/src/EMDUnifrac.py)):
```python
tree_str = '((B:0.1,C:0.2)A:0.3);'  # input taxonomic tree
(Tint, lint, nodes_in_order) = parse_tree(tree_str)  # parse the tree, getting nodes (Tint), edge lengths (lint), and node names (nodes_in_order)
# Create a toy environment
nodes_samples = {
	'C':{'sample1':1,'sample2':0},
	'B':{'sample1':1,'sample2':1},
	'A':{'sample1':0,'sample2':0},
	'temp0':{'sample1':0,'sample2':1}}  # temp0 is the root node, not named in Newick format, but included in nodes_in_order
(nodes_weighted, samples) = parse_envs(nodes_samples, nodes_in_order)  # Parse the environments.
(Z , F, diffab) = EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, nodes_weighted['sample1'], nodes_weighted['sample2'])  # Run the weighted version of EMDUnifrac that returns the flow
# Check to make sure results make sense
assert Z == 0.25  # This is the Unifrac distance
assert F[(0,3)] == 0.5  # F is the flow and is in a sparse matrix format: a dictionary with tuple keys using elements of Tint and values T[(i, j)] equal to amount of abundance moved from organism nodes_in_order(i) to nodes_in_order(j)
assert F[(1,1)] == 0.5
assert sum(F.values()) == 1
assert diffab == {(2, 3): 0.14999999999999999, (0, 2): 0.10000000000000001}  # diffab is the differential abundance vector, also in a sparse matrix format: a dictionary with tuple keys using elements of Tint and values diffab[(i, j)] equal to the signed difference of abundance between the two samples restricted to the sub-tree defined by nodes_in_order(i) and weighted by the edge length lint[(i,j)].
```

## Authors ##
David Koslicki <david.koslicki@math.oregonstate.edu>
Jason McClelland <mcclellj@science.oregonstate.edu>


## License ##
This project is released under the GPL-3 License. Please view the [LICENSE](https://github.com/dkoslicki/EMDUnifrac/blob/master/LICENSE) file for more details.