# EMDUnifrac
EMDUnifrac is a computational method that computes the Unifrac [1,2,3] distance, but also returns information about which organisms are differentially abundant. This information comes in the form of a *flow* which is a matrix that shows exactly which organism abundances need to be moved where in the computation of Unifrac. A summary of this information is the *differential abundance* vector which shows which organisms are over/under expressed in each sample.

This method utilizes the Earth Mover's distance and details about the algorithm (and proof of correctness) are contained in the [preprint](http://biorxiv.org/content/early/2016/11/11/087171).

## Requirements ##
+ [dendropy](http://www.dendropy.org/)
+ numpy 
+ matplotlib (for plotting)
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
import EMDUnifrac as EMDU
tree_str = '((B:0.1,C:0.2)A:0.3);'  # input taxonomic tree
(Tint, lint, nodes_in_order) = EMDU.parse_tree(tree_str)  # parse the tree, getting nodes (Tint), edge lengths (lint), and node names (nodes_in_order)
# Create a toy environment
envs = {
	'C':{'sample1':1,'sample2':0},
	'B':{'sample1':1,'sample2':1},
	'A':{'sample1':0,'sample2':0},
	'temp0':{'sample1':0,'sample2':1}}  # temp0 is the root node, not named in Newick format, but included in nodes_in_order
(envs_prob_dict, samples) = EMDU.parse_envs(envs, nodes_in_order)  # Parse the environments.
(Z , F, diffab) = EMDU.EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, envs_prob_dict['sample1'], envs_prob_dict['sample2'])  # Run the weighted version of EMDUnifrac that returns the flow
# Check to make sure results make sense
assert Z == 0.25  # This is the Unifrac distance
assert F[(0,3)] == 0.5  # F is the flow and is in a sparse matrix format: a dictionary with tuple keys using elements of Tint and values T[(i, j)] equal to amount of abundance moved from organism nodes_in_order(i) to nodes_in_order(j)
assert F[(1,1)] == 0.5
assert sum(F.values()) == 1
assert diffab == {(2, 3): 0.14999999999999999, (0, 2): 0.10000000000000001}  # diffab is the differential abundance vector, also in a sparse matrix format: a dictionary with tuple keys using elements of Tint and values diffab[(i, j)] equal to the signed difference of abundance between the two samples restricted to the sub-tree defined by nodes_in_order(i) and weighted by the edge length lint[(i,j)].
```

An example ([Example.py](https://github.com/dkoslicki/EMDUnifrac/blob/master/src/Example.py)) is included using real biological data (restricted to the phylum 
level for simplicity).

### Description of formats and syntax ###
#### parse_tree ####
`(Tint, lint, nodes_in_order) = parse_tree(tree_str)`.

This function will parse a Newick tree string and return the dictionary of ancestors `Tint`.
`Tint` indexes the nodes by integers, `Tint[i] = j` means `j` is the ancestor of `i`.
`lint` is a dictionary returning branch lengths: `lint[(i,j)]` is the weight of the edge connecting `i` and `j`.
`nodes_in_order` is a list of the nodes in the input `tree_str` such that `Tint[i]=j` means `nodes_in_order[j]` is an ancestor of `nodes_in_order[i]`. Nodes are labeled from the leaves up.

#### parse_tree_file ####
`(Tint, lint, nodes_in_order) = parse_tree(tree_str_file)`.

Same as `parse_tree` but for reading in a Newick tree from a file `tree_str_file` instead of string.

#### parse_envs ####
`(envs_prob_dict, samples) = parse_envs(envs, nodes_in_order)`.

This function takes an environment `envs` and the list of nodes `nodes_in_order` and will return a dictionary `envs_prob_dict`
with keys given by samples. `envs_prob_dict[samples[i]]` is a probability vector on the basis `nodes_in_order` denoting abundances on the taxonomic tree for `samples[i]`.

The input data structure `envs` is a dictionary of dictionaries. The keys are elements of `nodes_in_order` and the values are dictionaries with *exactly two keys* `samples[i]` and `samples[j]`.
The values are the raw (or normalized) abundance of `samples[i]` assigned to a tree node given in `nodes_in_order`.

For example, the following environment:
```python
tree_str = '((B:0.1,C:0.2)A:0.3);'
envs = {
	'C':{'sample1':1,'sample2':0},
	'B':{'sample1':1,'sample2':1},
	'A':{'sample1':0,'sample2':0},
	'temp0':{'sample1':0,'sample2':1}}
```
represents two samples `sample1` and `sample2` where in `sample1`, one read respectively has been classified to the tree nodes `C` and `B`.
In `sample2`, one read respectively has been classified to the tree node `B` and the root node `temp0` (representing an unclassified read). This root node
is not contained in the Newick string `tree_str` but is automatically created by `parse_tree` and returned in `nodes_in_order`.

#### EMDUnifrac_weighted_flow ####
Weighted version of Unifrac that returns the flow.

`(Z, F, diffab) = EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, P, Q)`.

This function takes the ancestor dictionary `Tint`, the lengths dictionary `lint`, the basis `nodes_in_order`
and two probability vectors `P` and `Q` (typically `P = envs_prob_dict[samples[i]]`, `Q = envs_prob_dict[samples[j]])`.
Returns the weighted Unifrac distance `Z` (a scalar), the flow `F`, and the differential abundance vector `diffab`.
The flow `F` is a dictionary with tuple keys of the form `(i, j)` for `i,j` in `Tint` where `F[(i, j)] == num` means that in the calculation of the
Unifrac distance, a total mass of `num` was moved from the node `nodes_in_order[i]` to the node `nodes_in_order[j]`.
The differential abundance vector `diffab` is a dictionary with tuple keys using elements of `Tint` and values
`diffab[(i, j)]` equal to the signed difference of abundance between the two samples restricted to the sub-tree
defined by `nodes_in_order[i]` and weighted by the edge length `lint[(i, j)]`.

#### EMDUnifrac_weighted ####
Weighted version of Unifrac, does not return the flow (faster execution time than when returning the flow).

`(Z, diffab) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)`.

This function takes the ancestor dictionary `Tint`, the lengths dictionary `lint`, the basis `nodes_in_order`
and two probability vectors `P` and `Q` (typically `P = envs_prob_dict[samples[i]]`, `Q = envs_prob_dict[samples[j]])`.
Returns the weighted Unifrac distance `Z` (a scalar), and the differential abundance vector `diffab`.
The differential abundance vector `diffab` is a dictionary with tuple keys using elements of `Tint` and values
`diffab[(i, j)]` equal to the signed difference of abundance between the two samples restricted to the sub-tree
defined by `nodes_in_order[i]` and weighted by the edge length `lint[(i, j)]`.

#### EMDUnifrac_unweighted_flow ####
Unweighted version of Unifrac, returns the flow.

`(Z, F, diffab) = EMDUnifrac_unweighted_flow(Tint, lint, nodes_in_order, P, Q)`.

This function takes the ancestor dictionary `Tint`, the lengths dictionary `lint`, the basis `nodes_in_order`
and two probability vectors `P` and `Q` (typically `P = envs_prob_dict[samples[i]]`, `Q = envs_prob_dict[samples[j]])`.
Returns the weighted Unifrac distance `Z` (a scalar), the flow `F`, and the differential abundance vector `diffab`.
The flow `F` is a dictionary with tuple keys of the form `(i, j)` for `i,j` in `Tint` where `F[(i, j)] == num` means that in the calculation of the
Unifrac distance, a total mass of `num` was moved from the node `nodes_in_order[i]` to the node `nodes_in_order[j]`.
The differential abundance vector `diffab` is a dictionary with tuple keys using elements of `Tint` and values
`diffab[(i, j)]` equal to the signed difference of abundance between the two samples restricted to the sub-tree
defined by `nodes_in_order[i]` and weighted by the edge length `lint[(i, j)]`.

#### EMDUnifrac_unweighted ####
Unweighted version of Unifrac, does not return the flow (faster execution time than when returning the flow).

`(Z, diffab) = EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)`.

This function takes the ancestor dictionary `Tint`, the lengths dictionary `lint`, the basis `nodes_in_order`
and two probability vectors `P` and `Q` (typically `P = envs_prob_dict[samples[i]]`, `Q = envs_prob_dict[samples[j]])`.
Returns the weighted Unifrac distance `Z` (a scalar), and the differential abundance vector `diffab`.
The differential abundance vector `diffab` is a dictionary with tuple keys using elements of `Tint` and values
`diffab[(i, j)]` equal to the signed difference of abundance between the two samples restricted to the sub-tree
defined by `nodes_in_order[i]` and weighted by the edge length `lint[(i, j)]`.


## Authors ##
David Koslicki <david.koslicki@math.oregonstate.edu>

Jason McClelland <mcclellj@science.oregonstate.edu>


## License ##
This project is released under the GPL-3 License. Please view the [LICENSE](https://github.com/dkoslicki/EMDUnifrac/blob/master/LICENSE) file for more details.

## References ##
*To do*