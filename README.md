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





## Usage ##
David Koslicki [david.koslicki@math.oregonstate.edu](david.koslicki@math.oregonstate.edu)
Jason McClelland [mcclellj@science.oregonstate.edu](mcclellj@science.oregonstate.edu)
