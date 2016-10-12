# FastRFS

## Running FastRFS 

Download the FastRFS executable from the Releases page (https://github.com/pranjalv123/FastRFS/releases). 

The simplest way to run FastRFS is like

    FastRFS -i sourcetrees.newick
    
where sourcetrees.newick is the path to a file containing your source trees in newick format.

Currently, this file must have the following properties:

 1. The path to the file must not have spaces
 2. The newick trees must not have branch lengths
 3. The newick trees must be compatible with ASTRAL-II

FastRFS will output the score (total RF distance to the input trees) and the output tree to the command line.

To save this tree to a file, run

    FastRFS -i sourcetrees.newick -o outputfile

To add extra trees to expand the search space (the set X) of your analysis, run
    
    FastRFS -i sourcetrees.newick -e extratrees.newick
    
The extratrees file should be in the same format as sourcetrees.newick, except that it is OK if it has branch lengths.


## Building FastRFS

You only need to do this if the FastRFS executable doesn't run on your system 
or you're making changes to the source code.

### Dependencies

You need to have the following installed on your computer (Ubuntu
package names or URLs are listed):

 - Boost (libboost-all-dev)
 - ASTRAL (https://github.com/smirarab/ASTRAL/)
 
### Building

Clone and cd into the repository, then do:

```
mkdir build/
cd build/
cmake ../src
make
```
Then copy the ASTRAL installation into the directory with FastRFS
