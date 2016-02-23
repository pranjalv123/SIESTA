# The DP algorithm

## Dependencies

You need to have the following installed on your computer (Ubuntu
package names or URLs are listed):

 - Boost (libboost-all-dev)
 - Python (python-dev) with dendropy installed
 - ASTRAL (https://github.com/smirarab/ASTRAL/)
 - A compiler that supports OpenMP - on OS X you might have to do http://stackoverflow.com/questions/35134681/installing-openmp-on-mac-os-x-10-11
 
 
## Building

Clone and cd into the repository, then do:

```
mkdir build/
cd build/
cmake ../src
make
```

## Running

### Common arguments

The three parameters you will need for almost any run of the software
are:

```
-a/--astral path/to/astral.jar
-g/--genetrees genetrees-file
-c/--criterion criterion
-o/--output output
--maximize/--minimize
```

The astral jar file can be from any recent version of ASTRAL; the
genetrees file should be a list of trees in newick format, and the
criterion is the scoring function used. One of --maximize and
--minimize should be provided to determine if high scores are good or
bad. -o is optional and provides a location to output the computed
tree.


#### Criteria

Valid criteria are:

```
BryantSteel
DP
RF
Python
```


##### Quartet support

BryantSteel and DP both optimize induced quartet weights. These both
expect an option
```
-q quartetsfile
```

where quartetsfile is a file with quartets either in newick format,
for example

```
((a,b),(c,d)); 4.3
```

or in wQMC format, for example

```
a,b|c,d:4.3
```

##### RF Supertrees

RF optimizes the SumRF distance to the source trees. It expects no
additional options, and uses the trees given in -g as the source trees.

##### Python scorers

-c Python allows criterion scorers written in Python. If the criterion
 scorer is provided in a file called my_scorer.py in folder
 /path/to/scorers, PYTHONPATH must contain /path/to/scorers, which can
 be done with

 ```
 export PYTHONPATH=PYTHONPATH:/path/to/scorers
 ```

 Then, the name of the python scorer should be given as the name of
 the file *without* the .py extension, e.g.

 ```
 -c Python -p my_scorer
 ```

Two python scorers are provided: quintet_scorer and
mulrf_scorer. quintet_scorer uses rooted quintet weights, provided as
`-q quintetfile` in format

```
(a, ((b, c), (d, e)): 4.3
```

mulrf_scorer is an extension of the RF supertree scorer that supports
trees with polytomies, and uses the `-g genetrees` argument that is
provided to get the search space.


### Scoring trees

Trees can be scored with the options

```
-s/--score tree-to-score 
--rootedscore rooted-tree-to-score
```

In this case, the tree will be scored according to the provided
criterion, and the score will be output to the file given in `-o
output`.


### Controlling the search space

A number of options can be used to control the search space.  The
first is `-x/--exact`, which runs the optimization algorithm in exact
mode. This is only practical for a small number of taxa (up to 15 or
so, depending on the optimization criterion).

Next, `-e/--extragenetrees moregenetrees` can be used to add
bipartitions from an additional set of trees. These trees are not
expanded with ASTRAL's heuristics, but if they are incomplete or have
polytomies they are completed and resolved using ASTRAL's default
methods.

If `--extraextra` is also provided as an option, the ASTRAL expansion
of the extra gene trees, as well as the ASTRAL expansion of the
combination of the extra gene trees and the original gene trees.

If `--limited` is provided as an option, the input gene tree set will
not be expanded by ASTRAL. However, incomplete trees or trees with
polytomies will still be completed or resolved.

Finally, `-X cladefile` can be used to add clades/clusters to the
search space. These should be comma-delimited lists surrounded by
curly braces, e.g.
```
{a,b,c,d,e,f,g}
```

## Writing your own scoring functions

You can also write your own scoring functions to solve other
problems. It's possible to do this in C++ or in Python.


### Python

The DP code interfaces with DendroPy so you can use the DendroPy tools
you're already used to using. You should create a Python file with
functions `init(ts, opts)` and `score(a1, a2, rest)`.

`init(ts, opts)` gets a `dendropy.TaxonNamespace` object `ts` and a
list of command-line arguments `opts`. You can use `argparse` to parse
the command line options, and you should pass `ts` as a
taxon_namespace argument to dendropy constructors you call. You should
not add more taxa to the taxon namespace. In particular, you need to
be careful when reading newick files where taxa have underscores - you
should call the readers with `preserve_underscores=True`.

`score(a1, a2, rest)` gets three bitmasks corresponding to a
tripartition, where `a1` and `a2` are below the node and `rest` is
above the node. You must return a `float` (returning an `int` will
result in undefined behavior) corresponding to the weight of the
tripartition.

In addition, ou may implement the method `adjust(n)`. This is used
once the algorithm is finished running, to modify the score of the
final tree before printing and saving.


### C++

Writing a scorer in C++ is more involved, but is typically
substantially faster than Python.

Create a class that inherits from `TripartitionScorer`. Within the
class, you must call DEC_SCORER(scorername), where scorername is the
name of your class. You must also call DEF_SCORER(scorername) outside
of your main class declaration (in a .cpp file or wherever you define
your functions).


The constructor of your class must take a single argument `TaxonSet&
ts`. You should use this when reading data, constructing bitmasks,
etc. 

You must implement the method `virtual double score(const
Tripartition& t)`, which should return the score corresponding to
`t`. 

You may implement the virtual method `virtual
double adjust_final_score(double score)`. This is used once the
algorithm is finished running, to modify the score before printing and
saving.

#### Utility classes

##### TaxonSet

This maps taxon labels (strings) to indices (integers). If you have
`TaxonSet t`, `t['foo']` will give you the index of taxon 'foo', and
`t[4]` will give you the label of the taxon with index 4. TaxonSets
given to tripartition scorers are frozen, meaning that trying to
access a taxon not in the set will cause the program to print an error
and fail.

##### Tripartition

A tripartition contains three Clades - `a1, a2, rest`. `a1` and `a2`
are the partitions below the root, and `rest` is the partition above
the root. Note that the constructor ot Tripartition takes `TaxonSet& ts,
Clade& clade, Clade& subclade`. `a1` corresponds to `clade -
subclade`, `a2` corresponds to `subclade`, and `rest` corresponds to
`ts - clade`.

##### Clade

A clade represents a subset of taxa. It stores its data in
`BitVectorFixed taxa` and supports the following operations:

 - `bool Clade::contains(const Clade& other), bool
Clade::contains(const Taxon taxon)`: returns `true` if `other` is a
strict subset of this clade, or `taxon` is in this clade.

 - `Clade Clade::overlap (const Clade& other)`: returns a clade
   containing the taxa that are in both this clade and `other`.

 - `void add(const Taxon taxon)`: Add `taxon` to this clade

 - `Clade complement()`: Return a clade containing only the taxa not
   in this clade.

 - `Clade minus(const Clade& other)`: Return a clade containing only
   the taxa in this clade but not in `other`.

 - `int size()`: Return the number of taxa in this clade


The clade class also supports iteration with `begin()` and `end()`,
and is hashable.

##### BitVectorFixed

The `BitVectorFixed` object is a datatype representing an binary
string of arbitrary length, but which has length fixed on
initialization. It supports the `~`, `|`, `^`, and `&` (bitwise NOT,
OR, XOR, and AND) operators. It also has functions `get(int i)`,
`set(int i)` and `unset(int i)`, which get, set and unset a particular
bit, as well as `ffs()`, which finds the lowest set bit, and
`popcount()`, which counts the number of bits set to
`1`. `BitVectorFixed` also supports (efficient) iteration over set
bits using the `begin()` and `end()` functions, and is hashable.

