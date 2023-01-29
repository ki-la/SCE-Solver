# SCE-Solver
A solver for the `Split Cluster Editing` problem introduced by Bruckner, Hüffner, and Komusiewicz [[BHK15]](https://link.springer.com/article/10.1186/s13015-015-0043-7), which I developed for my master thesis.


## Installation

The solver is written in [Rust](https://www.rust-lang.org/) and can be built with `cargo` by navigating into the `sce-solver` folder and excuting:
```
cargo build --release
```


### Libraries

The solver utilizes the following external libraries:

* [`VieCut`](https://github.com/VieCut/VieCut) is used in a lower bound to quickly compute the minimum cut.

* [`autocxx`](https://crates.io/crates/autocxx) is used to execute the C++ code from the `VieCut` library.

* [`rplex`](https://github.com/emallson/rplex) is used as an interface to [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer) in order to solve the problem via its ILP formulation. Hence, CPLEX has to be installed and the path variable CPLEX_LIB has to be set accordingly.


* [`petgraph`](https://crates.io/crates/petgraph) is used to represent and modify graphs as well as execute algorithms on them.

* further libraries from `crates.io`: `rand`, `indexmap`, `ahash`

## Test Dataset

The `test_instances` consist of 

* protein similarity graphs published [here](https://bio.informatik.uni-jena.de/data/#cluster_editing_data) (the originally weighted graphs were transformed to unweighted graphs by interpreting the 33/50/66% of edges with the highest weight as edges and the remaining ones as non-edges)

* and random graphs generated with the Erdős–Rényi model and the stochastic block model.

## Usage

The problem can be solved on _one_ graph by building the project as described above and executing:
```
./target/release/sce-solver path/to/graph/file config.txt 
```
The solver can be configured (in `config.txt`) to use different solution approaches (the main solver is a search-tree algorithm, but also an ILP solver using CPLEX and a heuristic solver can be used).


The experiments from my master thesis can be reproduced by executing:
``` 
bash ./run_all.sh
```
The script runs the solver on all graphs in `test_instances/final` for the search-tree solver and ILP solver once with and once without applying data reduction. The results are written to `results/final` and can be plotted by executing:
``` 
python3 evaluation.py
```

