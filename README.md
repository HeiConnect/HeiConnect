HeiConnect v1.0  [![Codacy Badge](https://app.codacy.com/project/badge/Grade/9d0d08ba6b2d42699ab74fe5f9697bb9)](https://www.codacy.com/gh/KaHIP/KaHIP/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=KaHIP/KaHIP&amp;utm_campaign=Badge_Grade)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
=====

# Weighted Connectivity Augmentation Algorithms
This is the code repository contains an exact ILP-based solver as well as several heuristic algorithms for the weighted connectivity augmentation problem.

## Building
The project can be build using CMake. Make sure that the [Gurobi Optimizer](https://www.gurobi.com/solutions/gurobi-optimizer) and the [Boost](https://boost.org) and [LEMON](https://lemon.cs.elte.hu/trac/lemon) graph libraries can be found on your system. Gurobi is used to solve ILPs and LPs, Boost is mainly used to compute maximum flows while LEMON is used to compute a maximum weighted matching due to bugs in the Boost graph library.
```sh
./build.sh
```

Verbose logging can be disabled at compile time with the CMake flag `-DENABLE_VERBOSE_LOG=OFF` for code optimization.

## Running
The main binary is `solver`, i.e. `deploy/solver`.

Required arguments are
- `--graph <path>`: The path to the original graph file (in the [METIS graph format](http://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html))
- `--cactus <path>`: The path to the cactus graph file (GraphML xml format, i.e. created by [VieCut](https://github.com/VieCut/VieCut), `mincut <graph> cactus -t <cactus>`)
- `--algorithm <algorithm>`: the algorithm to use

Optional arguments include:
- `--verbose`: Enable verbose logging
- `--output`: Specify an output file where the augmented graph should be written to

Examples of how to use our different algoirhtms can be found in `run_examples.sh`.

### Algorithms
|Name|Status|Description|Parameters|Limitations|
|----|------|-----------|----------|-----------|
|`greedy`|OK|Weight + Cut Heuristic (naive implementation)|none||
|`greedy-strong`|OK|Global Weight Heuristic (naive implementation)|none||
|`heuristic`|Deprecated / replaced by `dynamic`|Weight-Heuristic (naive implementation)|none||
|`dynamic`|OK|Weight-Coverage Heuristic with dynamic cactus data structure|`--sampling=<0,1>` (`1` enables sampling)||
|`dynamic-bounded`|OK|Weight-Coverage Heuristic with dynamic cactus and bounds|none|Weights must be in `[0,1]`|
|`heuristic-sampling`|OK|Naive implementation of Weight-Coverage Heuristic, but with sampling|`--sampling=<int>` (number of link buckets, `0` is `sqrt n`)||
|`ilp`|OK|Optimal ILP|`--use-initial` (Use `mst-flow` as initial solution)||
|`apx2lp`|OK|LP-based 2-approximation|none||
|`mst`|Deprecated / replaced by `mst-flow`|MST Algorithm (naive implementation)|none||
|`mst-ilp`|OK|Find optimum solution using only links of MSTs|`--trees=<int>` (number of minimum spanning trees to use)||
|`apx1ln2e`|OK||`--epsilon=<float>` (set epsilon, default 0.15)||
|`apx1.5e`|OK||`--epsilon=<float>` (set epsilon, default 0.15)||
|`mst-flow`|OK|MST Algorithm|none||
|`full-mst`|OK|Use whole MST as solution|none||
|`mst-ls`|Deprecated|Local Search (naive implementation)|`--depth=<int>` (max. length of alternating paths)||
|`mst-ls-flow`|OK|Local Search on MST algorithm solution|`--depth=<int>` (max. length of alternating paths), `--cache` (cache invalid alternating paths), `--trees=<int>` (number of MSTs to use)||
|`smc`|OK|SMC implementation using dynamic graph data structure|none||
|`fsm`|OK|FSM implementation using dynamic graph data structure|none||
|`hbd`|OK|HBD implementation using dynamic graph data structure|none||
|`mst-heuristic`|Failed|MST, but use weight coverage heuristic as weights. Significantly worse solutions|none||

### Graph instances
The generated graphs mentioned in the paper are included in `graphs/generated`. A small subset of the graphs from the 10th DIMACS Implementation Challenge is included as well, the complete set can be found at https://dimacs10.github.io/downloads.shtml.


### Licence
The program is licenced under MIT licence.
If you publish results using our algorithms, please acknowledge our work by quoting the following paper:
 *Fonseca, M., Großmann, E., Joos, F., Möller, T. and Schulz C. 2024. Engineering Weighted Connectivity Augmentation Algorithms. arXiv preprint [arXiv:1708.06127.](https://arxiv.org/abs/2402.07753)*


```bibtex
@inproceedings{DBLP:conf/wea/FarajGJM024,
  author       = {Marcelo Fonseca Faraj and
                  Ernestine Gro{\ss}mann and
                  Felix Joos and
                  Thomas M{\"{o}}ller and
                  Christian Schulz},
  editor       = {Leo Liberti},
  title        = {Engineering Weighted Connectivity Augmentation Algorithms},
  booktitle    = {22nd International Symposium on Experimental Algorithms, {SEA} 2024,
                  July 23-26, 2024, Vienna, Austria},
  series       = {LIPIcs},
  volume       = {301},
  pages        = {11:1--11:22},
  publisher    = {Schloss Dagstuhl - Leibniz-Zentrum f{\"{u}}r Informatik},
  year         = {2024},
  url          = {https://doi.org/10.4230/LIPIcs.SEA.2024.11},
  doi          = {10.4230/LIPICS.SEA.2024.11},
  timestamp    = {Fri, 12 Jul 2024 15:29:30 +0200},
  biburl       = {https://dblp.org/rec/conf/wea/FarajGJM024.bib},
  bibsource    = {dblp computer science bibliography, https://dblp.org}
}
```
