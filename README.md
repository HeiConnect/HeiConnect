HeiConnect v1.0 [![Codacy Badge](https://app.codacy.com/project/badge/Grade/9d0d08ba6b2d42699ab74fe5f9697bb9)](https://www.codacy.com/gh/KaHIP/KaHIP/dashboard?utm_source=github.com&utm_medium=referral&utm_content=KaHIP/KaHIP&utm_campaign=Badge_Grade)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
=====

## Weighted Connectivity Augmentation Algorithms

Increasing the connectivity of a graph is a pivotal challenge in robust network design. The weighted connectivity augmentation problem is a common version of the problem that takes link costs into consideration. The problem is then to find a minimum cost subset of a given set of weighted links that increases the connectivity of a graph by one when the links are added to the edge set of the input instance. Here, we give a first implementation of recently discovered better-than-2 approximations. Furthermore, we propose three new heuristic and one exact approach. These include a greedy algorithm considering link costs and the number of unique cuts covered, an approach based on minimum spanning trees and a local search algorithm that may improve a given solution by swapping links of paths. Our exact approach uses an ILP formulation with efficient cut enumeration as well as a fast initialization routine. 

## Download 
You can download HeiConnect with the following command line:

```console
git clone --recursive git@github.com:HeiConnect/HeiConnect.git
```

## Building

The project can be build using CMake. Make sure that the [Gurobi Optimizer](https://www.gurobi.com/solutions/gurobi-optimizer) (tested >= 11.0.2) and the [Boost](https://boost.org) (tested >= 1.83) and [LEMON](https://lemon.cs.elte.hu/trac/lemon) (tested >= 1.3.1) graph libraries can be found on your system. Gurobi is used to solve ILPs and LPs, Boost is mainly used to compute maximum flows while LEMON is used to compute a maximum weighted matching due to bugs in the Boost graph library.

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
- `--links`: Specify a file containing one link per line (`source target weight`)
- `--output`: Specify an output file where the augmented graph should be written to

Examples of how to use our different algorithms can be found in `run_examples.sh`.

### Algorithms

| Name             | Description                                                                     | Parameters                                                                                                                               | Recommended for instances with |
| ---------------- | ------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------|
| `gwc`            | Greedy Weight-Coverage Heuristic using dynamic cactus data structure and bounds | `--sampling=<0,1>` (`1` enables sampling)                                                                                                | small link costs (generated)   | 
| `mst-connect`    | MST-Connect Algorithm                                                           | none                                                                                                                                     | large link costs               |
| `mst-connect-ls` | Local Search on MST-Connect Algorithm solution                                  | `--depth=<int>` (max. length of alternating paths), `--cache` (cache invalid alternating paths), `--trees=<int>` (number of MSTs to use) | large link costs               |
| `eilp`           | Optimal ILP                                                                     | `--use-initial` (Use `mst-connect` as initial solution)                                                                                  | small link costs (real world)  |
| `apx2lp`         | LP-based 2-approximation                                                        | none                                                                                                                                     |
| `apx1ln2e`       |                                                                                 | `--epsilon=<float>` (set epsilon, default 0.15)                                                                                          |
| `apx1.5e`        |                                                                                 | `--epsilon=<float>` (set epsilon, default 0.15)                                                                                          |
| `smc`            | SMC implementation using dynamic graph data structure                           | none                                                                                                                                     |
| `fsm`            | FSM implementation using dynamic graph data structure                           | none                                                                                                                                     |
| `hbd`            | HBD implementation using dynamic graph data structure                           | none                                                                                                                                     |

### Graph instances

The generated graphs mentioned in the paper are included in `graphs/generated`. A small subset of the graphs from the 10th DIMACS Implementation Challenge is included as well, the complete set can be found at https://dimacs10.github.io/downloads.shtml.

### Licence

The program is licenced under MIT licence.
If you publish results using our algorithms, please acknowledge our work by quoting the following paper:
_Fonseca, M., Großmann, E., Joos, F., Möller, T. and Schulz C. 2024. Engineering Weighted Connectivity Augmentation Algorithms. arXiv preprint [arXiv:1708.06127.](https://arxiv.org/abs/2402.07753)_

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
