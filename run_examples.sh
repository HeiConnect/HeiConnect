#!/bin/sh

# path to instances
graph_path="graphs/misc"
cactus_path="graphs/misc"
g="cycle-100.graph"
c="cycle-100.xml"

# compute the cactus graph using VieCut:
./extern/VieCut/build/mincut ${graph_path}/${g} -c ${cactus_path}/${c}

# run MSTConnect 
echo "Running MSTConnect"
./deploy/solver -g ${graph_path}/${g} -c ${cactus_path}/${c} -a mst-connect -o augmented_graph.graph

# run MSTConnect+LS
echo "Running MSTConnect + LS"
./deploy/solver -g ${graph_path}/${g} -c ${cactus_path}/${c} -a mst-connect-ls --depth 3 -o augmented_graph.graph

# run GreedyWeightCoverage 
echo "Running GreedyWeightCoverage"
./deploy/solver -g ${graph_path}/${g} -c ${cactus_path}/${c} -a gwc -o augmented_graph.graph

# run eILP 
echo "Running eILP"
# ./deploy/solver -g ${graph_path}/${g} -c ${cactus_path}/${c} -a eilp -o augmented_graph.graph
