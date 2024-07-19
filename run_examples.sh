#! /bin/bash
#path to instances
graph_path="graphs/misc"
cactus_path="graphs/misc"
g="cycle-100.graph"
c="cycle-100.xml"

# compute the cactus graph using VieCut:
./extern/VieCut/build/mincut ${graph_path}/${g} -c ${cactus_path}/${c}

# run MSTConnect 
echo "Running MSTConnect"
./deploy/solver -g ${graph_path}/${g} -c ${cactus_path}/${c} -a mst-flow -o augmented_graph.graph

# run MSTConnect+LS
echo "Running MSTConnect + LS"
./deploy/solver -g ${graph_path}/${g} -c ${cactus_path}/${c} -a mst-ls -o augmented_graph.graph

# run GreedyWeightCoverage 
echo "Running GreedyWeightCoverage"
./deploy/solver -g ${graph_path}/${g} -c ${cactus_path}/${c} -a dynamic -o augmented_graph.graph --sampling=1

# run eILP 
echo "Running eILP"
# ./deploy/solver -g ${graph_path}/${g} -c ${cactus_path}/${c} -a ilp -o augmented_graph.graph
