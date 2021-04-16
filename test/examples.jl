# using Revise
# using Pkg; Pkg.activate(".")
using HyperModularity
using Test

## Reading in a dataset
dataset = "contact-high-school-classes"
maxsize = 5
minsize = 2
return_labels = true
H, L = read_hypergraph_data(dataset,maxsize,minsize,return_labels)


## Run graph Louvain on reduced graph

# weighted = true means that when clique expandin hyperedge e
# weight all edges by 1/(|e|-1)
# weighted = false means they are weighted by 1, but the
# resulting graph is still integer weighted
weighted = true

# binary = true means that after clique expansion, all edges
# are unweighted. If true, "weighted" flag doesn't matter
binary = false
A = CliqueExpansion(H,weighted,binary)

# Run modularity on clique expansion
gamma = 1.0
Z = CliqueExpansionModularity(H,gamma;weighted=true,randflag=false,binary=false,clusterpenalty=0.0,maxits=1000)

# Compute the MLE for the resolution parameter based on Z
(ωᵢ, ωₒ) = computeDyadicResolutionParameter(H, Z; mode = 0)
γ = (ωᵢ - ωₒ)/(log(ωᵢ) - log(ωₒ))

# log-likelihood
loglike = dyadicLogLikelihood(H, Z, ωᵢ, ωₒ)

# Above steps can be iteratively repeated to reach global max for log-likelihood
# Good place to test log likelihood computations
