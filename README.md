[![Build Status](https://github.com/nveldt/HyperModularity/workflows/CI/badge.svg)](https://github.com/nveldt/HyperModularity/actions)
[![Build Status](https://travis-ci.com/nveldt/HyperModularity.svg?branch=master)](https://travis-ci.com/github/nveldt/HyperModularity)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/nveldt/HyperModularity?svg=true)](https://ci.appveyor.com/project/nveldt/hypermodularity)
[![Coverage](https://codecov.io/gh/nveldt/HyperModularity/branch/master/graph/badge.svg)](https://codecov.io/gh/nveldt/HyperModularity.jl)

# HyperModularity

This package contains code for hypergraph modularity clustering algorithms, based on the generative model described in

**Generative hypergraph clustering: from blockmodels to modularity**  
Philip S. Chodrow, Nate Veldt, Austin R. Benson  
[preprint](https://arxiv.org/abs/2101.09611)

## Package Installation
##### To install package
```
using Pkg
Pkg.add("HyperModularity")
using HyperModularity
```

## Data available
##### For a full list of all datasets:
```
hypermodularity_datasets()
```
##### Loading data example:
```
dataset = "contact-high-school-classes"
maxsize = 5	# max hyperedge size
minsize = 2	# min hyperedge size
return_labels = true
H, L = read_hypergraph_data(dataset,maxsize,minsize,return_labels)
```

## Some Examples

#### Run graph Louvain algorithm on clique expanded graph

```
gamma = 1.0 # default resolution parameter
Z_g = CliqueExpansionModularity(H,gamma) # see code for default parameters

# Compute MLE resolution parameter given clustering Z_g
(ωᵢ, ωₒ) = computeDyadicResolutionParameter(H, Z_g; mode = 0)
γ_mle = (ωᵢ - ωₒ)/(log(ωᵢ) - log(ωₒ))
loglike = dyadicLogLikelihood(H, Z, ωᵢ, ωₒ)
```

#### Run all-or-nothing Louvain algorithm

```
n = size(H,2)
Z_ = collect(1:n) # trivial clustering

# all or nothing aggregator: p -> [length(p) == 1, sum(p)]
# This gives a starter estimate for Ω, from a trivial clustering Z_
Ω = estimateΩEmpirically(H, Z_; aggregator = p -> [length(p) == 1, sum(p)])

Z = AON_Louvain(H,Ω)

# Alternatively, one can learn Ω from graph Louvain solution Z_g
Ω = estimateΩEmpirically(H, Z_g; aggregator = p -> [length(p) == 1, sum(p)])
Z = AON_Louvain(H,Ω)
```

#### Additional examples

See demo notebooks in demos folder for other examples of how to use the code.