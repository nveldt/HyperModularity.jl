module HyperModularity

import Combinatorics #
import Distributions #
import LinearAlgebra #
import Parameters #
import Random #
import SparseArrays #
import StatsBase #
import Optim #
import SpecialFunctions #

include("omega.jl")     # uses Combinatorics
include("HSBM.jl")      # uses Parameters, Distributions
include("dict_ops.jl")
include("utils.jl")
include("objectives.jl")
include("inference.jl") # uses Optim
include("analysis_helpers.jl")
include("vol.jl")       # uses StatsBase
include("louvain_utils.jl")
include("graph_louvain.jl") # uses Special Functions
include("read_data.jl")
include("AON_hyperlouvain.jl") # uses Random
include("general_hypergraph_louvain.jl")
include("simple_synthetic.jl")

export partitionize

# from analysis_helpers.jl
export downSampleEdges!
export subhypergraph
export removeEdges!
export kcore
export projectedGraph
export mutualInformation

# from HSBM.jl
export sampleSBM
export hypergraph

# from omega.jl
export partitionsUpTo
export partitionIntensityFunction
export IntensityFunction
export allOrNothingIntensityFunction
export sumOfExteriorDegreesIntensityFunction
export empiricalIntensityFunction

# from objectives.jl
export first_term_eval
export modularity
export logLikelihood

# from interence.jl
export estimateΩEmpirically
export learnParamters
export formObjective
export formObjectives

# from vol.jl
export evalSums
export evalConstants
export increments
export addIncrements
export second_term_eval
export momentIncrements

# from read_data.jl
export read_hypergraph_label_names
export read_hypergraph_labels
export read_hypergraph_edges
export read_hypergraph_data
export hyperedges

# from AON_hypergraph_louvain.jl; uses Random
export SuperNode_PPLouvain

# from general_hypegraph_louvain.jl
export HyperLouvain
export SuperNodeLouvain
export evalCuts

# from graph_louvain.jl
export CliqueExpansion
export CliqueExpansionModularity
export StarExpansionModularity
export computeDyadicResolutionParameter
export dyadicModularity
export dyadicLogLikelihood
export dyadicLogLikelihoodDirect

# from simple_synthetic.jl
export SimpleSyntheticHypergraph


end
