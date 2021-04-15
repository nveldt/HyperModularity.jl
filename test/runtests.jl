#cd(homedir()*"/GitHubRepos/Working-directory/HyperModularity")
# using Revise
# using Pkg; Pkg.activate(".")
using HyperModularity
using Test
using StatsBase

@testset "utils" begin
    include("utils_tests.jl")
end


# @testset "Ω" begin
#     include("omega_tests.jl")
# end
#

#
# @testset "vol" begin
#     # create a .jl file for testing volume functions
#     # include("vol_tests.jl")
# end
#
# @testset "read_data" begin
#     # file for testing read data, if desired
# end
#
# @testset "hsbm" begin
#
# end
#
# @testset "inference" begin
#
# end
#
# @testset "analysis_helpers" begin
#
# end
#
# @testset "louvain_utils" begin
#
# end
#
# @testset "graph_louvain" begin
#
# end
#
# @testset "general_louvain" begin
#
# end
#
# @testset "aon_louvain" begin
#
# end
