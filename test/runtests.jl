#cd(homedir()*"/GitHubRepos/Working-directory/HyperModularity")
# using Revise
# using Pkg; Pkg.activate(".")
using HyperModularity
using Test
using StatsBase
using Combinatorics
using SparseArrays

include("../src/test_funs.jl")

@testset "utils" begin
    include("utils_tests.jl")
end

@testset "Ω" begin
    include("omega_tests.jl")
end

@testset "hsbm" begin
    n = 20
    Z = rand(1:5, n)
    Z = Int64.(Z)
    ϑ = dropdims(ones(1,n) + rand(1,n), dims = 1)
    μ = mean(ϑ)
    kmax = Int64(4)

    ω(x, α) = (x[1]+1)^(-α[x[2]]*x[2])

    α = repeat([2.0], kmax)

    Ω = sumOfExteriorDegreesIntensityFunction(ω, kmax)
    H = sampleSBM(Z, ϑ, Ω;α=α, kmax=kmax, kmin = 1)

    @test length(keys(H.E)) == kmax
end

@testset "louvain_utils" begin

    include("louvain_utils_tests.jl")
end

# @testset "vol" begin
#     # some tests here are erroring
#     # include("vol_tests.jl")
# end

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
