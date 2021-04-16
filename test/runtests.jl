#cd(homedir()*"/GitHubRepos/Working-directory/HyperModularity")
# using Revise
# using Pkg; Pkg.activate(".")
using HyperModularity
using Test
using StatsBase
using Combinatorics

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
    ϑ = dropdims(ones(1,n) + rand(1,n), dims = 1)
    μ = mean(ϑ)
    kmax = 4

    ω(x, α) = (x[1]+1)^(-α[x[2]]*x[2])

    α = repeat([2.0], kmax)

    Ω = sumOfExteriorDegreesIntensityFunction(ω, kmax)
    H = sampleSBM(Z, ϑ, Ω;α=α, kmax=kmax, kmin = 1)

    @test length(keys(H.E)) == kmax
end

@testset "read_data" begin
    dataset = "contact-high-school-classes"
    maxsize = 5
    minsize = 2
    return_labels = true
    H, Z = read_hypergraph_data(dataset,maxsize,minsize,return_labels)
    @test typeof(H) == hypergraph
    @test length(Z) == length(H.D)
end

@testset "vol" begin
    # some tests here are erroring
    # include("vol_tests.jl")
end

@testset "louvain_utils" begin

    # Nate will eventually write tests for all of these

end

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
