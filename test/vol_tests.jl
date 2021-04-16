n = 20
Z = rand(1:5, n)
ϑ = dropdims(ones(1,n) + rand(1,n), dims = 1)
μ = mean(ϑ)
kmax = 4

ω(x, α) = (x[1]+1)^(-α[x[2]]*x[2])

α = repeat([2.0], kmax)

Ω = sumOfExteriorDegreesIntensityFunction(ω, kmax)
H = sampleSBM(Z, ϑ, Ω;α=α, kmax=kmax, kmin = 1)
D = H.D
all_partitions = collect(Iterators.flatten([(partitions(i,j)) for i =1:kmax for j=1:i]))

p = [3, 1, 1]

s1 = evalSumNaive(p, Z, D)
s2 = evalSumNaive2(p, Z, D)

## d
@test s1 == s2

sumsNaive = evalSumsNaive(Z,D,kmax)
sumsPV    = evalSumsPV(Z,D,kmax)
sums      = evalSums(Z,D,kmax)[3]

@test all([sumsNaive[p] == sumsPV[p] for p in all_partitions])
@test all([sumsNaive[p] == sums[p] for p in all_partitions])

# voldiff

ℓ = maximum(Z) # number of total groups

# let's move all the nodes in group 4 to group 5
to_move = findall(==(4), Z)
t = 5
Z_ = copy(Z)

## first we'll step through, moving each node one-by-one
V, μ, M = evalSums(Z, D, kmax; constants=false, bigInt=false)

for i in to_move
    global V, μ, M
    ΔV, Δμ, ΔM = increments(V, μ, M, i, t, D, Z_)
    V, μ, M = HyperModularity.addIncrements(V, μ, M, ΔV, Δμ, ΔM)
    # carry out the change in membership
    Z_[i] = t
end
step_wise = M

# next we'll compute the complete set of increments at once
V, μ, M = evalSums(Z, D, kmax; constants=false, bigInt=false);
ΔV, Δμ, ΔM = increments(V, μ, M, to_move, t, D, Z)
V, μ, M = HyperModularity.addIncrements(V, μ, M, ΔV, Δμ, ΔM)
batch = M
#
@test all([step_wise[p] == batch[p] for p in keys(M)])
