c = rand(1:10,10)
cnew = HyperModularity.renumber(c)
@test maximum(cnew) == length(unique(cnew))
@test sort(unique(cnew)) == collect(1:maximum(cnew))
for j = 1:maximum(cnew)
    S = findall(x->x==j,cnew)
    orig = c[S]
    t = orig[1]
    @test all(orig .== t)
end

n = 10
m = rand(5:15) # number edges
Edges = Vector{Vector{Int64}}()
for i = 1:m
    k = rand(2:5)
    E = sample(1:n,k,replace = false)
    push!(Edges,sort(E))
end
He2n = elist2incidence(Edges,n)
Hn2e = sparse(He2n')
node2edge = incidence2elist(Hn2e)
@test size(He2n,1) == m
@test m == length(Edges)
@test size(He2n,2) == n
@test n == length(node2edge)

for i = 1:m
    @test findall(x->x==1,He2n[i,:]) == Edges[i]
    @test Edges[i] == HyperModularity.getnodes(Hn2e,i)
end

N = NeighborList(node2edge,Edges)
N2 = NeighborList(He2n,Hn2e)
@test N == N2
