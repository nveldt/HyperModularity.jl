c = rand(1:20,10)
cnew = renumber(c)
@test maximum(cnew) == length(unique(cnew))
@test sort(unique(cnew)) == collect(1:maximum(cnew))
for j = 1:length(cnew)
    S = findall(x->x==j,cnew)
    @test ~notsame(c[S])
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
Hn2e = transpose(He2n)
node2edge = incidence2elist(Hn2e)
@test size(He2n,1) == m == length(Edges)
@test size(He2n,2) == n == length(node2edge)

for i = 1:m
    @test findall(x->x==1,He2n[i,:]) == Edges[i]
    @test Edges[i] == getnodes(He2n,i)
end

N = NeighborList(node2edge,Edges)
N2 = NeighborList(He2n,Hn2e)
@test N == N2
