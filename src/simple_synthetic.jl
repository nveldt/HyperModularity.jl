function SimpleSyntheticHypergraph(n,m,K,pvals,rmin,rmax,cluster_sizes=ones(K),r_sizes=ones(rmax-rmin),cluster_prefs=ones(K))
    """
    Generate a synthetic hypergraph, store in incidence matrix
    Output several other things about hypergraph like edge sizes and degrees
    and edge list.

    This can be viewed as a simplified approximation to the HSBM, that is
    easier and faster to generate in practice.

    Procedure:
    * Set a fixed number of hyperedges m and clusters K (you can randomize this outside this function)
    * Set a distribution determining the relative proportion of hyperedges of each size from r_min to r_max
    * For each hyperedge, with probability p, assign it to one community at random,
        selecting a r nodes uniformly at random for a hyperedge of size r
    * Otherwise, with probability 1-p, select r nodes uniformly at random from across the entire hypergraph
    """

    # Normalize to get proportions/probabilities
    r_sizes = Weights(r_sizes/sum(r_sizes))
    cluster_sizes = Weights(cluster_sizes/sum(cluster_sizes))         # proportion of nodes in each cluster
    cluster_prefs = Weights(cluster_prefs/sum(cluster_prefs))        # relative proportion of interior hyperedge in each cluster

    # Generate the ground truth clustering of nodes
    # ground_truth = GenerateTruth(n,cluster_sizes)
    Clusters = Vector{Vector{Int64}}()
    for j = 1:K
        push!(Clusters,Vector{Int64}())
        @assert(n*cluster_sizes[j] > 2*rmax)
    end

    ground_truth = zeros(Int64,n)
    for i = 1:n
        c = sample(1:K,cluster_sizes)
        ground_truth[i] = c
        push!(Clusters[c],i)
    end

    U = Vector{Int64}()
    E = Vector{Int64}()
    EdgeList = Vector{Vector{Int64}}()
    E_lengths = zeros(Int64,m)
    for enum = 1:m
        e = SampleEdge(pvals,n,rmin,rmax,cluster_prefs,r_sizes,Clusters)
        push!(EdgeList,e)
        E_lengths[enum] = length(e)
        for node in e
            push!(U,node)
            push!(E,enum)
        end
    end
    He2n = SparseArrays.sparse(E,U,ones(length(U)),m,n)
    deg = vec(sum(He2n,dims = 1))
    return He2n, EdgeList, E_lengths, deg, ground_truth

end

function SampleEdge(pvals,n,rmin,rmax,cluster_prefs,r_sizes,Clusters)
    replace = false
    K = length(Clusters)
    r = sample(rmin:rmax,r_sizes)       # select a hyperedge size
    if rand(1)[1] < pvals[r-rmin+1]

        c = sample(1:K,cluster_prefs)   # select a random cluster to put the edge in
        return sort(sample(Clusters[c],r,replace = replace))
    else
        return sort(sample(1:n,r,replace = replace))
    end
end
