function read_hypergraph_label_names(dataname::String)
    names = String[]
    pathname = joinpath(dirname(dirname(@__FILE__)),"data")
    # open("data/$dataname/label-names-$dataname.txt") do f
    open("$pathname/$dataname/label-names-$dataname.txt") do f
        for line in eachline(f)
            push!(names, line)
        end
    end
    return names
end

function read_hypergraph_labels(dataname::String)
    labels = Int64[]
    pathname = joinpath(dirname(dirname(@__FILE__)),"data")
    # open("data/$dataname/node-labels-$dataname.txt") do f
    open("$pathname/$dataname/node-labels-$dataname.txt") do f
        for line in eachline(f)
            push!(labels, parse(Int64, line))
        end
    end
    return labels
end

function read_hypergraph_edges(dataname::String, maxsize::Int64=25, minsize::Int64=2)
    E = Dict{Integer, Dict}()
    pathname = joinpath(dirname(dirname(@__FILE__)),"data")
    # open("data/$dataname/hyperedges-$dataname.txt") do f
    open("$pathname/$dataname/hyperedges-$dataname.txt") do f
        for line in eachline(f)
            edge = [parse(Int64, v) for v in split(line, ',')]
            sort!(edge)
            if minsize <= length(edge) <= maxsize
                sz = length(edge)
                if !haskey(E, sz)
                    E[sz] = Dict{}()
                end
                E[sz][edge] = 1
            end
        end
    end
    return E
end

function read_hypergraph_data(dataname::String, maxsize::Int64=25, minsize::Int64=2, return_labels=true)
    E = read_hypergraph_edges(dataname, maxsize, minsize)

    n = maximum([maximum(e) for k in keys(E) for e in keys(E[k])])
    D = zeros(Int64, n)
    for (sz, edges) in E
        for (e, _) in edges
            D[e] .+= 1
        end
    end

    maxedges = maximum(keys(E))
    for k in 1:maxedges
        if !haskey(E, k)
            E[k] = Dict{}()
        end
    end

    N = 1:n

    if return_labels
        labels = read_hypergraph_labels(dataname)
        return hypergraph(N, E, D), labels[N]
    end

    return hypergraph(N, E, D)
end

function hypermodularity_datasets()
    println("HyperModularity Package Datasets \n")
    pathname = joinpath(dirname(dirname(@__FILE__)),"data")
    foreach(readdir(pathname)) do f
       if isdir("$pathname/"*f)
           println("\t $f")
       end
    end
end
