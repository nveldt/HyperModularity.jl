# Testing partitionize
for n = 5:20
    for r = 5:10
        p = rand(1:n,r)
        d = countmap(p)
        v = sort(collect(values(d)),rev = true)
        @test HyperModularity.partitionize(p) == v
    end
end

v = vec([3 4 5 0 0])
@test HyperModularity.sortedremovezeros(v) == [3;4;5]
