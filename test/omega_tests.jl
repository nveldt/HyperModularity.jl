
kmax = 10

# general partitionize-based intensity function
Ω = partitionIntensityFunction(p->sum(p), kmax)
@test typeof(Ω) == IntensityFunction

# all-or-nothing cut
Ω = allOrNothingIntensityFunction(x->(x[1]+1)^(-x[2]), kmax)
@test typeof(Ω) == IntensityFunction

# number of unique elements in z: generalizes all-or-nothing
Ω = sumOfExteriorDegreesIntensityFunction((x, α)->(x[1]+1)^(-α[x[2]]*x[2]), kmax)
@test typeof(Ω) == IntensityFunction
