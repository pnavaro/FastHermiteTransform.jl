using Test
using ImportMacros
import Base.Libc:rand, srand

@using FastHermiteTransform as FHT


n = FHT.LITN
N = FHT.BIGN

RAND_MAX = 2^31-1

data = zeros(Float64, 2n+1)

open("data.txt") do file
    for i in eachindex(data)
        data[i] = parse(Float64, readline(file))
    end
end

@test all(data .== [0.852779, 0.000549, 0.136753, 0.529954, 0.823420, 0.388714, 0.261452, 0.719810, 0.266854, 0.903989, 0.524461])


fancyResult = zeros(N)
naiveResult = zeros(N)

srand(0)
@show daConstant = -1^(rand()%2)*rand()/RAND_MAX

diag = zeros(Float64, (2n+1,2n+1))

for i in 1:2n+1
    diag[i,i] = exp(-FHT.xk(i)^2/2)
    data[i] *= diag[i,i]
end

@show data

#FHT.initRns(N)
#
#FHT.oneDTransform(data,fancyResult);
#
#FHT.naiveTransform(data,naiveResult);
