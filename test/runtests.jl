using Test
using ImportMacros
import Base.Libc:rand, srand

@using FastHermiteTransform as FHT


n = FHT.LITN
N = FHT.BIGN

RAND_MAX = 2^31-1

data = zeros(Float64, 2n+1)

data[n÷2+1:3n÷2+1] .= 1.0

fancyResult = zeros(N)
naiveResult = zeros(N)

diag = zeros(Float64, (2n+1,2n+1))

for i in 1:2n+1
    diag[i,i] = exp(-FHT.xk(i-1)^2/2)
    data[i] *= diag[i,i]
end

@show data

rns = FHT.Rns(N)

#FHT.oneDTransform(data,fancyResult);
#
#FHT.naiveTransform(data,naiveResult);
