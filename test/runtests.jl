using Test
using ImportMacros
import Base.Libc.rand

@using FastHermiteTransform as FHT

@test true

n = FHT.LITN
N = FHT.BIGN

RAND_MAX = 2^31-1

data = zeros(Float64, 2n+1)

fancyResult = zeros(N)
naiveResult = zeros(N)

daConstant = -1^(rand()%2)*rand()/RAND_MAX

diag = zeros(Float64, (2n+1,2n+1))

for i in 1:2n+1
    diag[i,i] = exp(-FHT.xk(i)^2/2)
    data[i] *= diag[i,i]
end

FHT.initRns(N)

FHT.oneDTransform(data,fancyResult);

FHT.naiveTransform(data,naiveResult);
