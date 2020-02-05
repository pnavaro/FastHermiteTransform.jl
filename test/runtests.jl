using Test
using FastHermiteTransform

@test true

int n = LITN;
int N = BIGN;
int i,j,x;
tester = fopen("doc/testing.txt", "wt+");

data = zeros(n*2+1)

fancyResult = zeros(N)
naiveResult = zeros(N)

daConstant = pow(-1,rand()%2)*rand()/RAND_MAX;

diag = zeros((2*n+1)*(2*n+1));

for i in eachindex(diag)
    diag[(2*n+1)*i+i] = exp(0-pow(xk(i),2)/2);
    data[i] *= diag[(2*n+1)*i+i];
end

initFastFouriers(N);
initRns(N);

oneDTransform(data,fancyResult);

naiveTransform(data,naiveResult);
