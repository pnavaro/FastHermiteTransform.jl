# -*- coding: utf-8 -*-
# +
using Polynomials

function coefs( n::Int )

    if n == 0 return [0.0] end
    if n == 1 return [1.0] end
    
    m = zeros(Float64, (n, n))
    m[1, 1] = 1.0
    m[2, 2] = 1.0
      
    for l in 3:n
        m[1, l] = - (l - 2) * m[1, l-2]
        for i in 2:n
            m[i, l] = m[i-1, l-1] - (l - 2) * m[i, l-2]
        end
    end

    m[:,n]

end

Hermite(n::Int) = Poly(coefs(n))

"""
    isht( ht, x)

Inverse Hermite transform

```math
H^{-1}\\{f_H(n)\\} = F(x) = \\sum_{n=0}^\\infty \\frac{1}{\\sqrt\\pi 2^n n!} f_H(n) H_n(x)
```

```math
H_{n+1}(x) = x H_{n}(x) - n H_{n-1}(x)
```
"""
function isht( ht, x )

    result = 0.0
    
    N = length(ht)
    
    for n in 0:N-1
    
        h_n = Hermite(n) 
    
        result += 1 / sqrt(π * 2^n * factorial(n)) * ht[n+1] * h_n(x)
    
    end 

    return real(result)


end
# -

using Plots

# +
BIGN = 8
LITN = BIGN*5
BIGC = sqrt(2*BIGN+1)
Al(l) = sqrt(2.0/(l+1.0))
Bl(l) = 0.0
Cl(l) = (-1.0)*sqrt(l/(l+1.0))
D0 = 1.0/pi^(1/4)

xk(k :: Int) = ((k-LITN)/LITN)*BIGC



# +
n = LITN
N = BIGN
data = zeros(2n+1)
for i in eachindex(data)
    if i > n÷2 && i <= 3n÷2+1
        data[i] = 1.0
    else
        data[i] = 0.0
    end
end

diag = zeros((2n+1,2n+1));

for i in 1:2n+1
    diag[i,i] = exp(-xk(i-1)^2/2)
    data[i] *= diag[i, i]
end

# +
PI = 3.14159265
cur_val = 1.0/(PI^(0.25))
results = zeros(Float64, N)
for x in 1:2n+1
    lminus1 = 0.0
    
    cur_val = 1.0/(PI^(0.25))
    for y in 1:N
        results[y] += data[x] * cur_val
        lminus2 = lminus1
        lminus1 = cur_val
        cur_val = (Al(y-1)*(xk(x-1)) + Bl(y-1))*lminus1 + Cl(y-1)*lminus2;
    end
    
end

ht = results .* (2 * BIGC/n)

# +
x = LinRange(0,10,100)
y = zeros(100)
for (i,px) in enumerate(x)
    for n in 1:N
        h_n = Hermite(n-1)
        hn = 2^(n-1) * factorial(n) * sqrt(π)
        y[i] += hn^(-1/2) * ht[n] * h_n(px)
    end
end 

plot(x, y)
# -


