# -*- coding: utf-8 -*-
# +
using Polynomials

"""
    Hermite(n)

Compute the nth "physicists' Hermite polynomials" 
using recurrence relation
"""
function Hermite(n::Int) 
    
    if( n == -1 ) return Poly([0.0]) end
    if( n == 0 ) return Poly([1.0]) end
    
    return Poly([0.0,2.0]) * Hermite(n-1) - 2*(n-1) * Hermite(n-2)
end

# +
"""
    hermite_coeffs( n )

For the physicists polynomials, assuming
```math
H_n(x) = \\sum^n_{k=0} a_{n,k} x^k,
```
we have:
```math
H_{n+1}(x) = 2x H_n(x) - 2n H_{n-1}(x).
```

Individual coefficients are related by the following recursion formula:

```math 
a_{n+1,k} = \\begin{cases}
 - a_{n,k+1} & k = 0, \\\\ 
 2 a_{n,k-1} - (k+1)a_{n,k+1} & k > 0,
\\end{cases}
```

and ``a_{0,0} = 1, a_{1,0} = 0, a_{1,1} = 2``.
"""
function hermite_coeffs( n )
    
    a = zeros(Float64, n+1, n+2)
    a[1, 1] = 1.0
    if n == 0 return [1.0] end
    a[2, 2] = 2.0

    for i in 2:n
        a[i+1, 1] = - a[i, 2]
        for k in 2:n+1
            a[i+1, k] = 2 * a[i, k-1] - k * a[i,k+1]
        end
    end
    
    a[n+1, :]
    
end

for i in 0:11
    println(i, ":", Poly(hermite_coeffs(i))-Hermite(i))
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
# -

# The Hermite transform expresses a function $f(x)$ in the complete orthonormal 
# basis of Hermite functions $\{\psi_n(x)\}^\infty_{n=0}=0$ where
#
# $$
# \psi_n(x)=(h_n)^{\frac{1}{2}} e^{-\frac{x^2}{2}} H_n(x)
# $$
#
# and $h_n = 2^n n! \sqrt{\pi}$ provided $H_n(x)$ is the nth Hermite polynomial. 
#
# The Hermite polynomials are determined by the three-term recurrence relation
#
# $$
# H_{n+1} = 2xH_n(x) - 2n H_{n-1}(x)
# $$
#
# with $H_{-1}(x) = 0$ and $H_0(x) = 1$

typeof(Hermite(1))

# +
struct Psi
    
    n :: Int
    H :: Poly{Float64}
    
    function Psi(n)
        
        new(Hermite(n))
        
    end
    
end

function (psi::Psi)(x :: Float64)
    
    h = 2^n * factorial(n) * sqrt(π)
    return 
    
end
    
    
# -

# $f$ can be represented as
#
# $$
# f(x) = \sum_{n=0}^\infty \hat{f}(n) \psi_n(x)
# $$
#
# where $\hat{f}(n)$ is the nth Hermite coefficient and is defined by
#
# $$
# \hat{f}(n) = <f,\psi_n> = \int_{-\infty}^{\infty} f(x) \psi_n(x) dx
# $$

# +
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




