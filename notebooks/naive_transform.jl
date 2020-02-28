# -*- coding: utf-8 -*-
using FastGaussQuadrature, FastTransforms
using LinearAlgebra
using Polynomials
using Plots

# The Hermite transform expresses a function $f(x)$ in the complete orthonormal 
# basis of Hermite functions $\{\psi_n(x)\}^\infty_{n=0}=0$ where
#
# $$
# \psi_n(x)=(h_n)^{-\frac{1}{2}} e^{-\frac{x^2}{2}} H_n(x)
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
#
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

"""
    Hermite(n)

Compute the nth "physicists' Hermite polynomials" 
using recurrence relation
"""
function Hermite(n::Int) 
    
    if n == -1  return Poly([0.0]) end
    if n == 0   return Poly([1.0]) end
    
    H_n = Poly([0.0,2.0]) * Hermite(n-1) - 2*(n-1) * Hermite(n-2)
    
    return H_n
end

# # Efficient transforms
#
# Orthogonal polynomial bases:
# $$
# \psi_n(x) = \sqrt{w(x)} p_n(x)
# $$
#
# with $w(x)$ the *weight function* and the $\{p_n\}^\infty_{n=0}$ are a collection of orthogonal polynomials:
#
# $$
# p_{n+1} = (a_nx + b_n) p_n(x) - c_n p_{n-1}(x)
# $$
# with $p_{-1} = 0, p_0(x) = d_0 = \frac{1}{\pi^{1/4}} $
#
# The absorbed form is 
# $$
# p_{n+1} = x \sqrt{\frac{2}{n+1}} p_n(x) - \sqrt{\frac{n}{n+1}} p_{n-1}(x)
# $$
#

# +
const BIGN = 2
const LITN = BIGN*5
const BIGC = sqrt(2*BIGN+1)
const D0 = 1.0/π^(1/4)

xk(k :: Int) = ((k-LITN)/LITN)*BIGC

# +
using QuadGK
n = 1000
y(x) = exp.(-x^2/2)

x, w = QuadGK.gauss(n)
x .*= 5
sum( w )
# -

quadgk( y, -10, 10)

n = 10
x, w = FastGaussQuadrature.gausshermite(n)
sum(w)

x, w = FastGaussQuadrature.unweightedgausshermite(n)
sum(w)

# $$
# xk(n) = \frac{k - n}{n}  \sqrt{2N+1}
# $$



plot(x, w)

# +
n, N = LITN, BIGN

function generate_data( n )
    data = ones(2n+1)
    x = [xk(i) for i in 0:2n]
    data .*= exp.(-x.^2/2)
    x, data
end

x, data = generate_data(n)
plot(x, data; m=:o, mc=:red, ms=2, label=:data)


# +
function naive_transform( data, N, n )
    
    a(n) = sqrt(2.0/(n+1.0))
    b(n) = 0.0
    c(n) = (-1.0)*sqrt(n/(n+1.0))
    
    Lminus1 :: ComplexF64 = 0.0
    Lminus2 :: ComplexF64 = 0.0
    curVal  :: ComplexF64 = 0.0
    
    results = zeros(Float64, N)
    
    for x in eachindex(data)
        lminus1 = 0.0
        curVal = D0
        for y in 1:N
            results[y] += data[x] * curVal
            Lminus2 = Lminus1
            Lminus1 = curVal
            curVal = (a(y-1)*(xk(x-1)) + b(y-1))*Lminus1 + c(y-1)*Lminus2;
        end
    
    end
    BIGC = sqrt(2N+1)
    results .* (2 * BIGC/n)
end



# -

using FastGaussQuadrature, FastTransforms
hermitepoints(n) = FastGaussQuadrature.unweightedgausshermite( n )[1]
weightedhermitetransform(exp.(-hermitepoints(2).^2/2))

# ## Compute weights
#
# ### Gauss quadrature transform
#
# we compute 
#
# $$
# \hat{f}_g(n) = \sum^{M-1}_{k=0} f(x_k) w_k p(x_k)}{w(x_k)}
# $$
#
#  
# $x_k$ are the abscissae of the root $p_M(x)$ where $M$ is the bandlimit.
#

# +
struct Psi
    n :: Int
    h :: Float64
    H :: Poly{Float64}
    function Psi(n :: Int)
        H = Hermite(n)
        h = 1. / sqrt(2^n * factorial(n) * sqrt(π))
        new( n, h, H)
    end
end

function (psi :: Psi)(x :: Float64)
    psi.h * psi.H(x) * exp(-x^2/2)
end




# +
f̂ = naive_transform( data, N, n )

y = zero(x)
for (i,px) in enumerate(x)
    for n in 1:N
        ψ = Psi(n-1)
        y[i] += f̂[n] * ψ(px)
    end
end 

plot([x x], [data y])
# -


