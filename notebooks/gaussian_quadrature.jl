# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.3
#   kernelspec:
#     display_name: Julia 1.3.1
#     language: julia
#     name: julia-1.3
# ---

# $$
# \int_a^b W(x) f(x) dx = \sum_{j=0}^{N-1} w_j f(x_j)
# $$
#
# $$
# g(x) = W(x) f(x)
# $$
#
# $$
# v_j = \frac{w_j}{W(x_j)}
# $$
#
# $$
# \int_a^b g(x) dx = \sum_{j=0}^{N-1} v_j g(x_j)
# $$

# +
using FastGaussQuadrature

a , b = - 2, 2
x, w = gausslegendre(4)
xr = 0.5 * (b-a)
xm = 0.5 * (b+a)
xm .+ xr .* x

# +
"""
Returns the integral of the function f between a and b, by ten-point Gauss- Legendre integration: 
the function is evaluated exactly ten times at interior points in the range of integration.
"""
function qgauss(f :: Function , a, b)
    @assert b > a
    n = 100
    x, w = gausslegendre(n)
    xr = 0.5 * (b-a)
    xm = 0.5 * (b+a)
    xr * sum( w[i] * f(xm + xr * x[i]) for i in 1:n )
    
end
# -

f(x) = exp( - x^2/2)
qgauss( f, -10, 10)

sqrt(2π)

using LinearAlgebra
x, w = gausslegendre(10)
dot( w,(x.^2))

# $$
# \int_{-\infty}^\infty e^{-x^2} f(x) dx = \sum^{N-1}_{j=0} w_jf(x_j)
# $$

gausshermite(10)
