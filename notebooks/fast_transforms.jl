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

# +
using FastTransforms, FastGaussQuadrature, Plots

hermitepoints(n) = FastGaussQuadrature.unweightedgausshermite( n )[1]

weightedhermitetransform(exp.(-hermitepoints(2).^2/2))
# -

weightedhermitetransform(exp.(-hermitepoints(3).^2/2))

# $$
#     f(x) = e^{-\frac{x^2}{2}} \frac{ 2x}{\sqrt(2)}
# $$

x = hermitepoints(100)
f = exp.(-x.^2 ./ 2) .* 2x/sqrt(2)
f̂ = iweightedhermitetransform([0.0; 1.0; zeros(98)])
plot( [ x x], [f f̂])
savefig("in.svg")

# $$
#     f(x) = e^{-\frac{x^2}{2}} \frac{4x^2 - 2}{2\sqrt{2}}
# $$

f = exp.(-x.^2 ./ 2) .* (4x.^2 .- 2)/(sqrt(2)*2)
f̂ = iweightedhermitetransform([0.0; 0; 1.0; zeros(97)])
plot( [ x x], [f f̂])

# $$
#     f(x) = e^{-\frac{x^2}{2}} \frac{-12x + 8x^3}{4\sqrt{3}}
# $$

f = exp.(-x.^2 ./ 2) .* (-12x + 8x.^3) / (sqrt(2*3)*2^(3/2))
f̂ = iweightedhermitetransform([0.0; 0; 0; 1.0; zeros(96)])
plot( [ x x], [f f̂])

# $$
#     f(x) = e^{-\frac{x^2}{2}} \frac{12 - 48x^2 + 16x^4}{8\sqrt{6}}
# $$

# +
f = exp.(-x.^2 ./ 2) .* (12 .- 48x.^2 .+ 16x.^4) / (sqrt(2*3*4)*2^(4/2))
f̂ = iweightedhermitetransform([0.0; 0; 0; 0; 1.0; zeros(95)])
plot( [ x x], [f f̂])



# +
import FastTransforms: ForwardWeightedHermitePlan
import FastTransforms: BackwardWeightedHermitePlan

n = 100
fht = ForwardWeightedHermitePlan(n)
iht = BackwardWeightedHermitePlan(n)

# +
using BenchmarkTools
f = rand(n)

@btime iht * ( fht * f );

# +
import FastTransforms: weightedhermitetransform
import FastTransforms: iweightedhermitetransform

@btime iweightedhermitetransform( weightedhermitetransform(f) );
# -


