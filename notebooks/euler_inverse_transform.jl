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

# Points equally spaced
# $$
# x_k - x_{k-1} = \frac{C}{M_1}
# $$
#
# $$
# \hat{f}_{c,M_1}(n) = \hat{f}_c(n) = \sum_{k=0}^{M_1-1} f(x_k) c_n(x_k) \qquad 0 \leq n < M,
# $$
#
# $$
# c_n(x_k) = \frac{c_k C}{M} \sqrt{w(x_k)} p_n(x_k)
# $$
# $c_k$ depends on the integration technique
#
#


