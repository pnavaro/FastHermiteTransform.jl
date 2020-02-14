# FastHermiteTransform.jl

The Hermite transform expresses a function ``f(x)`` in the complete orthonormal 
basis of Hermite functions ``\{\psi_n(x)\}^\infty_{n=0}=0`` where

```math
\psi_n(x)=(h_n)^{\tfrac{1}{2}} − \euler^{-x^2/2} H_n(x)
```

and ``h_n = 2^n n! \sqrt{\pi}`` provided ``H_n(x)`` is the nth Hermite polynomial. 
The Hermite polynomials are determined by the three-term recurrence relation

```math
H_{n+1} = 2xH_n(x) - 2n H_{n-1}(x)
```

with ``H_{-1}(x) = 0`` and ``H_0(x) = 1``

``f`` can be represented as

```math
f(x) = \sum_{n=0}^\infty \hat{f}(n) \psi_n(x)
```

where ``\hat{f}(n)`` is the nth Hermite coefficient and is defined by

```math
\hat{f}(n) = <f,\psi_n> = \int_{-\infty}^{\infty} f(x) \psi_n(x) dx
```
