# FastHermiteTransform.jl

The Hermite transform of a function ``F(x)`` is

```math
H\{F(x)\} = f_H(n) = \int_{-\infty}^\infty e^{-x^2} \ H_n(x)\  F(x) \ dx
```

The inverse Hermite transform is given by

```math
H^{-1}\{f_H(n)\} = F(x) = \sum_{n=0}^\infty \frac{1}{\sqrt\pi 2^n n!} f_H(n) H_n(x)
```
