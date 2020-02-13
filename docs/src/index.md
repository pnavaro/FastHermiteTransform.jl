# FastHermiteTransform.jl

The Hermite transform expresses a function `f(x)` in the complete orthonormal `L2 (|R)`-basis of Hermite functions `{\psi_n(\xi)}^\infty_{n=0}=0` where

```math
\psi(\xi)=(h_n)^{\frac{1}{2}}− e^{-x^2/2} H_n(\xi)
```
and `h_n = 2^n n! \sqrt{\pi}` provided `H_n(x)` is the nth Hermite polynomial. The Hermite polynomials are determined by the three-term recurrence relation

$$
H_{n+1} = 2xH_n(x) - 2n H_{n-1}(x)
$$

$H_{-1}(x) = 0$ and $H_0(x) = 1$

f can be represented as
$$
f(x) = \sum_{n=0}^\infty \hat{f}(n) \psi_n(x)
$$

where $\hat{f}(n)$ is the nth Hermite coefficient and is defined by

$$
\hat{f}(n) = <f,\psi_n> = \int_{-\infty}^{\infty} f(x) \psi_n(x) dx
$$

The Hermite transform of a function ``F(x)`` is

```math
H\{F(x)\} = f_H(n) = \int_{-\infty}^\infty e^{-x^2} \ H_n(x)\  F(x) \ dx
```

The inverse Hermite transform is given by

```math
H^{-1}\{f_H(n)\} = F(x) = \sum_{n=0}^\infty \frac{1}{\sqrt\pi 2^n n!} f_H(n) H_n(x)
```
