using Polynomials

"""
    hermite_coeff( n )

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
function hermite_coefs( n )
    
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

Hermite(n::Int) = Poly(hermite_coefs(n))

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
    
    for n in eachindex(ht)
    
        h_n = Hermite(n) 
    
        result += 1 / sqrt(π * 2^n * factorial(n)) * ht[n] * h_n(x)
    
    end 

    return result

end

"""
    sht( data )

 performs a slow hermite transform in the most naive way possible directly from the
     data points given in xl()
"""
function sht(data, n)

    results = zeros(Float64, n)

    al(l) = sqrt(2.0/(l+1.0))
    bl(l) = 0.0
    cl(l) = (-1.0)*sqrt(l/(l+1.0))

    litn = length(data)
    bigc = sqrt(2n+1)
    xk(k) = ((k-litn)/litn)*bigc
    d0 = 1.0/pi^(1/4)

    for x in eachindex(data)
        lminus1 = 0
        cur_val = d0
        for y in eachindex(results)
            results[y] += data[x] * cur_val
            lminus2 = lminus1
            lminus1 = cur_val
            cur_val = (al(y-1)*(xk(x-1)) + bl(y-1))*lminus1 + cl(y-1)*lminus2;
        end
    end

    return results

end
