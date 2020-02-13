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
