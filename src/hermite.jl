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
    
        result += 1 / sqrt(π * 2^n * factorial(n)) * ht(n) * h_n(x)
    
    end 

    return result


end
