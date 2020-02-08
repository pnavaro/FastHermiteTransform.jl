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
