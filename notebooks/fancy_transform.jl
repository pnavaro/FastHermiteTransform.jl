# -*- coding: utf-8 -*-
using FFTW

const BIGN = 8
const LITN = (BIGN * 5)
const BIGC = sqrt(2 * BIGN + 1)

const ALPHA = (2.0 / BIGC)
const BETA = 0
const GAMMA = (-1.0)
const D0 = 1.0 / pi^(1.0 / 4.0)


AL(l) = sqrt(2.0 / (l + 1.0))
BL(l) = 0
CL(l) = ((-1.0) * sqrt(l / (l + 1.0)))

UL(l) = (AL(l) / ALPHA)
VL(l) = (BL(l) + (-1) * (UL(l) * BETA))
WL(l) = ((-1) * GAMMA * UL(l))

xk(k) = ((k - LITN) / LITN) * BIGC;

daRns = zeros(Float64, (BIGN + 1) * BIGN * 8 * BIGN);
daRnsSize = BIGN;

# +
"""
give me the first column of a circulant matrix in M.
"""
function circulantVcMatrixMultiply(c, v, result)

    result .= ifft(fft(c) .* fft(v))

end

# +
"""
 A B  *  E F  = AE+BG AF+BH
 C D     G H    CE+DG CF+DH
"""
function fourBcirculantSqMatrixMultiply(M1, M2, n, result)

    @debug "doing a multiplication of size $n"

    temp1 = zeros(Float64, n ÷ 2)
    temp2 = zeros(Float64, n ÷ 2)

    # fill up the columns
    A = @view M1[:]
    E = @view M2[:]
    B = @view M1[n÷2:end]
    F = @view M2[n÷2:end]
    C = @view M1[n:end]
    G = @view M2[n:end]
    D = @view M1[3*n÷2:end]
    H = @view M2[3*n÷2:end]

    # A*E+B*G top left
    circulantVcMatrixMultiply(A, E, temp1)
    circulantVcMatrixMultiply(B, G, temp2)

    # Add em up
    for i in eachindex(temp1)
        result[i] = temp1[i] + temp2[i]
    end

    # A*F+B*H top right
    circulantVcMatrixMultiply(A, F,  temp1)
    circulantVcMatrixMultiply(B, H,  temp2)
    # Add em up
    for i in eachindex(temp1)
        result[i+n÷2] = temp1[i] + temp2[i]
    end

    # C*E+D*G bottom left
    circulantVcMatrixMultiply(C, E, temp1)
    circulantVcMatrixMultiply(D, G, temp2)

    # Add em up
    for i in eachindex(temp1)
        result[i+n] = temp1[i] + temp2[i]
    end

    # C*F+D*H bottom right
    circulantVcMatrixMultiply(C, F,  temp1)
    circulantVcMatrixMultiply(D, H,  temp2)
    # Add em up
    for i in eachindex(temp1)
        result[i+3*n÷2] = temp1[i] + temp2[i]
    end


end
# -


"""
define An(l) as defined in the paper
result is a 8*n sized vector with each 2n representing a circulant block
"""
function createAn(n, l, result)
    # top left is all zeros (we'll zero out everything else while we're at it.
    fill!(result, 0)
    # top right is I2n
    result[2*n+1] = 1
    # bottom left is cl*I2n
    result[4*n+1] = CL(l)
    # bottom right is Cn(wl,vl,ul);
    result[6*n+1] = VL(l)
    result[6*n+2] = WL(l)
    result[8*n] = UL(l)
end

# +
"""
compute desired Rn 
needed columns of circulant matricies listed vertically as 
top left, top right, bottom left, bottom right.
"""
function precomputeRn(n, l, result)

    # now it's not so we open our write file and start computing

    temp = zeros(Float64, 8 * n)

    temp2 = zeros(Float64, 8 * n)

    createAn(n, n / 2 + l, result)

    for i = l+n÷2-1:-1:l
        createAn(n, i, temp)
        fourBcirculantSqMatrixMultiply(result, temp, 4 * n, temp2)
        result .= temp2
    end


end
# -

"""
multiply Z by the Rn that was precomputed at n,l
"""
function preFourBcirculantVcMatrixMultiply(n, l, Vec, result)
    n *= 4
    temp = zeros(Float64, n ÷ 2)
    temp2 = zeros(Float64, n ÷ 2)

     # top left (want a column of the top left times first half of Vec)
    circulantVcMatrixMultiply(
        daRns[daRnsSize * 8 * (daRnsSize * n / 4 + l)],
        Vec,
        n ÷ 2,
        result,
    )

    # top right (want a column of the top right times second half of Vec)
    circulantVcMatrixMultiply(
        daRns[daRnsSize*8*(daRnsSize * n / 4 + l)+n/2],
        Vec[n÷2+1:end],
        n ÷ 2,
        temp,
    )

    # add top left and top right
    for i = 1:n÷2
        result[i] += temp[i]
    end

    # bottom left (want a column of the bottom left times first half of Vec)
    circulantVcMatrixMultiply(
        daRns[daRnsSize*8*(daRnsSize * n ÷ 4 + l)+n],
        Vec,
        n ÷ 2,
        result + n ÷ 2,
    )

    # bottom right (want a column of the bottom right times second half of Vec)
    circulantVcMatrixMultiply(
        daRns[daRnsSize*8*(daRnsSize * n ÷ 4 + l)+3*n÷2],
        Vec[n÷2+1:end],
        n ÷ 2,
        temp2,
    )

    # add bottom left and bottom right
    for i = n÷2+1:n
        result[i] += temp2[i-n÷2]
    end

    n ÷= 4
end

# +
"""
a recursive function which performs a Hermite transform in O(n(logn)^2) time
Z0 and Z1 must be precomputed as defined in the paper.  l should be first
  set to 1
you must precompute all the necessary Rns.
"""
function performTransform(Z0, Z1, n, l, result)
    result[l-1] = Z0[n-1]
    result[l] = Z1[n-1]

    if (n < 3)
        return
    end

     # temp to store the new data
    temp = zeros(Float64, 4 * n)


     # combine Z0 and Z1 into Z to get ready for the matrix multiply
    Z = vcat(Z0, Z1)
    preFourBcirculantVcMatrixMultiply(n, l - 1, Z, temp)

    nover2 = n ÷ 2
    performTransform(Z0[nover2+1:end], Z1[nover2+1:end], nover2, l, result)
    performTransform(temp[nover2+1:end], temp[5*nover2+1:end], nover2, l + nover2, result)

end


# +
function main()

    n = LITN
    N = BIGN

    data = zeros(Float64, (n * 2 + 1))

    result = zeros(Float64, N)

    for i in eachindex(data)
        data[i] = exp(-xk(i)^2 / 2)
    end


    # precompute the necessary Rns
    i = N
    while i >= 4
        j = 1
        while j <= N
            precomputeRn(i, j, daRns[8 * N * (N * i + j)])
            j += i
        end
        i /= 2
    end


    Z0 = zeros(ComplexF64, 2 * N)
    dblZ0 = zeros(ComplexF64, 2 * N)
    Z1 = zeros(ComplexF64, 2 * N)
    dblZ1 = zeros(ComplexF64, 2 * N)

    # perform a chebyshev transform in the most naive way possible directly from the
    # data points defined by xl()

    for i in eachindex(data)
        Lminus1 = xk(x) / BIGC
        Lminus2 = 2.0 * Lminus1^2 - 1  # 2*x^2-1
        for j = 1:BIGN
            curVal = (ALPHA * (xk(i)) + BETA) * Lminus1 + GAMMA * Lminus2
            Z0[j+N-1] += data[i] * curVal
            Lminus2 = Lminus1
            Lminus1 = curVal
        end
    end


    Z0[2*N-1] = 0

    # we only want the real parts
    for i = 1:N
        dblZ0[N-1+i] = real(Z0[N-1+i])
    end

    # expand the data
    for i = 1:N
        dblZ0[i] = dblZ0[2*N-i-2]
    end


    # find the next data point
    # use Zl to calculate Zl+1.  Z0 = Zl.  Z1 = Zl+1.
    # n is the length  of Z0
    # normalize the chebyshev result

    for i = 1:2N
        dblZ0[i] *= D0
    end
        # now figure out Z1
    for i = 2:2N-1
        dblZ1[i] = AL(0) * (1 / ALPHA * dblZ0[i+1] - BETA / ALPHA * dblZ0[i] -
                    GAMMA / ALPHA * dblZ0[i-1]) + BL(0) * dblZ0[i]
    end


    # do the second part
    performTransform(dblZ0, dblZ1, N, 1, result)

    for i = 1:N
        result[i] *= (2 * BIGC) / n
        println(result[i])
    end

end
