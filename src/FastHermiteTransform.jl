module FastHermiteTransform

    using FFTW
    include("hermite.jl")
    
    const BIGN = 8
    const LITN = BIGN*5
    const BIGC = sqrt(2.0*BIGN+1)
    
    const ALPHA = 2.0/BIGC
    const BETA  = 0
    const GAMMA = -1.0
    const D0 = 1.0/pi^(1/4)
    
    
    AL(l) = sqrt(2.0/(l+1.0))
    BL(l) = 0.0
    CL(l) = (-1.0)*sqrt(l/(l+1.0))
    
    UL(l) = AL(l)/ALPHA
    VL(l) = BL(l)+(-1)*(UL(l)*BETA)
    WL(l) = (-1)*GAMMA*UL(l)
    
    include("rns.jl")
    
    """
    defines the data points
    """
    xk(k) = ((k-LITN)/LITN)*BIGC
    
    
    """
    use Zl to calculate Zl+1.  Z0 = Zl.  Z1 = Zl+1.
    n is the length  of Z0
    """
    function calculateFirstZ(Z0, n)
    
        Z1 = zeros(Float64, n)
    
        # normalize the chebyshev result
        for i in eachindex(Z0)
            Z0[i] *= D0;
        end
    
        # now figure out Z1
    
        for i=2:n-2
            Z1[i] = AL(0)*(1/ALPHA*Z0[i+1]-BETA/ALPHA*Z0[i]-GAMMA/ALPHA*Z0[i-1]) + BL(0)*Z0[i];
        end

        return Z1
    
    end
    
    """
     performs a hermite transform in the most naive way possible directly from the
         data points given in xl()
    """
    function naiveTransform(data, results)
    
        results = zeros(BIGN)
    
        for x=1:2*LITN # for each data point
            Lminus1 = 0;
            curVal = D0;
            for y in 1:BIGN # go through the n hermites
                results[y] += (data[x])*curVal;
                Lminus2 = Lminus1;
                Lminus1 = curVal;
                curVal  = (AL(y)*(xk(x)) + BL(y))*Lminus1 + CL(y)*Lminus2;
            end
        end
    end
    
    """
    perform a chebyshev transform in the most naive way possible 
    directly from the
    data points defined by xl()
    """
    function chebyshev(data)
    
       results = zeros(ComplexF64, BIGN)
       
       for x in eachindex(data)
           Lminus1 = xk(x-1)/BIGC;          # x
           Lminus2 = 2.0 * Lminus1^2 - 1;   # 2*x^2-1
           for y in 1:BIGN                  # go through the n chebyshevs
               curVal = (ALPHA*(xk(x-1)) +BETA)*Lminus1 + GAMMA*Lminus2;
               results[y] += data[x] * curVal;
               Lminus2=Lminus1;
               Lminus1=curVal;
           end
       end

       return results
    
    end

    """
    a recursive function which performs a Hermite transform in O(n(logn)^2) time
    Z0 and Z1 must be precomputed as defined in the paper.  l should be first
      set to 1
    you must precompute all the necessary Rns.
    """
    function performTransform( rns, Z0, Z1, n, l, result)

         result[l-1] = Z0[n];
         result[l]   = Z1[n];
    
         if (n<3) return end
    
         #temp to store the new data
         temp = zeros(Float64, 4n);
    
         # combine Z0 and Z1 into Z to get ready for the matrix multiply
         
         Z = vcat(Z0, Z1)

         preFourBcirculantVcMatrixMultiply(rns, n, l-1, Z, temp);
    
         nover2 = n÷2;

         performTransform(rns, Z0+nover2,Z1+nover2,nover2,l,result);

         performTransform(rns, temp+nover2,temp+5*nover2,nover2,l+nover2,result);
         
         return

    end
    
    
    function oneDTransform(data, result)
    
        n = BIGN;
    
        Z0    = zeros( ComplexF64, 2n)
        dblZ0 = zeros( ComplexF64, 2n)
        Z1    = zeros( ComplexF64, 2n)
        dblZ1 = zeros( ComplexF64, 2n)
    
        Z0[n:end-1] .= chebyshev(data);

        Z0[end] = 0.0
    
        # we only want the real parts
        for i=0:n-1
            dblZ0[n+i] = real(Z0[n+i])
        end
    
        # expand the data
        for i=0:n-1
            dblZ0[i+1] = dblZ0[2n-i-1];
        end

        
        # find the next data point
        dblZ1 .= calculateFirstZ(dblZ0, 2n);

        rns = Rns(n)

        # do the second part
        performTransform(rns, dblZ0, dblZ1, n, 2, result);

    
    end


end 
