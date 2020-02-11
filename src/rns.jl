struct Rns

    data :: Vector{Float64}
    size :: Int

    function Rns( n :: Int)
    
        da = zeros(Float64, (n+1) * n * 8n);
        da_size = n;
        # precompute the necessary Rns

        i = n
        while i>=4 
            j = 0
            while j<n
                precompute(i, j, da[8n*(n*i+j)]);
                j += i 
            end
            i = i ÷ 2  
        end


        new( da, da_size)

    end

end

"""
define An(l) as defined in the paper
result is a 8*n sized vector with each 2n representing a circulant block
- top left is all zeros (we'll zero out everything else while we're at it.
- top right is I2n
- bottom left is cl*I2n
- bottom right is Cn(wl,vl,ul);
"""
function createAn(n, l) 

    result = zeros(Float64, 8n)
    result[2n+1] = 1;
    result[4n+1] = CL(l);
    result[6n+1] = VL(l);
    result[6n+2] = WL(l);
    result[8n]   = UL(l);

    return result

end

"""
 computes desired Rn
 needed columns of circulant matrices listed vertically as top left, top right, bottom left, bottom right.
"""
function precompute( n, l, result) 

    @show n, l
    temp  = zeros(8*n);
    temp2 = zeros(8*n);

    result = createAn(n, n÷2+l );


    for i=l+n÷2-1:-1:l
         temp = createAn(n,i);
         fourBcirculantSqMatrixMultiply(result, temp, 4n, temp2);
         result .= temp2
    end

    println("***  createAn result *** ")
    for i in eachindex(result) 
        if (result[i] != 0) 
            println(i, " ", round(result[i], digits=7)) 
        end
    end
    exit()
    
   
end



# give me the first column of a circulant matrix in M.
function circulantVcMatrixMultiply( c, VecCpy, n, result)

     println("circulantVcMatrixMultiply")

     @show size(c)
     @show size(VecCpy)
     @show n
     @show size(result)
     fftc     = rfft(c)
     fftVec   = rfft(VecCpy)
     multiply = fftc .* fftVec
     @show size(multiply)
     result  .= irfft(multiply, length(c))

end

"""
A B  *  E F  = AE+BG AF+BH
C D     G H    CE+DG CF+DH
"""
function fourBcirculantSqMatrixMultiply( M1, M2, n, result) 

    println("doing a multiplication of size $n")

    temp1 = zeros(n÷2);
    temp2 = zeros(n÷2);

    #fill up the columns
    A = M1
    E = M2
    B = M1[n÷2+1:end]
    F = M2[n÷2+1:end]
    C = M1[n+1:end]
    G = M2[n+1:end]
    D = M1[3n÷2+1:end]
    H = M2[3n÷2+1:end]

    # A*E+B*G top left
    circulantVcMatrixMultiply(A,E,n÷2,temp1);
    circulantVcMatrixMultiply(B,G,n÷2,temp2);
    # Add em up
    for i=1:n÷2
        result[i]=temp1[i]+temp2[i];
    end

    # A*F+B*H top right
    circulantVcMatrixMultiply(A,F,n÷2,temp1);
    circulantVcMatrixMultiply(B,H,n÷2,temp2);
    # Add em up
    for i=1:n÷2
        result[i+n÷2]=temp1[i]+temp2[i];
    end

    # C*E+D*G bottom left
    circulantVcMatrixMultiply(C,E,n÷2,temp1);
    circulantVcMatrixMultiply(D,G,n÷2,temp2);
    # Add em up
    for i=1:n÷2
        result[i+n]=temp1[i]+temp2[i];
    end

    # C*F+D*H bottom right
    circulantVcMatrixMultiply(C,F,n÷2,temp1);
    circulantVcMatrixMultiply(D,H,n÷2,temp2);
    # Add em up
    for i=1:n÷2
        result[i+3n÷2]=temp1[i]+temp2[i];
    end


end




# multiply Z by the Rn that was precomputed at n,l
function preFourBcirculantVcMatrixMultiply( rns, n, l, Vec, result)

    n *= 4;

    # top left (want a column of the top left times first half of Vec)
    circulantVcMatrixMultiply(rns.data[rns.size*8*(rns.size*n/4+l)],Vec,n/2,result);

    # top right (want a column of the top right times second half of Vec)
    circulantVcMatrixMultiply(rns.data[rns.size*8*(rns.size*n/4+l)+n/2],Vec+n/2,n/2,temp);

    # add top left and top right
    for x=0:n/2
        result[x] += temp[x];
    end

    # bottom left (want a column of the bottom left times first half of Vec)
    circulantVcMatrixMultiply(rns.data[rns.size*8*(rns.size*n/4+l)+n],Vec,n/2,result+n/2);

    # bottom right (want a column of the bottom right times second half of Vec)
    circulantVcMatrixMultiply(rns.data[rns.size*8*(rns.size*n/4+l)+3*n/2],Vec+n/2,n/2,temp2);

    # add bottom left and bottom right
    for x=n/2:n
        result[x] += temp2[x-n/2];
    end

    n /= 4

end
