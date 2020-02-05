module FastHermiteTransform

const BIGN = 128
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


"""
defines the data points
"""
fftw_complex xk(k) = ((k-LITN)/LITN)*BIGC


"""
use Zl to calculate Zl+1.  Z0 = Zl.  Z1 = Zl+1.
n is the length  of Z0
"""
function calculateFirstZ(Z0, Z1, n)

    Z1 = zeros(Float64, n)

    # normalize the chebyshev result
    for i in eachindex(Z0)
	Z0[i] *= D0;
    end

    # now figure out Z1

    for i=2:n-1
        Z1[i] = AL(0)*(1/ALPHA*Z0[i+1]-BETA/ALPHA*Z0[i]-GAMMA/ALPHA*Z0[i-1]) + BL(0)*Z0[i];
    end
}
end

"""
define An(l) as defined in the paper
result is a 8*n sized vector with each 2n representing a circulant block
"""
function createAn(n, l, result) 

    # top left is all zeros (we'll zero out everything 
    # else while we're at it.
    result = zeros(Float64, 8*n)
    # top right is I2n
    result[2*n]=1;
    # bottom left is cl*I2n
    result[4*n]=CL(l);
    # bottom right is Cn(wl,vl,ul);
    result[6*n]   = VL(l);
    result[6*n+1] = WL(l);
    result[8*n-1] = UL(l);

end


"""
//stores desired Rn in a file /Rns/xxxxx_xxxxx.dat if file does not already exist;
//needed columns of circulant matrcies listed vertically as top left, top right, bottom left, bottom right.
"""
function precomputeRnAndStore( n, l, result) 

    temp  = zeros(8*n);
    temp2 = zeros(8*n);

    createAn(n,n/2+l,result);
    for i=l+n÷2-1:-1;l
         createAn(n,i,temp);
         fourBcirculantSqMatrixMultiply(result,temp,4*n,temp2);
         result .= temp2
    end
    
   
end

"""
perform a chebyshev transform in the most naive way possible 
directly from the
data points defined by xl()
"""
function naiveChebyshev(data, results)

results = zeros(ComplexF64, BIGN)

#fftw_complex Lminus1;
#fftw_complex Lminus2;
#fftw_complex curVal;
#int x,y;

for(x=0; x<=2*LITN; ++x) #for each data point
    Lminus1 = xk(x)/BIGC;	# x
    Lminus2 = 2.0*pow(Lminus1,2)-1;  # 2*x^2-1
    for(y=0; y<BIGN; ++y) # go through the n chebyshevs
        curVal = (ALPHA*(xk(x)) +BETA)*Lminus1 + GAMMA*Lminus2;
results[y] += ((fftw_complex)data[x])*curVal;
Lminus2=Lminus1;
Lminus1=curVal;
}
}
end
#=

//a recursive function which performs a Hermite transform in O(n(logn)^2) time
//Z0 and Z1 must be precomputed as defined in the paper.  l should be first
//  set to 1
//you must precompute all the necessary Rns.
void performTransform(double* Z0, double* Z1, int n, int l, double* result) {
     result[l-1] = Z0[n-1];
     result[l]   = Z1[n-1];

     if(n<3) return;

     //temp to store the new data
     double* temp = (double *) fftw_malloc(sizeof(double) * 4*n);


     //combine Z0 and Z1 into Z to get ready for the matrix multiply
     double *Z;
     Z = (double *) fftw_malloc(sizeof(double) * 4*n);
     memcpy(Z,Z0,sizeof(double)*2*n);
	 memcpy(Z+n*2,Z1,sizeof(double)*2*n);
     preFourBcirculantVcMatrixMultiply(n,l-1,Z,temp);

	 fftw_free(Z);
	 int nover2 = n/2;
     performTransform(Z0+nover2,Z1+nover2,nover2,l,result);
     performTransform(temp+nover2,temp+5*nover2,nover2,l+nover2,result);
	 fftw_free(temp);
     return;
}

//performs a hermite transform in the most naive way possible directly from the
//	data points given in xl()
void naiveTransform(double *data, double *results) {
	memset(results,0,sizeof(double)*BIGN);
	fftw_complex Lminus1;
	fftw_complex Lminus2;
	fftw_complex curVal;
//	double *tj = (double*) fftw_malloc(sizeof(double) * BIGN * (2*LITN+1));
	int x,y;
	for(x=0; x<=2*LITN; ++x) {//for each data point
	    Lminus1 = 0;
	    curVal = D0;
		for(y=0; y<BIGN; ++y) {//go through the n hermites
//			tj[y*(2*LITN+1)+x]=curVal;
			results[y] += (data[x])*curVal;
			Lminus2=Lminus1;
			Lminus1=curVal;
			curVal = (AL(y)*(xk(x)) + BL(y))*Lminus1 + CL(y)*Lminus2;
		}
	}
//    printNonSqMatrix(tj,BIGN,2*LITN+1);
}

void oneDTransform(double *data, double* result) {
	int i, j;
	int n=BIGN;

	fftw_complex *Z0 =  (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (2 * n));
	double *dblZ0 	 = 	(double *) fftw_malloc(sizeof(fftw_complex) * (2 * n));
	fftw_complex *Z1 = 	(fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2 * n);
	double *dblZ1 	 = 	(double *) fftw_malloc(sizeof(fftw_complex) * (2 * n));

 //do a chebyshev Transform
	naiveChebyshev(data,Z0+n-1);
//	printVector("cheby",Z0+n-1,n);
	Z0[2*n-1]=0;

	//we only want the real parts
	for(i=0; i<n; ++i) {
		dblZ0[n-1+i] = creal(Z0[n-1+i]);
	}

	//expand the data
	for(i=0; i<n; ++i)
		dblZ0[i]=dblZ0[2*n-i-2];

	//find the next data point
	calculateFirstZ(dblZ0,dblZ1,2*n);

	//do the second part
	performTransform(dblZ0,dblZ1,n,1,result);

	fftw_free(Z0);
	fftw_free(dblZ0);
	fftw_free(Z1);
	fftw_free(dblZ1);
}

int main(int argc, char *argv[])
{
    int n = LITN;
	int N = BIGN;
    int i,j,x;
    tester = fopen("doc/testing.txt", "wt+");


 	double *data   = (double *) fftw_malloc(sizeof(double) * (n*2+1));

    double *fancyResult = (double *) fftw_malloc(sizeof(double) * N);
	double *naiveResult = (double *) fftw_malloc(sizeof(double) * N);

    FILE *inputFile;
    inputFile = fopen("doc/data.txt", "wt+");

	srand(time(0));
	double daConstant = pow(-1,rand()%2)*rand()/RAND_MAX;
	for(i=0; i<=2*n; ++i) {
        fprintf(inputFile,"%Lf	",(double)rand()/RAND_MAX);//daConstant*1/pow(i+1,2));
//		for(j=0; j<(n/4); ++j)
//	        fprintf(inputFile,"%Lf	",(double)0);
//		for(j=n/4; j<(3*n/4); ++j)
//	   		fprintf(inputFile,"%Lf	",(double)1);
//		for(j=(3*n/4);j<n;++j)
//		    fprintf(inputFile,"%Lf	",(double)0);
		fprintf(inputFile,"\n");
	}

    fclose(inputFile);
	oneDFileReader("doc/data.txt",(2*n+1),data);

    double* diag = (double *) fftw_malloc(sizeof(double) * (2*n+1)*(2*n+1));
   	memset(diag,0,sizeof(double) * (2*n+1)*(2*n+1));
	for(i=0;i<=2*n; ++i) {
		diag[(2*n+1)*i+i] = exp(0-pow(xk(i),2)/2);
        data[i] *= diag[(2*n+1)*i+i];
	}


    initFastFouriers(N);
	initRns(N);
	fclock=clock();
	oneDTransform(data,fancyResult);
	fclock-=clock();

	nclock=clock();
	naiveTransform(data,naiveResult);
	nclock-=clock();

	FILE *results;
	results=fopen("doc/results.txt","wt+");
	for(i=0;i<N;++i){
		fancyResult[i]*=(2*BIGC)/n;
        naiveResult[i]*=(2*BIGC)/n;
		fprintf(results,"%.16Lf\n",fancyResult[i]);
//        fprintf(results,"%.16Lf\n\n",naiveResult[i]);
	}
	fclose(results);
    fprintf(stdout,"Fancy took %i \nNaive took %i \n to a degree of %ld per second\n\n\n",-fclock,-nclock, CLOCKS_PER_SEC);

    system("PAUSE");
    return 0;
}

=#
end 
