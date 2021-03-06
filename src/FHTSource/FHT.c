//******************************************
//NAME: Robert Taintor
//DATE: 8/3/07
//PURPOSE: Performs a fast polynomial transform on a group of orthogonal
//            polynomials which satisfy a 3 term recurrance
//            from uniformly sampled points.

#include <complex.h> //I is now defined as sqrt(-1)
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <testingstuff.h>
#include <lalgebra.h>
#include <string.h>

#define BIGN 8
#define LITN (BIGN*5)
#define BIGC sqrt((double)2*BIGN+1)

#define ALPHA (2.0/BIGC)
#define BETA 0
#define GAMMA (-1.0)
#define PI 3.141592654
#define D0 (double)1.0/pow(PI,(double)1.0/4.0)

//define al, bl, cl, ul, vl, and wl
#define AL(l) (double)sqrt((double)2.0/(l+1.0))
#define BL(l) (double)0
#define CL(l) (double)((-1.0)*sqrt((double)l/(l+1.0)))

#define UL(l) (double)(AL(l)/ALPHA)
#define VL(l) (double)(BL(l)+(-1)*(UL(l)*BETA))
#define WL(l) (double)((-1)*GAMMA*UL(l))


clock_t fclock, nclock;

//defines the data points
fftw_complex xk(int k) {
    return ((double)(k-LITN)/(double)LITN)*BIGC;
}

// use Zl to calculate Zl+1.  Z0 = Zl.  Z1 = Zl+1.
// n is the length  of Z0
void calculateFirstZ(double *Z0, double *Z1,int n)
{
    memset(Z1,0,sizeof(double) * n);
    int i;
    //normalize the chebyshev result
    for(i=0; i<n; ++i)
        Z0[i]*=D0;
    //now figure out Z1
    for(i=1; i<(n-2); ++i) {
        Z1[i] = AL(0)*(1/ALPHA*Z0[i+1]-BETA/ALPHA*Z0[i]-GAMMA/ALPHA*Z0[i-1]) + BL(0)*Z0[i];
    }
}

//define An(l) as defined in the paper
//result is a 8*n sized vector with each 2n representing a circulant block
void createAn(int n, int l, double *result) {
    int i;
    //top left is all zeros (we'll zero out everything else while we're at it.
    memset(result,0,sizeof(double)*8*n);
    //top right is I2n
    result[2*n]   = 1;
    //bottom left is cl*I2n
    result[4*n]   = CL(l);
    //bottom right is Cn(wl,vl,ul);
    result[6*n]   = VL(l);
    result[6*n+1] = WL(l);
    result[8*n-1] = UL(l);
}

//stores desired Rn in a file /Rns/xxxxx_xxxxx.dat if file does not already exist;
//needed columns of circulant matrcies listed vertically as top left, top right, bottom left, bottom right.
void precomputeRnAndStore(int n, int l, double* result) {
    int i;

    //now it's not so we open our write file and start computing
    double *temp = (double *) fftw_malloc(sizeof(double) * 8*n);
    double *temp2 = (double *) fftw_malloc(sizeof(double) * 8*n);

    //printf(" l = %d \n ", l);

    createAn(n,n/2+l,result);


    for(i=l+n/2-1; i>l; --i) {
         createAn(n,i,temp);
         fourBcirculantSqMatrixMultiply(result,temp,4*n,temp2);
         memcpy(result,temp2,sizeof(double)*8*n);
    }
    //for(i=0;i<=8*n;++i){
        //if (result[i] != 0.0) printf("%d %.16lf\n",i+1, result[i]);
    //}
    fftw_free(temp);
    fftw_free(temp2);

}


//perform a chebyshev transform in the most naive way possible directly from the
//    data points defined by xl()
void naiveChebyshev(double *data, fftw_complex *results) {
    memset(results,0,sizeof(fftw_complex)*BIGN);
    fftw_complex Lminus1;
    fftw_complex Lminus2;
    fftw_complex curVal;
    int x,y;
    for(x=0; x<=2*LITN; ++x) {//for each data point
        Lminus1 = (double)xk(x)/BIGC;    //x
        Lminus2 = 2.0*pow(Lminus1,2)-1;  //2*x^2-1
        for(y=0; y<BIGN; ++y) {//go through the n chebyshevs
            curVal = (ALPHA*(xk(x)) +BETA)*Lminus1 + GAMMA*Lminus2;
            results[y] += ((fftw_complex)data[x])*curVal;
            Lminus2=Lminus1;
            Lminus1=curVal;
        }
    }

    //for(int i; i<BIGN; ++i) printf("%lf+%lf\n", creal(results[i]), cimag(results[i]));
}

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
//    data points given in xl()
void naiveTransform(double *data, double *results) {
    memset(results,0,sizeof(double)*BIGN);
    fftw_complex Lminus1;
    fftw_complex Lminus2;
    fftw_complex curVal;
//    double *tj = (double*) fftw_malloc(sizeof(double) * BIGN * (2*LITN+1));
    int x,y;
    for(x=0; x<=2*LITN; ++x) {//for each data point
        Lminus1 = 0;
        curVal = D0;
        for(y=0; y<BIGN; ++y) {//go through the n hermites
//            tj[y*(2*LITN+1)+x]=curVal;
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

    fftw_complex *Z0 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (2 * n));
    double *dblZ0 = (double *) fftw_malloc(sizeof(fftw_complex) * (2 * n));
    fftw_complex *Z1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2 * n);
    double *dblZ1 = (double *) fftw_malloc(sizeof(fftw_complex) * (2 * n));

    memset(Z0,    0, sizeof(fftw_complex) * (2*n));
    memset(dblZ0, 0, sizeof(fftw_complex) * (2*n));
    memset(Z1,    0, sizeof(fftw_complex) * (2*n));
    memset(dblZ1, 0, sizeof(fftw_complex) * (2*n));

    //do a chebyshev Transform
    naiveChebyshev(data,Z0+n-1);


    //printVector("cheby",Z0+n-1,n);
    //exit(0);
    Z0[2*n-1]=0;

    //printf("\n Z0 \n");
    //for(int i; i < 2*n; ++i) printf("%d %lf+%lf \n", i+1, creal(Z0[i]), cimag(Z0[i]));

    //we only want the real parts
    for(i=0; i<n; ++i) {
        dblZ0[n-1+i] = creal(Z0[n-1+i]);
    }

    //expand the data
    for(i=0; i<n; ++i)
        dblZ0[i]=dblZ0[2*n-i-2];


    //find the next data point
    calculateFirstZ(dblZ0,dblZ1,2*n);

    // printf("\n dblZ1 \n");
    // for(int i; i < 2*n; ++i) printf("%d %lf+%lf \n", i+1, creal(dblZ1[i]), cimag(dblZ1[i]));
    // exit(0);

    //do the second part
    performTransform(dblZ0,dblZ1,n,1,result);

    fftw_free(Z0);
    fftw_free(dblZ0);
    fftw_free(Z1);
    fftw_free(dblZ1);
}

void twoDTransform(double *data, int n, double* result) {
    int i, j;

    //initialize the plans for the Fourier Transforms
    initFastFouriers(n);
    initRns(n);
    //transform the rows
    for(i=0; i<=2*n; ++i) {
        oneDTransform(&data[n*i],&result[n*i]);
    }
    //transform the columns
    double* column = (double *) fftw_malloc(sizeof(double) * n);
    double* temp = (double *) fftw_malloc(sizeof(double) * n);
    for(i=0; i<n; ++i) {
        getColumnFromSq(result,column,i,n);
        oneDTransform(column,temp);
        setColumnToSq(temp,result,i,n);
    }
    fftw_free(column);
    fftw_free(temp);

    destroyFastFouriers(n);
}

int main(int argc, char *argv[])
{
    int n = LITN;
    int N = BIGN;
    int i,j,x;

    double *data   = (double *) fftw_malloc(sizeof(double) * (n*2+1));

    double *fancyResult = (double *) fftw_malloc(sizeof(double) * N);
    double *naiveResult = (double *) fftw_malloc(sizeof(double) * N);

    // for(int j=0;      j<n/2;   ++j) data[j] = 0.0;
    // for(int j=n/2;    j<=3*n/2; ++j) data[j] = 1.0;
    // for(int j=(3*n/2+1);j<2*n+1; ++j) data[j] = 0.0;

    for(int j=0;j<2*n+1; ++j) data[j] = 1.0;

    FILE *inputFile;
    inputFile = fopen("doc/data.txt", "wt+");
    for(i=0; i<=2*n; ++i) {
        fprintf(inputFile, "%lf", data[i]);
	fprintf(inputFile,"\n");
    }
    fclose(inputFile);

    double* diag = (double *) fftw_malloc(sizeof(double) * (2*n+1)*(2*n+1));
    memset(diag,0,sizeof(double) * (2*n+1)*(2*n+1));
    for(i=0;i<=2*n; ++i) {
        diag[(2*n+1)*i+i] = exp(0-pow(xk(i),2)/2);
        data[i] *= diag[(2*n+1)*i+i];
    }

    initFastFouriers(N);
    initRns(N);
    fclock=clock();
    oneDTransform(data, fancyResult);
    fclock-=clock();

    nclock=clock();
    naiveTransform(data,naiveResult);
    nclock-=clock();
    for(i=0;i<N;++i){
        fancyResult[i]*=(2*BIGC)/n;
        naiveResult[i]*=(2*BIGC)/n;
        //printf("%16.3lf\n",fabs(fancyResult[i]-naiveResult[i]));
        printf("%16.7lf\n",naiveResult[i]);
    }
    printf("Fancy took %i \nNaive took %i \n to a degree of %ld per second\n\n\n",-fclock,-nclock, CLOCKS_PER_SEC);

    return 0;
}
