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
#include <string.h>

#define BIGN 4
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

//****************************************************************************
//==========================LINEAR ALGEBRA TOOLS==============================
//****************************************************************************
FILE *tester;
fftw_plan *daPlans;
double  *daRns;
int daRnsSize;

clock_t fclock, nclock;

//defines the data points
fftw_complex xk(int k) {
    return ((double)(k-LITN)/(double)LITN)*BIGC;
}


void initFastFouriers(int n) {
    daPlans = (fftw_plan*) fftw_malloc(sizeof(fftw_plan) * 2*n+1);

    int i;
    for(i=8; i<=2*n; i*=2) {
        double* in = (double*) fftw_malloc(sizeof(double) * i);
        fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (i/2+1));
        daPlans[i] = fftw_plan_dft_r2c_1d(i,in,out,FFTW_ESTIMATE);
        fftw_free(in);
        fftw_free(out);
    }
}

//define An(l) as defined in the paper
//result is a 8*n sized vector with each 2n representing a circulant block
void createAn(int n, int l, double *result) {
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

//give me the first column of a circulant matrix in M.
void circulantVcMatrixMultiply(double* c, double* VecCpy, int n, double* result) {
     int i;
     fftw_complex* fftc = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n);
     fftw_complex* fftVec = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n);

     //fast fourier c
     //     fftw_plan cplan = fftw_plan_dft_r2c_1d(n,c,fftc,FFTW_ESTIMATE);
     fftw_execute_dft_r2c(daPlans[n],c,fftc);

     //fast fourier Vec
     fftw_execute_dft_r2c(daPlans[n],VecCpy,fftVec);

     fftw_complex* multiply = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * n);
     for(i=0; i<n; ++i) {
         multiply[i] = fftc[i]*fftVec[i];
     }

     fftw_plan Finalplan  = fftw_plan_dft_c2r_1d(n,multiply,result,FFTW_ESTIMATE);
     fftw_execute(Finalplan);

     for(i=0; i<n; ++i) {
              result[i]/=n;
     }

     fftw_destroy_plan(Finalplan);
     fftw_free(fftVec);
     fftw_free(multiply);
     fftw_free(fftc);

     return;
}

// A B  *  E F  = AE+BG AF+BH
// C D     G H    CE+DG CF+DH
void fourBcirculantSqMatrixMultiply(double* M1, double* M2, int n, double* result) {
    printf("doing a multiplication of size %i\n",n);
    double *A,*B,*C,*D,*E,*F,*G,*H;
    double* temp1      = (double *) fftw_malloc(sizeof(double) * n/2);
    double* temp2      = (double *) fftw_malloc(sizeof(double) * n/2);
    int i;

    //fill up the columns
    A = M1;
    E = M2;
    B = M1+n/2;
    F = M2+n/2;
    C = M1+n;
    G = M2+n;
    D = M1+3*n/2;
    H = M2+3*n/2;

    //A*E+B*G top left
    circulantVcMatrixMultiply(A,E,n/2,temp1);
    circulantVcMatrixMultiply(B,G,n/2,temp2);
    //Add em up
    for(i=0;i<n/2;++i)
        result[i]=temp1[i]+temp2[i];

    //A*F+B*H top right
    circulantVcMatrixMultiply(A,F,n/2,temp1);
    circulantVcMatrixMultiply(B,H,n/2,temp2);
    //Add em up
    for(i=0;i<n/2;++i)
        result[i+n/2]=temp1[i]+temp2[i];

    //C*E+D*G bottom left
    circulantVcMatrixMultiply(C,E,n/2,temp1);
    circulantVcMatrixMultiply(D,G,n/2,temp2);
    //Add em up
    for(i=0;i<n/2;++i)
        result[i+n]=temp1[i]+temp2[i];

    //C*F+D*H bottom right
    circulantVcMatrixMultiply(C,F,n/2,temp1);
    circulantVcMatrixMultiply(D,H,n/2,temp2);
    //Add em up
    for(i=0;i<n/2;++i)
        result[i+3*n/2]=temp1[i]+temp2[i];

    fftw_free(temp1);
    fftw_free(temp2);
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
    fftw_free(temp);
    fftw_free(temp2);

}



void initRns(int n) {
    int i,j;
    daRns         = (double*) fftw_malloc(sizeof(double) * (n+1) * n * 8*n);
    daRnsSize     = n;
    //precompute the necessary Rns
    for(i=n; i>=4;i/=2) {
        for(j=0; j<n; j+=i) {
            precomputeRnAndStore(i, j, &daRns[8*n*(n*i+j)]);
        }
    }

}


void destroyFastFouriers(int n) {
    int i;
    for(i=8; i<=2*n; i*=2) {
        fftw_destroy_plan(daPlans[i]);
    }
    fftw_free(daPlans);
}



//multiply Z by the Rn that was precomputed at n,l
void preFourBcirculantVcMatrixMultiply(int n, int l, double* Vec, double* result) {
    n*=4;
    double* temp = (double *) fftw_malloc(sizeof(double) * n/2);
    double* temp2 = (double *) fftw_malloc(sizeof(double) * n/2);

    int x;
    //top left (want a column of the top left times first half of Vec)
    circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*n/4+l)],Vec,n/2,result);

    //top right (want a column of the top right times second half of Vec)
    circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*n/4+l)+n/2],Vec+n/2,n/2,temp);

    //add top left and top right
    for(x=0; x<n/2; ++x)
        result[x] += temp[x];

    //bottom left (want a column of the bottom left times first half of Vec)
    circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*n/4+l)+n],Vec,n/2,result+n/2);

    //bottom right (want a column of the bottom right times second half of Vec)
    circulantVcMatrixMultiply(&daRns[daRnsSize*8*(daRnsSize*n/4+l)+3*n/2],Vec+n/2,n/2,temp2);

    //add bottom left and bottom right
    for(x=n/2; x<n; ++x)
        result[x] += temp2[x-n/2];

    fftw_free(temp);
    fftw_free(temp2);
    n/=4;
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


int main(int argc, char *argv[])
{
    int n = LITN;
    int N = BIGN;

    double *data   = (double *) fftw_malloc(sizeof(double) * (n*2+1));

    double *result = (double *) fftw_malloc(sizeof(double) * N);

    for(int i=0;i<2*n+1; ++i) data[i] = exp(0-pow(xk(i),2)/2);


    initFastFouriers(N);
    initRns(N);

    fftw_complex *Z0 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * (2 * N));
    double *dblZ0 = (double *) fftw_malloc(sizeof(fftw_complex) * (2 * N));
    fftw_complex *Z1 = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * 2 * N);
    double *dblZ1 = (double *) fftw_malloc(sizeof(fftw_complex) * (2 * N));

    memset(Z0,    0, sizeof(fftw_complex) * (2*N));
    memset(dblZ0, 0, sizeof(fftw_complex) * (2*N));
    memset(Z1,    0, sizeof(fftw_complex) * (2*N));
    memset(dblZ1, 0, sizeof(fftw_complex) * (2*N));

    //perform a chebyshev transform in the most naive way possible directly from the
    //    data points defined by xl()

    fftw_complex Lminus1;
    fftw_complex Lminus2;
    fftw_complex curVal;
    for(int x=0; x<=2*LITN; ++x) {//for each data point
        Lminus1 = (double)xk(x)/BIGC;    //x
        Lminus2 = 2.0*pow(Lminus1,2)-1;  //2*x^2-1
        for(int y=0; y<BIGN; ++y) {//go through the n chebyshevs
            curVal = (ALPHA*(xk(x)) +BETA)*Lminus1 + GAMMA*Lminus2;
            Z0[y+N-1] += ((fftw_complex)data[x])*curVal;
            Lminus2=Lminus1;
            Lminus1=curVal;
        }
    }

    // end chebyshev


    Z0[2*N-1]=0;

    //we only want the real parts
    for(int i=0; i<N; ++i) {
        dblZ0[N-1+i] = creal(Z0[N-1+i]);
    }

    //expand the data
    for(int i=0; i<N; ++i)
        dblZ0[i]=dblZ0[2*N-i-2];


    //find the next data point
    calculateFirstZ(dblZ0,dblZ1,2*N);

    //do the second part
    performTransform(dblZ0,dblZ1,N,1,result);

    fftw_free(Z0);
    fftw_free(dblZ0);
    fftw_free(Z1);
    fftw_free(dblZ1);

    for(int i=0;i<N;++i){
        result[i]*=(2*BIGC)/n;
        printf("%16.7lf\n", result[i]);
    }

    return 0;
}

