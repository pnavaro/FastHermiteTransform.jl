# The FHT package

A Fast Hermite Transform on the uniform grid

FHT is a c program that computes a 1D or 2D Fast Hermite Transform. The algorithim used is a modifcation of the 1995 work Fast Discrete Polynomial Transforms with Applications to Data Analysis for distance transitive graphs by Driscoll, Healy, and Rockmore.

The program was compiled under a windows enviroment using Dev C++ Version 4.9.9.2 and the GNU C compiler. I have not tested this under a Linux or Macintosh enviroment, but I believe it can be easily adapted to either. This code is by no means a polished package, but a source code meant to be adapted for a particular purpose.

The code as given performs a 1D transform from 1281 sample points to the first 128 Hermite polynomials. It is easily adaptable to perform both the 1D and 2D case for many different sizes.

Also included is a testing package called testingStuff which I used to print matrices and vectors in the format used by the program and the FFTW fast fourier software package that was utilized in performing the fast vector multiplications.

# Copyright

Dan Rockmore Rockmore@cs.dartmouth.edu, Gregory LeibonGregory.Leibon@dartmouth.edu, or Robert Taintor Robert.C.Taintor@gmail.com.
