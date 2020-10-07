#ifndef DATA_HPP
#define DATA_HPP
#include <iostream>

#include <algorithm>
#include <time.h>
#include <vector>
#include <math.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <sys/time.h>



double squared(double* x, double* expectation, double** variance,int nvar);
double rvtdt_brv(int** haps, int** O, int** P, int nfam, int nvar);
void gtdt(int** haps, int** O, int** P, int nfam, int nvar, double* s_ad, double* s_dom, double* s_ch);
double sv_fbat(double* x, double* expectation, double** variance,int nvar, int index);
double burden(double* x, double* expectation, double** variance,int nvar);
double hc(double* x, double* expectation, double** variance,int nvar);
double max_sv(double* x, double* expectation, double** variance,int nvar);

void	i_destroyMatrix(int** matrix, int rows);
int**	i_createMatrix(int sizeX, int sizeY);
void	d_destroyMatrix(double** matrix, int rows);
double**	d_createMatrix(int sizeX, int sizeY);
void	i_destroyVector(int* vector);
int*	i_createVector(int size);
void	d_destroyVector(double* vector);
double*	d_createVector(int size);

struct MersenneTwister
{
  /* Period parameters */
  const static int N = 624;
  const static int M = 397;
  const static unsigned long long MATRIX_A = 0x9908b0dfUL;  /* constant vector a */
  const static unsigned long long UPPER_MASK = 0x80000000UL;  /* most significant w-r bits */
  const static unsigned long long LOWER_MASK = 0x7fffffffUL;  /* least significant r bits */

  unsigned long mt[N]; /* the array for the state vector  */
  int mti; /* mti==N+1 means mt[N] is not initialized */

  MersenneTwister()
  {
    mti = N + 1;
    unsigned long init[4] = { 0x123, 0x234, 0x345, 0x456 }, length = 4;
    init_by_array( init, length );
  }

  MersenneTwister( unsigned int seed )
  {
    mti = N + 1;
    init_genrand( seed );
  }

private:
  /* initializes mt[N] with a seed */
  void init_genrand( unsigned long s )
  {
    mt[0] = s & 0xffffffffUL;
    for( mti = 1; mti < N; mti++ ) {
      mt[mti] =
        (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      mt[mti] &= 0xffffffffUL;
      /* for >32 bit machines */
    }
  }

  /* initialize by an array with array-length */
  /* init_key is the array for initializing keys */
  /* key_length is its length */
  /* slight change for C++, 2004/2/26 */
  void init_by_array( unsigned long init_key[], int key_length )
  {
    int i, j, k;
    init_genrand( 19650218UL );
    i = 1; j = 0;
    k = (N > key_length ? N : key_length);
    for( ; k; k-- ) {
      mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525UL))
        + init_key[j] + j; /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++; j++;
      if( i >= N ) { mt[0] = mt[N - 1]; i = 1; }
      if( j >= key_length ) j = 0;
    }
    for( k = N - 1; k; k-- ) {
      mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941UL))
        - i; /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      if( i >= N ) { mt[0] = mt[N - 1]; i = 1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
  }

  /* generates a random number on [0,0xffffffff]-interval */
  unsigned long randUInt()
  {
    unsigned long y;
    static unsigned long mag01[2] = { 0x0UL, MATRIX_A };
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if( mti >= N ) { /* generate N words at one time */
      int kk;

      if( mti == N + 1 )   /* if init_genrand() has not been called, */
        init_genrand( 5489UL ); /* a default initial seed is used */

      for( kk = 0; kk < N - M; kk++ ) {//227 iterations
        y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
        mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      for( ; kk < N - 1; kk++ ) {//goes the 397 remaining iterations (624 iterations total)
        y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
        mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
      }
      y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

      mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
  }

public:
  /* generates a random number on [0,0x7fffffff]-interval */
  // (0 and 2147483647)
  inline int randInt() { return randUInt() >> 1; }
  inline float randFloat() { return randInt()*(1.f / 4294967295.f); }
  inline double randDouble() { return randInt()*(1.0 / 4294967295.0); }
};
#endif
