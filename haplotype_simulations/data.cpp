#include "data.hpp"

#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

#ifndef HUGEDBL
#define HUGEDBL 1e38
#endif



/* Normal CDF */
/*
   Function normbase() implements approximations from
   Cody, W. J. Rational Chebyshev approximations for the error function
     Math. Comp. 23 (1969) 631-637
   Also in W. J. Kennedy and J. E. Gentle, Statistical Computing,
      Marcel Dekker 1980, pp. 90-92
*/

/* coefficients from Cody Table II, n=3*/
#define P1_0 242.66795523053175
#define P1_1 21.979261618294152
#define P1_2 6.9963834886191355
#define P1_3 -.035609843701815385
#define Q1_0 215.05887586986120
#define Q1_1 91.164905404514901
#define Q1_2 15.082797630407787
#define Q1_3 1.0

/* coefficients from Cody Table III, n=7*/
#define P2_0 300.4592610201616005
#define P2_1 451.9189537118729422
#define P2_2 339.3208167343436870
#define P2_3 152.9892850469404039
#define P2_4 43.16222722205673530
#define P2_5 7.211758250883093659
#define P2_6 .5641955174789739711
#define P2_7 -.0000001368648573827167067
#define Q2_0 300.4592609569832933
#define Q2_1 790.9509253278980272
#define Q2_2 931.3540948506096211
#define Q2_3 638.9802644656311665
#define Q2_4 277.5854447439876434
#define Q2_5 77.00015293522947295
#define Q2_6 12.78272731962942351
#define Q2_7 1.0

/* coefficients from Cody Table IV, n=4*/
#define P3_0 -.00299610707703542174
#define P3_1 -.0494730910623250734
#define P3_2 -.226956593539686930
#define P3_3 -.278661308609647788
#define P3_4 -.0223192459734184686
#define Q3_0 .0106209230528467918
#define Q3_1 .191308926107829841
#define Q3_2 1.05167510706793207
#define Q3_3 1.98733201817135256
#define Q3_4 1.0

#define SQRT2 1.414213562373095049
#define SQRTPI 1.772453850905516027

void normbase(double * x, double * phi)
{
	long            sn;
	double          r1, r2, y, y2, erfVal, erfcVal, z;

	y = *x / SQRT2;
	if (y < 0)
	{
		y = -y;
		sn = -1;
	}
	else
	{
		sn = 1;
	}

	y2 = y * y;
	if (y < 0.46875)
	{
		r1 = ((P1_3 * y2 + P1_2) * y2 + P1_1) * y2 + P1_0;
		r2 = ((Q1_3 * y2 + Q1_2) * y2 + Q1_1) * y2 + Q1_0;
		erfVal = y * r1 / r2;
		if (sn == 1)
		{
			*phi = 0.5 + 0.5 * erfVal;
		}
		else
		{
			*phi = 0.5 - 0.5 * erfVal;
		}
	} /*if (y < 0.46875)*/
	else
	{		
		if (y < 4.0)
		{
			r1 = ((((((P2_7 * y + P2_6) * y + P2_5) * y + P2_4) * y + P2_3) * 
				   y + P2_2) * y + P2_1) * y + P2_0;
			r2 = ((((((Q2_7 * y + Q2_6) * y + Q2_5) * y + Q2_4) * y + Q2_3) *
				   y + Q2_2) * y + Q2_1) * y + Q2_0;
			erfcVal = exp(-y2) * r1 / r2;
		} /*if (y < 4.0)*/
		else
		{
			z = y2*y2;
			r1 = (((P3_4 * z + P3_3) * z + P3_2) * z + P3_1) * z + P3_0;
			r2 = (((Q3_4 * z + Q3_3) * z + Q3_2) * z + Q3_1) * z + Q3_0;
			erfcVal = (exp(-y2) / y) * (1.0 / SQRTPI + r1 / (r2 * y2));
		} /*if (y < 4.0){}else{}*/
		if (sn == 1)
		{
			*phi = 1.0 - 0.5 * erfcVal;
		}
		else
		{
			*phi = 0.5 * erfcVal;
		}
	} /*if (y < 0.46875){}else{}*/	
} /*normbase()*/

/* Normal inverse */

#define SPLIT 0.42
#define A0 2.50662823884
#define A1 -18.61500062529
#define A2 41.39119773534
#define A3 -25.44106049637

#define B1 -8.47351093090
#define B2 23.08336743743
#define B3 -21.06224101826
#define B4 3.13082909833

#define C0 -2.78718931138
#define C1 -2.29796479134
#define C2 4.85014127135
#define C3 2.32121276858

#define D1 3.54388924762
#define D2 1.63706781897
/*
  
  Algorithm as 111 Applied Statistics (1977), Vol 26 No 1 Page 121
  Produces normal deviate corresponding to lower tail area of p.
  The hash sums are the sums of the moduli of the coefficients;
  they have no inherent meanings but are incuded for use in
  checking transcriptions.  Functions fabs, log and sqrt are used.  

  Derived from AS111 Fortran version
*/

double squared(double* x, double* expectation, double** variance,int nvar)
{
	int i,j;
	
    double stat=0.0;
	for(i=1;i<=nvar;i++){
	  stat+=(x[i]-expectation[i])*(x[i]-expectation[i]);	
	}
    return stat;

}

double hc(double* x, double* expectation, double** variance,int nvar)
{
	  int i,j,k;
	  double* q=d_createVector(nvar);
	  int nnvar=0;
	  double xt,phi1,phi2;
	  double stat_temp;
	  double temp;
	  for(i=1;i<=nvar;i++)
	  {
		  if(variance[i][i]>0){
			  nnvar++;
			  xt=fabs(x[i]-expectation[i])/sqrt(variance[i][i]);
			  normbase(&xt,&phi1);
			  xt=-fabs(x[i]-expectation[i])/sqrt(variance[i][i]);
			  normbase(&xt,&phi2);
			  q[nnvar]=fmin(1-phi1+phi2,0.9999);
		  }
	  }
	  for (k = 1; k <= nnvar-1; k++){       
			for (j = 1; j <= nnvar-k; j++)  {
					if (q[j] > q[j+1]) {
						temp = q[j+1]; 
						q[j+1] = q[j]; 
						q[j] = temp; 
			   }
		    }			  
	  }
	  stat_temp=((double)(1/((double)nnvar))-q[1])/sqrt(q[1]*(1-q[1]));
	  for (k = 1; k <= (int)nnvar/2; k++){ 
		   temp=((double)((double)k/((double)nnvar))-q[k])/sqrt(q[k]*(1-q[k]));
		   if(temp>=stat_temp) stat_temp=temp;
		   //cout << q[k]<<endl;
	  }
	  delete[] q;
	  return stat_temp;
}
double rvtdt_brv(int** haps, int** O, int** P, int nfam, int nvar)
{
	int i,j; //compare with RV-TDT BRV, He et al. AJHG 2014
	double ctr_minor=0.0;
	double ctr_major=0.0;
	for(i=1;i<=nfam;i++){
		for(j=1;j<=nvar;j++){
			if(haps[O[i][1]][j]==1){ctr_minor+=1.0;}
			if(haps[O[i][2]][j]==1){ctr_minor+=1.0;}
			if(haps[O[i][1]][j]==0){ctr_major+=1.0;}
			if(haps[O[i][2]][j]==0){ctr_major+=1.0;}
		}
	}
    return (ctr_minor-ctr_major)*(ctr_minor-ctr_major)/(ctr_minor+ctr_major);

}
void gtdt(int** haps, int** O, int** P, int nfam, int nvar, double* s_ad, double* s_dom, double* s_ch)
{
    int i,j,k;
	double stat_ad=0.0;
	double stat_dom=0.0;
	double stat_ch=0.0;
	double tmp_sum1,tmp_sum2,tmp_sum3,tmp_sum4,tmp_sum5,tmp_sum6;
	for(i=1;i<=nfam;i++){
		tmp_sum1=tmp_sum2=tmp_sum3=tmp_sum4=tmp_sum5=tmp_sum6=0.0;
		
		for(j=1;j<=nvar;j++){
		  
		  
		  tmp_sum1+=(double)(haps[O[i][1]][j]);
		  tmp_sum2+=(double)(haps[O[i][2]][j]);
		  tmp_sum3+=(double)(haps[P[i][1]][j]);
		  tmp_sum4+=(double)(haps[P[i][2]][j]);
		  tmp_sum5+=(double)(haps[P[i][3]][j]);
		  tmp_sum6+=(double)(haps[P[i][4]][j]);
		}
		//compare with Chen et al. Bioinformatics 2015
		stat_ad+=tmp_sum1+tmp_sum2-0.5*(tmp_sum3+tmp_sum4+tmp_sum5+tmp_sum6);
		if((tmp_sum1 + tmp_sum2)>0.0){ stat_dom+=1.0;}
		if((tmp_sum3 + tmp_sum5)>0.0){ stat_dom-=0.25;}
		if((tmp_sum3 + tmp_sum6)>0.0){ stat_dom-=0.25;}
		if((tmp_sum4 + tmp_sum5)>0.0){ stat_dom-=0.25;}
		if((tmp_sum4 + tmp_sum6)>0.0){ stat_dom-=0.25;}
		
		if(tmp_sum1>=1.0 && tmp_sum2>=1.0){ stat_ch+=1.0;}
		if(tmp_sum3>=1.0 && tmp_sum5>=1.0){ stat_ch-=0.25;}
		if(tmp_sum3>=1.0 && tmp_sum6>=1.0){ stat_ch-=0.25;}
		if(tmp_sum4>=1.0 && tmp_sum5>=1.0){ stat_ch-=0.25;}
		if(tmp_sum4>=1.0 && tmp_sum6>=1.0){ stat_ch-=0.25;}
		
	}
	*s_ad=stat_ad;
	*s_dom=stat_dom;
	*s_ch=stat_ch;

}
double max_sv(double* x, double* expectation, double** variance,int nvar)
{
	 int i;
	 double stat=0.0;
	 double tmp;
	 for(i=1;i<=nvar;i++)
	 {
		 tmp=sv_fbat(x,expectation,variance,nvar,i);
		 if(tmp*tmp>=stat){stat=tmp*tmp;} 
	 }
	 return stat;
	
}
double sv_fbat(double* x, double* expectation, double** variance,int nvar, int index)
{
    double stat;
    if(variance[index][index]>0.0){ stat=(x[index]-expectation[index])/sqrt(variance[index][index]);}
	else { stat=-0.0;}
    return stat;
}
double burden(double* x, double* expectation, double** variance,int nvar)
{
	int i,j;
	double var=0.0;
	for(i=1;i<=nvar;i++){
		for(j=1;j<=nvar;j++){
			var+=variance[i][j];
		}
	}
    double stat=0.0;
	for(i=1;i<=nvar;i++){
	  stat+=(x[i]-expectation[i]);	
	}
    return stat/sqrt(var);

}

double*	d_createVector(int size)
{
    double *Array=0;
    
	Array=new double[size+1];
	
	
	return Array;
}
void	d_destroyVector(double* vector){
	delete[] vector;
}

int*	i_createVector(int size)
{
	int* Vector;
	Vector=new int[size+1];
	return Vector;
}
void	i_destroyVector(int* vector)
{
	delete[] vector;
}

double**	d_createMatrix(int sizeX, int sizeY){
	int i;
	double **Array=new double*[sizeX+1];
	for(i=0; i<=sizeX; i++){
	        Array[i]=new double[sizeY+1];
	}
	return Array;
}



void	d_destroyMatrix(double** matrix, int rows){
	    int i;
		for(i=0; i<=rows; i++){
			       delete[] matrix[i];
		}
		delete[] matrix;
}
int**	i_createMatrix(int sizeX, int sizeY){
	int i;
	int **Array=new int*[sizeX+1];
	for(i=0; i<=sizeX; i++){
	        Array[i]=new int[sizeY+1];
	}
	return Array;
}



void	i_destroyMatrix(int** matrix, int rows){
	    int i;
		for(i=0; i<=rows; i++){
			       delete[] matrix[i];
		}
		delete[] matrix;
}
