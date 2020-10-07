#include <iostream>
#include <algorithm>
#include <time.h>
#include <vector>
#include <math.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <sys/time.h>
#include "data.hpp"
using namespace std;

void get_hap(int** haps,int nrow,int ncol,const char* file_haps);
void get_sum_minor_alleles(int** haps, int** O, int nfam,int nvar,double* x);
void get_moments(int** haps, int** O, int** P, int nfam,int nvar,double* expectation,double** variance);
int main(int argc,char *argv[])
{
   int i,j,k;
   
   MersenneTwister mt(atoi(argv[1]));
   
   double correct_factor=1.0-1.0*pow(10,-10);
   int nhaps=1006;
   int nvar=1000;
   int nreps=100;
   int nsim=1000;
   int nfam=atoi(argv[2]);
   int ncausal=atoi(argv[3]);
   
   int** haps=i_createMatrix(nhaps,nvar);
   int** P=i_createMatrix(nfam,4);
   int** O=i_createMatrix(nfam,2);
  
   int* causal_inds=i_createVector(ncausal);
   double* beta=d_createVector(ncausal);
   for(i=1;i<=ncausal;i++){
	   causal_inds[i]=atoi(argv[3+i]);
	   beta[i]=atof(argv[3+ncausal+i]);
   }
  
   get_hap(haps,nhaps,nvar,"haps");
   
   double* x=d_createVector(nvar);
   double* expectation=d_createVector(nvar);
   double** variance=d_createMatrix(nvar,nvar);
   
   double tmpsum;
   double prob=0;
   
   double fbat_burden;
   double fbat_squared;
   double fbat_max;
   double fbat_hc;
   double gtdt_ad;
   double gtdt_dom;
   double gtdt_ch;
   double rvtdtbrv;
   double* fbat_svs=d_createVector(ncausal);
   
   double fbat_burden_t;
   double fbat_squared_t;
   double fbat_max_t;
   double fbat_hc_t;
   double gtdt_ad_t;
   double gtdt_dom_t;
   double gtdt_ch_t;
   double rvtdtbrv_t;
   double* fbat_svs_t=d_createVector(ncausal);
  
   int ctr_fbat_burden=0;
   int ctr_fbat_squared=0;
   int ctr_fbat_max=0;
   int ctr_fbat_hc=0;
   int ctr_gtdt_ad=0;
   int ctr_gtdt_dom=0;
   int ctr_gtdt_ch=0;
   int ctr_rvtdtbrv=0;
   int* ctr_fbat_svs=i_createVector(ncausal);
   for(i=1;i<=ncausal;i++){ctr_fbat_svs[i]=0;}

   int ctr_fam;
   int ind1,ind2,ind3,ind4,off1,off2;
   for(i=1;i<=nreps;i++)
   {
	       ctr_fam=0;
		   while(ctr_fam<nfam){
			   
			   ind1=(int)(nhaps*2.0*mt.randDouble())+1;
			   ind2=(int)(nhaps*2.0*mt.randDouble())+1;
			   ind3=(int)(nhaps*2.0*mt.randDouble())+1;
			   ind4=(int)(nhaps*2.0*mt.randDouble())+1;
			   
			   off1=ind2;
			   off2=ind4;
			   if(2.0*mt.randDouble()<=0.5){ off1=ind1;}
			   if(2.0*mt.randDouble()<=0.5){ off2=ind3;}
			   
			   tmpsum=0.0;
			   for(k=1;k<=ncausal;k++){
				   tmpsum+=beta[k]*(double)(haps[off1][causal_inds[k]]+haps[off2][causal_inds[k]]);
			   }
			   prob=0.1*exp(tmpsum);
			   if(2.0*mt.randDouble()<=prob){
				   ctr_fam++;
				   P[ctr_fam][1]=ind1; P[ctr_fam][2]=ind2; P[ctr_fam][3]=ind3; P[ctr_fam][4]=ind4;
				   O[ctr_fam][1]=off1; O[ctr_fam][2]=off2;
			   }
		   }
           
		   get_sum_minor_alleles(haps,O,nfam,nvar,x);
		   get_moments(haps,O,P, nfam,nvar,expectation,variance);
		   
		   
		   ///////////////////////////////////////////////////
		   fbat_burden=burden(x,expectation,variance,nvar);
		   fbat_squared=squared(x,expectation,variance,nvar);
		   fbat_hc=hc(x,expectation,variance,nvar);
		   fbat_max=max_sv(x,expectation,variance,nvar);
		   for(j=1;j<=ncausal;j++)
		   {
			   fbat_svs[j]=sv_fbat(x,expectation,variance,nvar,causal_inds[j]);
		   }
		   gtdt(haps,O,P,nfam,nvar,&gtdt_ad,&gtdt_dom,&gtdt_ch);
		   rvtdtbrv=rvtdt_brv(haps,O,P,nfam,nvar);
		   
		   ////////////////////////////////////////////////////
		   ctr_fbat_burden=0;
		   ctr_fbat_squared=0;
           ctr_fbat_max=0;
           ctr_fbat_hc=0;
		   ctr_gtdt_ad=0;
		   ctr_gtdt_dom=0;
		   ctr_gtdt_ch=0;
		   ctr_rvtdtbrv=0;
   
		   for(k=1;k<=ncausal;k++){ctr_fbat_svs[k]=0;}
		   for(k=1;k<=nsim;k++){
			  for(j=1;j<=nfam;j++){
				  O[j][1]=P[j][2];
			      if(2.0*mt.randDouble()<=0.5){ O[j][1]=P[j][1];}
				  O[j][2]=P[j][4];
			      if(2.0*mt.randDouble()<=0.5){ O[j][2]=P[j][3];}
			  }
			  //cout << k<<endl;
			  get_sum_minor_alleles(haps,O,nfam,nvar,x);
			  
			  fbat_burden_t=burden(x,expectation,variance,nvar);
			  fbat_squared_t=squared(x,expectation,variance,nvar);
		      fbat_hc_t=hc(x,expectation,variance,nvar);
		      fbat_max_t=max_sv(x,expectation,variance,nvar);
		      for(j=1;j<=ncausal;j++)
		      {
			     fbat_svs_t[j]=sv_fbat(x,expectation,variance,nvar,causal_inds[j]);
		      }
		      gtdt(haps,O,P,nfam,nvar,&gtdt_ad_t,&gtdt_dom_t,&gtdt_ch_t);
		      rvtdtbrv_t=rvtdt_brv(haps,O,P,nfam,nvar);
		   
			  if(correct_factor*fabs(fbat_burden)<=fabs(fbat_burden_t)){ ctr_fbat_burden++;}
			  if(correct_factor*fabs(fbat_squared)<=fabs(fbat_squared_t)){ ctr_fbat_squared++;}
			  if(correct_factor*fabs(fbat_max)<=fabs(fbat_max_t)){ ctr_fbat_max++;}
			  if(correct_factor*fabs(fbat_hc)<=fabs(fbat_hc_t)){ ctr_fbat_hc++;}
			  
			  if(correct_factor*fabs(gtdt_ad)<=fabs(gtdt_ad_t)){ ctr_gtdt_ad++;}
			  if(correct_factor*fabs(gtdt_dom)<=fabs(gtdt_dom_t)){ ctr_gtdt_dom++;}
		      if(correct_factor*fabs(gtdt_ch)<=fabs(gtdt_ch_t)){ ctr_gtdt_ch++;}
			  
			  if(correct_factor*fabs(rvtdtbrv)<=fabs(rvtdtbrv_t)){ ctr_rvtdtbrv++;}
			  
			  for(j=1;j<=ncausal;j++){
				  if(fabs(fbat_svs[j])<=fabs(fbat_svs_t[j])){ ctr_fbat_svs[j]++;}
			  }
		   }
		   cout << "Burden: ";
		   cout << (double)ctr_fbat_burden/(double)nsim<<" SKAT: ";
		   cout << (double)ctr_fbat_squared/(double)nsim<<" MAX: ";
		   cout << (double)ctr_fbat_max/(double)nsim<<" HC: ";
		   cout << (double)ctr_fbat_hc/(double)nsim<<" gTDT-AD: ";
		   
		   cout << (double)ctr_gtdt_ad/(double)nsim<<" gTDT-DOM: ";
		   cout << (double)ctr_gtdt_dom/(double)nsim<<" gTDT-CH: ";
		   cout << (double)ctr_gtdt_ch/(double)nsim<<" RV-TDT-BRV: ";
		   
		   cout << (double)ctr_rvtdtbrv/(double)nsim;
		   cout<<" MAX-FBAT-STATISTIC: ";
		   
		   cout << fbat_max<<" ";
		   cout <<endl;
   }
   i_destroyMatrix(haps,nhaps);
   i_destroyMatrix(P,nfam);
   i_destroyMatrix(O,nfam);
   
   i_destroyVector(causal_inds);
   d_destroyVector(beta);
   d_destroyVector(x);
   d_destroyVector(expectation);
   d_destroyMatrix(variance, nvar);
   d_destroyVector(fbat_svs);
   d_destroyVector(fbat_svs_t);
   i_destroyVector(ctr_fbat_svs);
}
void get_moments(int** haps, int** O, int** P, int nfam,int nvar,double* expectation,double** variance)
{
	int i,j,k;
	for(i=1;i<=nvar;i++){
		 expectation[i]=0.0;
		 for(j=1;j<=nfam;j++){
			  expectation[i]+=0.5*(double)(haps[P[j][1]][i]+haps[P[j][2]][i]+haps[P[j][3]][i]+haps[P[j][4]][i]);
	     }	
	}
	double exp_i,exp_j;
	for(i=1;i<=nvar;i++){
		 
		 for(j=i;j<=i;j++){ // just need the diagonal for the current purposes
			     
		      	 variance[i][j]=0.0;
				 for(k=1;k<=nfam;k++){
					 exp_i=0.5*(double)(haps[P[k][1]][i]+haps[P[k][2]][i]+haps[P[k][3]][i]+haps[P[k][4]][i]);
					 exp_j=0.5*(double)(haps[P[k][1]][j]+haps[P[k][2]][j]+haps[P[k][3]][j]+haps[P[k][4]][j]);
					 variance[i][j]+=0.25*((double)haps[P[k][1]][i]+(double)haps[P[k][3]][i]-exp_i)*((double)haps[P[k][1]][j]+(double)haps[P[k][3]][j]-exp_j);
				     variance[i][j]+=0.25*((double)haps[P[k][1]][i]+(double)haps[P[k][4]][i]-exp_i)*((double)haps[P[k][1]][j]+(double)haps[P[k][4]][j]-exp_j);
					 variance[i][j]+=0.25*((double)haps[P[k][2]][i]+(double)haps[P[k][3]][i]-exp_i)*((double)haps[P[k][2]][j]+(double)haps[P[k][3]][j]-exp_j);
					 variance[i][j]+=0.25*((double)haps[P[k][2]][i]+(double)haps[P[k][4]][i]-exp_i)*((double)haps[P[k][2]][j]+(double)haps[P[k][4]][j]-exp_j);
				 }
				 
				 variance[j][i]=variance[i][j];
	     } 
	}
}
void get_sum_minor_alleles(int** haps, int** O, int nfam,int nvar,double* x)
{
	int i,j;
	for(i=1;i<=nvar;i++){
	   x[i]=0.0;
       for(j=1;j<=nfam;j++){
		   x[i]+=(double)haps[O[j][1]][i]+(double)haps[O[j][2]][i];
       }		
	}
}
void get_hap(int** haps,int nrow,int ncol,const char* file_haps)
{
    FILE* fp;
    fp=fopen(file_haps,"r");
    int i,j;
    int f;
    for(i=1;i<=nrow;i++)
    {
        for(j=1;j<=ncol;j++)
		{
			fscanf(fp,"%i",&f);
			haps[i][j]=f;
		}
    }
	fclose(fp);
}



 

