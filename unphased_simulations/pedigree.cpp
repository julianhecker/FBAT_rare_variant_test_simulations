#include <iostream>

#include <algorithm>
#include <time.h>
#include <vector>
#include <math.h>
#include <string.h>



#include <stdio.h>
#include <stdlib.h>




#include "pedigree.hpp"
#include <sys/time.h>
using namespace std;

void family::initialize(int np, int op,MersenneTwister* rand,int statesp,int nvarp, int ncausalp) //x
{
	   int i;
	   F=np;
	   O=op;
	   ntotal=F+O;
	   r=rand;
	   ncausal=ncausalp;
	   
	   founder_states=statesp;
	   
	   nmarker=nvarp;
	   founder_distr=d_createMatrix(F,founder_states);
	   
	   observed=i_createVector(ntotal);
	   gender=i_createVector(ntotal);
	   affection=i_createVector(ntotal);
	   po=i_createMatrix(ntotal,3);
	   
	   ind_csl=i_createVector(ncausal);
	   beta_csl=d_createVector(ncausal);
	   //cout<< nmarker<<endl;
	   members=new member[ntotal+1];
	   for(i=1;i<=ntotal;i++){ members[i].nmarker=nmarker;members[i].init(); members[i].r=rand; }
}

void member::inherit(int** p1, int** p2, int ind1, int ind2) 
{
     int i; // DSL index 1, marker index 2 to nmarker+1
	 
	 
	 
	for(i=1;i<=nmarker;i++)
	{
		  x[1][i]=p1[ind1][i];
	}
	 
	
	 
	for(i=1;i<=nmarker;i++)
	{
		  x[2][i]=p2[ind2][i];
	} 
}

void member::init() //x
{
     
     x=i_createMatrix(2,nmarker);
     
}
void member::draw_founders(double* probs,int** haps)
{
        int i;
        double alpha=2.0*r->randDouble();
		int index_ctr=1;
		double cum_prob=0.0;
		while(cum_prob<alpha)
		{
			cum_prob+=probs[index_ctr];
			index_ctr++;
			
		}
		index_ctr=index_ctr-1;
		
		for(i=1;i<=nmarker;i++)
		{
		   x[1][i]=haps[index_ctr][i];
		   
		}
		
        alpha=2.0*r->randDouble();
		index_ctr=1;
		cum_prob=0.0;
		while(cum_prob<alpha)
		{
			cum_prob+=probs[index_ctr];
			index_ctr++;
		}
		index_ctr=index_ctr-1;
		for(i=1;i<=nmarker;i++)
		{
		   x[2][i]=haps[index_ctr][i];
		}		
		
}

int family::draw_d() //x
{
	    int i,j;
		double prob;
		int geno;
		double alpha;
		int aff;
		int check=1;
		for(i=1;i<=ntotal;i++)
		{
		   
		   prob=0.0;
		   for(j=1;j<=ncausal;j++)
		   {
			    geno=members[i].x[1][ind_csl[j]]+members[i].x[2][ind_csl[j]]; //DSL geno
			    prob+=geno*beta_csl[j];
		   }
		   prob=0.1*exp(prob);
		
		   alpha=2.0*r->randDouble();
		   aff=0;
		   if(alpha<=prob){aff=1;}
		   if(affection[i]==1 && aff==0){ check=0; }
		   if(affection[i]==0 && aff==1){ check=0; }
		}
		
		return check;
}


void family::set_fam(int* aff,int** pop,double** f_d,int* ind_cslp,double* beta_cslp,int* obs,int* gen) //x
{
			int i,j;
		    
			
			
			for(i=1;i<=ntotal;i++)
			{
			  affection[i]=aff[i];
			  
			}
			for(i=1;i<=O;i++)
			{
			  
			  for(j=1;j<=3;j++)
			  {
				 po[i][j]=pop[i][j];
				 
			  }
			}
			for(i=1;i<=F;i++)
			{
			  
			  for(j=1;j<=founder_states;j++)
			  {
				 founder_distr[i][j]=f_d[i][j];
				 
			  }
			}
			for(i=1;i<=ncausal;i++)
			{
			  ind_csl[i]=ind_cslp[i];
			  beta_csl[i]=beta_cslp[i];
			  
			}
			for(i=1;i<=ntotal;i++)
			{
			  observed[i]=obs[i];
			  
			  
			}
		
			for(i=1;i<=ntotal;i++)
			{
			  gender[i]=gen[i];
			  
			}
			
			
}

int family::draw() //x
{
           int ret;
	       int i,j;
		   for(i=1;i<=F;i++)
		   {
		      
			  members[i].draw_founders(founder_distr[i],haps);
			 
			  
		   }
		   draw_v();
		   ret=draw_d();
		   return ret;
}
void family::draw_v() //x
{
	       int i;
		   
	       double alpha;
		   int parent1,parent2,offspring;
		   int ind1,ind2;
		   for(i=1;i<=O;i++)
		   {
		      parent1=po[i][1];
			  parent2=po[i][2];
			  offspring=po[i][3];
			  alpha=2.0*r->randDouble();
			  ind1=2;
			  if(alpha<=0.5){ ind1=1;}
			  alpha=2.0*r->randDouble();
			  ind2=2;
			  if(alpha<=0.5){ ind2=1;}
			  members[offspring].inherit(members[parent1].x,members[parent2].x,ind1,ind2);
			  
			  
		   }
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