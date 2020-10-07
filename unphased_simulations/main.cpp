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
#include "pedigree.hpp"
#include "data.hpp"
using namespace std;

void read_in(family* f, MersenneTwister* r,const char* filename);
void get_hap(int** haps,int nrow,int ncol,const char* file_haps);

int main(int argc,char *argv[])
{
   int i,j,k;
   // argv[1]: output folder output/ 
   // argv[2]: scenario name
   // argv[3]: seed
   // argv[4]: config_file of scenario
   // argv[5]: haplotype_file
   // argv[6]: number of variants
   // argv[7]: number of families
   // argv[8]: number of replications
   // argv[9]: scenario name without number of variants (for RV-GDT)
   // argv[10]: number of haplotypes
   
   
   MersenneTwister mt(atoi(argv[3]));
   
   
   int nreps=atoi(argv[8]);
   int nfams=atoi(argv[7]);
   int nvar=atoi(argv[6]);
   int nhaps=atoi(argv[10]);
   family f;
   int** haps=i_createMatrix(nhaps,nvar);
   
  
   
   
   
  
   get_hap(haps,nhaps,nvar,argv[5]);
   
   read_in(&f,&mt,argv[4]);
   f.haps=haps;
   int** genomat=i_createMatrix(nfams*f.ntotal,nreps*nvar);
   
   // gTDT
   std::ostringstream gtdt_ped;
   gtdt_ped <<argv[1] << "GTDT/gtdt_"<<argv[2] <<".ped";
   std::ostringstream gtdt_vcf;
   gtdt_vcf <<argv[1] << "GTDT/gtdt_"<<argv[2] <<".vcf";
   std::ostringstream gtdt_grp;
   gtdt_grp <<argv[1] << "GTDT/gtdt_"<<argv[2] <<".grp";
   
   FILE* fp_gtdt_ped;
   fp_gtdt_ped=fopen((gtdt_ped.str()).c_str(),"w");
   FILE* fp_gtdt_vcf;
   fp_gtdt_vcf=fopen((gtdt_vcf.str()).c_str(),"w");
   FILE* fp_gtdt_grp;
   fp_gtdt_grp=fopen((gtdt_grp.str()).c_str(),"w");
   
   // FBAT
   std::ostringstream fbat_ped;
   fbat_ped <<argv[1] << "FBAT/fbat_"<<argv[2] <<".ped";
   FILE* fp_fbat;
   fp_fbat=fopen((fbat_ped.str()).c_str(),"w");
   
   // RV_GDT
   std::ostringstream rvgdt_ped;
   rvgdt_ped <<argv[1] << "RVGDT/"<<argv[9]<< "/rvgdt_"<<argv[2] <<".ped";
   FILE* fp_rvgdt_ped;
   fp_rvgdt_ped=fopen((rvgdt_ped.str()).c_str(),"w");
   
   
   
   
   
   int ctr=0;
   int ctr_tot=0;
   
   for(i=1;i<=nreps;i++)
   {
           ctr=0;
		   ctr_tot=0;
		   while(ctr<nfams)
		   {
			   if(f.draw()==1)
			   {
					
					for(j=1;j<=f.F;j++)
					{
					   for(k=1;k<=nvar;k++)
					   {
						 genomat[ctr_tot+j][(i-1)*nvar+k]=f.members[j].x[1][k]+f.members[j].x[2][k];
						 //cout << genomat[ctr_tot+j][(i-1)*nvar+k]<< " ";
					   }
					   //cout << endl;
					   
					}
					for(j=1;j<=f.O;j++)
					{
					   for(k=1;k<=nvar;k++)
					   {
						 genomat[ctr_tot+f.F+j][(i-1)*nvar+k]=f.members[f.F+j].x[1][k]+f.members[f.F+j].x[2][k];
						 //cout << genomat[ctr_tot+f.F+j][(i-1)*nvar+k]<<" ";
					   }
					   //cout <<endl;
					}
					ctr++;
					ctr_tot=ctr_tot+f.ntotal;
					
			   }
		   }
		   cout<< "rep "<< i<<endl;
   }
   

   //FBAT
   write_fbat_header(fp_fbat, nreps*nvar);
   for(i=1;i<=nfams;i++)
   {
		write_fbat_ped(&f,fp_fbat, genomat,i,nvar*nreps); 
   }

   //GTDT
   for(i=1;i<=nfams;i++)
   {
     write_gtdt_ped(&f,fp_gtdt_ped,i);
   }
   write_gtdt_vcf(&f,fp_gtdt_vcf,genomat,nfams,nreps*nvar);
   write_gtdt_grp(fp_gtdt_grp,nreps,nvar);
   
   //RV-GDT
   for(i=1;i<=nfams;i++)
   {
	   write_rvgdt_ped(&f,fp_rvgdt_ped,i); 
   }
   for(i=1;i<=nreps;i++)
   {
	   write_rvgdt_geno(&f,argv,i,genomat,nvar,nfams);
   }
}

void read_in(family* f, MersenneTwister* r,const char* filename)
{
    FILE* fp;
    
    fp=fopen(filename,"r");
   
    int i,j;
	int tmp_np,tmp_op,tmp_int,states, nvar, ncausal;
	double tmp_double;
    fscanf(fp,"%i",&tmp_np);
	
    fscanf(fp,"\n");
	fscanf(fp,"%i",&tmp_op);
	
  
    fscanf(fp,"\n");
	fscanf(fp,"%i",&states);
	
	fscanf(fp,"\n");
	fscanf(fp,"%i",&nvar);
	
    fscanf(fp,"\n");
	
	fscanf(fp,"%i",&ncausal);
	
    fscanf(fp,"\n");
	
	f->initialize(tmp_np,tmp_op,r,states,nvar,ncausal);
	
	///-----------------
	int* a=i_createVector(f->ntotal);
	for(i=1;i<=f->ntotal;i++)
	{
	  fscanf(fp,"%i",&tmp_int);
	  a[i]=tmp_int;
      //cout << a[i]<<endl;
	}
	fscanf(fp,"\n");
	int** pp=i_createMatrix(f->O,3);
	
    for(i=1;i<=f->O;i++)
	{
	    fscanf(fp,"%i",&tmp_int);
        pp[i][1]=tmp_int;
        fscanf(fp,"%i",&tmp_int);
        pp[i][2]=tmp_int;
        fscanf(fp,"%i",&tmp_int);
        pp[i][3]=tmp_int;	
        //cout << pp[i][1]<< "  "<<pp[i][2]<<"  "<<pp[i][3]<<endl;		
	}
	fscanf(fp,"\n");
	
	
	
	double** pf=d_createMatrix(f->F,states);
	
	for(i=1;i<=f->F;i++)
	{
	   for(j=1;j<=states;j++)
	   {
	    fscanf(fp,"%le",&tmp_double);
        pf[i][j]=tmp_double;

	   }
	   fscanf(fp,"\n");
	}
    
	int* ind_causal=i_createVector(ncausal);
    double* beta_causal=d_createVector(ncausal);
	
	for(i=1;i<=ncausal;i++)
	{
		fscanf(fp,"%i",&tmp_int);
		ind_causal[i]=tmp_int;
		
	}
	fscanf(fp,"\n");
	for(i=1;i<=ncausal;i++)
	{
		fscanf(fp,"%le",&tmp_double);
		beta_causal[i]=tmp_double;
	}
	fscanf(fp,"\n");
	int* obs=i_createVector(f->ntotal);
	
	for(i=1;i<=f->ntotal;i++)
	{
	    fscanf(fp,"%i",&tmp_int);
        obs[i]=tmp_int;
		//cout << obs[i]<< " ";
	}
	//cout <<endl;
	fscanf(fp,"\n");
	int* gen=i_createVector(f->ntotal);
	
	for(i=1;i<=f->ntotal;i++)
	{
	    fscanf(fp,"%i",&tmp_int);
        gen[i]=tmp_int;
		//cout << gen[i]<< " ";
	}
	//cout <<endl;
	fscanf(fp,"\n");
    f->set_fam(a,pp,pf,ind_causal,beta_causal,obs,gen);
    
	fclose(fp);	
	i_destroyVector(a); 
	i_destroyVector(ind_causal); 
	d_destroyVector(beta_causal);
	d_destroyMatrix(pf,f->F);
	i_destroyMatrix(pp,f->O);
	i_destroyVector(obs);
	i_destroyVector(gen);
	
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





