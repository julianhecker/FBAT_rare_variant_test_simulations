#include "data.hpp"

#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

void write_gtdt_vcf(family* f, FILE* fp,int** genomat,int nfams,int nvar)
{
   int i,j,k;
   fprintf(fp,"##fileformat=VCFv4.0\n");
   fprintf(fp,"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
   for(j=1;j<=nfams;j++)
   {
	   for(i=1;i<=f->F;i++)
	   {
		  if(f->observed[i]==1)
		  {
			 fprintf(fp,"\t%i:%i",j,i);
		  }
		 
	   }
	   for(i=1;i<=f->O;i++)
	   {
			if(f->observed[i+f->F]==1)
			{
			 
				  fprintf(fp,"\t%i:%i",j,f->po[i][3]);	 
			}   
	   }
   }
   fprintf(fp,"\n");
   for(j=1;j<=nvar;j++)
   {
	   fprintf(fp,"1\t%i\tm%i\tA\tG\t30\tPASS\tDP=20\tGT",j,j);
	   for(k=1;k<=nfams;k++)
	   {
			   for(i=1;i<=f->F;i++)
			   {
				  if(f->observed[i]==1)
				  {
					   if(genomat[(k-1)*f->ntotal+i][j]==2){ fprintf(fp,"\t1/1");}
					   if(genomat[(k-1)*f->ntotal+i][j]==1){ fprintf(fp,"\t0/1");}
					   if(genomat[(k-1)*f->ntotal+i][j]==0){ fprintf(fp,"\t0/0");}
				  }
				 
			   }
			   for(i=1;i<=f->O;i++)
			   {
					if(f->observed[i+f->F]==1)
					{
					   if(genomat[(k-1)*f->ntotal+f->F+i][j]==2){ fprintf(fp,"\t1/1");}
					   if(genomat[(k-1)*f->ntotal+f->F+i][j]==1){ fprintf(fp,"\t0/1");}
					   if(genomat[(k-1)*f->ntotal+f->F+i][j]==0){ fprintf(fp,"\t0/0");}
					}   
			   }
	   }
	   fprintf(fp,"\n");
   }
}
void write_gtdt_ped(family* f, FILE* fp, int id)
{
   int i,j;
   //cout << id<<"  "<<f->O<<"  "<<f->F<<"  " << f->observed[1]<<"  "<< f->observed[2]<<"  "<< f->observed[3]<<"  "<<endl;
   for(i=1;i<=f->F;i++)
   {
      if(f->observed[i]==1)
	  {
		 fprintf(fp,"%i %i:%i 0 0 %i %i",id,id,i,f->gender[i],f->affection[i]+1);
		 fprintf(fp,"\n");
	  }
     
   }
   for(i=1;i<=f->O;i++)
   {
		if(f->observed[i+f->F]==1)
		{
		 
			  fprintf(fp,"%i %i:%i %i:%i %i:%i %i %i",id,id,f->po[i][3],id,f->po[i][1],id,f->po[i][2],f->gender[i+f->F],f->affection[i+f->F]+1);
			  
			  fprintf(fp,"\n");
		}
      
   }
}
void write_gtdt_grp(FILE* fp, int nreps,int nvar)
{
   int i,j;
   fprintf(fp,"#CHR\tPOS\tGroup\n");
   for(i=1;i<=nreps;i++)
   {
      for(j=1;j<=nvar;j++)
	  {
		 fprintf(fp,"1\t%i\tGene%i\n",(i-1)*nvar+j,i);
		 
	  }
     
   }
  
}


void write_rvgdt_ped(family* f, FILE* fp, int id)
{
   int i,j;
   for(i=1;i<=f->F;i++)
   {
      if(f->observed[i]==1)
	  {
		 fprintf(fp,"%i %i:%i 0 0 %i %i",id,id,i,f->gender[i],f->affection[i]+1);
		 fprintf(fp,"\n");
	  }
     
   }
   for(i=1;i<=f->O;i++)
   {
		if(f->observed[i+f->F]==1)
		{
		 
			  fprintf(fp,"%i %i:%i %i:%i %i:%i %i %i",id,id,f->po[i][3],id,f->po[i][1],id,f->po[i][2],f->gender[i+f->F],f->affection[i+f->F]+1);
			  
			  fprintf(fp,"\n");
		}
      
   }
}
void write_rvgdt_geno(family* f, char *argv[], int id,int** geno_mat, int nvars, int nfams)
{
   int i,j,k;
   
   
   std::ostringstream rvgdt_geno;
   rvgdt_geno <<argv[1] <<"RVGDT/"<<argv[9]<<"/rvgdt_"<<argv[2] <<"_"<<id<<".geno";
   FILE* fp;
   fp=fopen((rvgdt_geno.str()).c_str(),"w");
   
   for(k=1;k<=nfams;k++)
   {
		   for(i=1;i<=f->F;i++)
		   {
			  if(f->observed[i]==1)
			  {
				 fprintf(fp,"%i:%i",k,i);
				 for(j=1;j<=nvars;j++)
				 {
				   if(geno_mat[(k-1)*f->ntotal+i][(id-1)*nvars+j]==2){ fprintf(fp," 2");}
				   if(geno_mat[(k-1)*f->ntotal+i][(id-1)*nvars+j]==1){ fprintf(fp," 1");}
				   if(geno_mat[(k-1)*f->ntotal+i][(id-1)*nvars+j]==0){ fprintf(fp," 0");}
				 }
				 fprintf(fp,"\n");
			  }
			 
		   }
		   for(i=1;i<=f->O;i++)
		   {
				if(f->observed[i+f->F]==1)
				{
				 
					  fprintf(fp,"%i:%i",k,f->po[i][3]);
					  for(j=1;j<=nvars;j++)
					  {
						if(geno_mat[(k-1)*f->ntotal+i+f->F][(id-1)*nvars+j]==2){ fprintf(fp," 2");}
						if(geno_mat[(k-1)*f->ntotal+i+f->F][(id-1)*nvars+j]==1){ fprintf(fp," 1");}
						if(geno_mat[(k-1)*f->ntotal+i+f->F][(id-1)*nvars+j]==0){ fprintf(fp," 0");}
					  }
					  fprintf(fp,"\n");
				}
			  
		   }
   }
   fclose(fp);
}
void write_gdt_dat(FILE* fp, int nvars)
{
   int i;
   fprintf(fp,"A trait",i);
   fprintf(fp,"\n");
   for(i=1;i<=nvars;i++)
   {
      fprintf(fp,"M m%i",i);
	  fprintf(fp,"\n");
   }
   
}
void write_gdt_map(FILE* fp, int nvars)
{
   int i;
   for(i=1;i<=nvars;i++)
   {
      fprintf(fp,"24 m%i 0.000",i);
	  fprintf(fp,"\n");
   }
   
}
void write_gdt_ped(FILE* fp, int id, int nvars,family* f, int** geno_mat)
{
       int i,j;
	   for(i=1;i<=f->F;i++)
	   {
		  if(f->observed[i]==1)
		  {
			 fprintf(fp,"%i %i 0 0 %i %i",id,i,f->gender[i],f->affection[i]+1);
			 for(j=1;j<=nvars;j++)
			 {
			   if(geno_mat[(id-1)*f->ntotal+i][j]==2){ fprintf(fp," 2 2");}
			   if(geno_mat[(id-1)*f->ntotal+i][j]==1){ fprintf(fp," 1 2");}
			   if(geno_mat[(id-1)*f->ntotal+i][j]==0){ fprintf(fp," 1 1");}
			 }
			 fprintf(fp,"\n");
		  }
		  if(f->observed[i]==0)
		  {
			 fprintf(fp,"%i %i 0 0 %i %i",id,i,f->gender[i],f->affection[i]+1);
			 for(j=1;j<=nvars;j++)
			 {
			   fprintf(fp," 0 0");
			 }
			 fprintf(fp,"\n");
		  }
		 //cout <<id<<" "<<i<<" 0 0 "<< f->gender[i]<<" "<<f->affection[i]<<endl;
	   }
	   for(i=1;i<=f->O;i++)
	   {
			if(f->observed[i+f->F]==1)
			{
			 
				  fprintf(fp,"%i %i %i %i %i %i",id,f->po[i][3],f->po[i][1],f->po[i][2],f->gender[i+f->F],f->affection[i+f->F]+1);
				  for(j=1;j<=nvars;j++)
				  {
					if(geno_mat[(id-1)*f->ntotal+i+f->F][j]==2){ fprintf(fp," 2 2");}
					if(geno_mat[(id-1)*f->ntotal+i+f->F][j]==1){ fprintf(fp," 1 2");}
					if(geno_mat[(id-1)*f->ntotal+i+f->F][j]==0){ fprintf(fp," 1 1");}
				  }
				  fprintf(fp,"\n");
			}
			if(f->observed[i+f->F]==0)
			{
			 
				  fprintf(fp,"%i %i %i %i %i %i",id,f->po[i][3],f->po[i][1],f->po[i][2],f->gender[i+f->F],f->affection[i+f->F]+1);
				  for(j=1;j<=nvars;j++)
				  {
					fprintf(fp," 0 0");
				  }
				  fprintf(fp,"\n");
			}
		  //cout <<id<<" "<<f->po[i][3] <<" "<< f->po[i][1]<<" "<<f->po[i][2]<< " "<< f->gender[i+f->F]<<" "<<f->affection[i+f->F]<<endl;
	   }
}
void write_fbat_header(FILE* fp, int nvars)
{
   int i;
   for(i=1;i<=nvars;i++)
   {
      fprintf(fp,"m%i ",i);
   }
   fprintf(fp,"\n");
}
/*void write_fbat_ped(FILE* fp,FILE* fp_batch,int nreps,int nvars)
{
   fprintf(fp,"load %s ",i,fp_batch);
   fprintf(fp,"log  \n");
   
}*/
void write_fbat_ped(family* f, FILE* fp, int** geno_mat,int id,int nvars)
{
   int i,j;
   for(i=1;i<=f->F;i++)
   {
      if(f->observed[i]==1)
	  {
		 fprintf(fp,"%i %i 0 0 %i %i",id,i,f->gender[i],f->affection[i]+1);
		 for(j=1;j<=nvars;j++)
		 {
		   if(geno_mat[(id-1)*f->ntotal+i][j]==2){ fprintf(fp," 2 2");}
		   if(geno_mat[(id-1)*f->ntotal+i][j]==1){ fprintf(fp," 1 2");}
		   if(geno_mat[(id-1)*f->ntotal+i][j]==0){ fprintf(fp," 1 1");}
		 }
		 fprintf(fp,"\n");
	  }
     
   }
   for(i=1;i<=f->O;i++)
   {
		if(f->observed[i+f->F]==1)
		{
		 
			  fprintf(fp,"%i %i %i %i %i %i",id,f->po[i][3],f->po[i][1],f->po[i][2],f->gender[i+f->F],f->affection[i+f->F]+1);
			  for(j=1;j<=nvars;j++)
			  {
				if(geno_mat[(id-1)*f->ntotal+i+f->F][j]==2){ fprintf(fp," 2 2");}
				if(geno_mat[(id-1)*f->ntotal+i+f->F][j]==1){ fprintf(fp," 1 2");}
				if(geno_mat[(id-1)*f->ntotal+i+f->F][j]==0){ fprintf(fp," 1 1");}
			  }
			  fprintf(fp,"\n");
		}
      
   }
}

