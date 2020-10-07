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
#include "pedigree.hpp"

void write_rvgdt_ped(family* f, FILE* fp, int id);
void write_rvgdt_geno(family* f, char *argv[], int id,int** genomat, int nvar, int nfams);

void write_gdt_dat(FILE* fp, int nvars);
void write_gdt_map(FILE* fp, int nvars);
void write_gdt_ped(FILE* fp, int id, int nvars,family* f, int** geno_mat); // geno_mat f.ntotal x nvars

void write_fbat_header(FILE* fp, int nvars);
void write_fbat_ped(family* f, FILE* fp, int** geno_mat,int id,int nvars); // geno_mat f.ntotal x nvars



void write_gtdt_vcf(family* f, FILE* fp,int** genomat,int nfams,int nvar);
void write_gtdt_ped(family* f, FILE* fp, int id);
void write_gtdt_grp(FILE* fp, int nreps,int nvar);
#endif
