//
//  sweepdefs.h
//  
//
//  Created by Lindsay Martin on 10/22/18.
//

#ifndef sweepdefs_h
#define sweepdefs_h

#define DIM1 2
#define DIM2 800
#include <iostream>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <stdio.h>
using namespace std;
#include <fstream>

int fsm(double *u, double *f, double *BC, double *obs, int iter, int N, double h);

void coarseinitial(double *Ucoarse[DIM1][DIM2],double * R, double *BC,double *obs, int N,int K, int p, int z, double *wind[DIM1][DIM2], int sweepiter,int *totalsweeps);


void coarseupdate(double *Ucoarse[DIM1][DIM2], double *Ufine, double * UK[DIM1][DIM2], double  * R, double *BC,double *obs, int N,int K,int p, int z, double *wind[DIM1][DIM2], double *fwind, double T, int sweepiter,int *totalsweeps);

void coarseupdatenewtheta(double *Ucoarse[DIM1][DIM2], double *Ucoarse1[DIM1][DIM2], double *Ufine1, double *Ufine2,  double *UK1[DIM1][DIM2], double *UK2[DIM1][DIM2], double *UK3[DIM1][DIM2],  double * R, double *BC,double *obs, int N,int K, int p, int z, double *wind[DIM1][DIM2], double *fwind, int sweepiter, double epsilon, double delta, double x0, double *theta[DIM1][DIM2], int wt1, int wt2, int wt3,int *totalsweeps);

void causalsweeping(double *Ucoarse[DIM1][DIM2], double *BC,double *obs, int N,int K,double *wind[DIM1][DIM2],double* R);

void fineupdate(double* FINE, double *Ucoarse[DIM1][DIM2], double *f, double *BC,double *obs, int N,int K, double *wind[DIM1][DIM2], double *fwind,int *totalsweeps);




#endif /* sweepdefs_h */
