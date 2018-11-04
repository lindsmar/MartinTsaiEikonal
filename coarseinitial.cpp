//
// Created by Lindsay Martin on 10/29/16.
//

//#ifndef EIKONALPROJECT_COARSEINITIAL_H
//#define EIKONALPROJECT_COARSEINITIAL_H
//
//  coarseupdate.h
//  Eikonal Project 3
//
//  Created by Lindsay Martin on 9/29/16.
//  Copyright © 2016 Lindsay Martin. All rights reserved.
//




//#define DIM1 2
//#define DIM2 80
//#include <algorithm>
//#include <cmath>
//using namespace std;



/*void coarseinitial(double *Ucoarse[DIM1][DIM2],double * R, double *BC,double *obs, int N,int K, int p, int z, double *wind[DIM1][DIM2], int sweepiter);*/

#include "sweepdefs.h"

/*inline void cartto2CG(int i, int j, int N, int K,int *DIRkLIN){
    
    int DIR=-2;
    int k=-3;
    int LIN=-4;
    
    DIRkLIN[0]=DIR;
    DIRkLIN[1]=k;
    DIRkLIN[2]=LIN;
    
    if ( i%K==0 || j%K==0)
    {
        if (i%K==0){
            if (j%K==0){
                
                DIR=0;
                k=0;
                LIN=i/K+j/K*N;
            }
            
            else {
                DIR=0;
                k=j%K;
                LIN=i/K+j/K*N;
                
            }
            DIRkLIN[0]=DIR;
            DIRkLIN[1]=k;
            DIRkLIN[2]=LIN;
        }
        else{
            if (j%K==0){
                DIR=1;
                k=i%K;
                LIN=i/K+j/K*N;
                
                DIRkLIN[0]=DIR;
                DIRkLIN[1]=k;
                DIRkLIN[2]=LIN;
            }
            else{
                cout << "Invalid i,j " << endl;
            }
            
            
        }
        if (DIR<0 || DIR>=2 || k<0 || k>=K || LIN<0 || LIN>=N*N){
            cerr<<"Invalid i,j: i="<<i<<", j="<<j<<endl;
        }
    }
    
    
    else{ cerr<<"Invalid i,j: i="<<i<<", j="<<j<<endl;}
    
    return;
}*/

inline void cartto2CG(int i, int j, int N, int K,int *DIRkLIN){

    if (i%K==0){
        DIRkLIN[0]=0;
        DIRkLIN[1]=j%K;
        DIRkLIN[2]=i/K+j/K*N;
    }
    else{
        DIRkLIN[0]=1;
        DIRkLIN[1]=i%K;
        DIRkLIN[2]=i/K+j/K*N;
    }
    
    return;
}


void coarseinitial(double *Ucoarse[DIM1][DIM2], double *R, double *BC,double *obs, int N,int K, int p, int z, double *wind[DIM1][DIM2], int sweepiter, int *totalsweeps)
{
    int M;
    double H,h,old,m1,m2,discr,hx,hy,unew,F;
    int tempwindx, tempwindy,istart1,istart2,jstart1,jstart2;;
    double currError;

    H=(double)1/(double) N;
    h=(double)1/(double)(N*K);
    M=N*K+1;
    
    int s1,s2,i,j;

    if (p==0){
        istart1=M-1;
        istart2=0;
    }
    else{
        istart1=M-1-(K-p);
        istart2=p;
    }

    if (z==0){
        jstart1=M-1;
        jstart2=0;
    }
    else{
        jstart1=M-1-(K-z);
        jstart2=z;
    }
  
    double carttime=0;
    double tic,toc;
    
    
    int DIRkLIN[3];
    int dir,k,lin;
    int DIRkLINnbx1[3],DIRkLINnbx2[3];
    int dirnbx1,knbx1,linnbx1,dirnbx2,knbx2,linnbx2;
    int DIRkLINnby1[3],DIRkLINnby2[3];
    int dirnby1,knby1,linnby1,dirnby2,knby2,linnby2;
    for (int iter=0; iter<sweepiter; iter++) {
        double maxError=0;
        for (s1 = 1; s1 >= -1; s1 -= 2){
            for (s2 = 1; s2 >= -1; s2 -= 2){
                for (i = (s1 < 0 ? istart1 : istart2); (s1 < 0 ? i >= 0 : i <= M - 1); i += K * s1){
                    for (j = (s2 < 0 ? jstart1 : jstart2); (s2 < 0 ? j >= 0 : j <= M - 1); j += K * s2) {
                        if (obs[i + j * M] > 1) {
                            //don't do anything; point is in obstacle
                        } else if (BC[i + j * M] >= 0) {
                            //tic=omp_get_wtime();
                            cartto2CG(i,j,N+1,K,DIRkLIN);
                            //toc=omp_get_wtime()-tic;
                            //toc;
                            dir=DIRkLIN[0];
                            k=DIRkLIN[1];
                            lin=DIRkLIN[2];
                            Ucoarse[dir][k][lin] = BC[i + j * M];//Ucoarse[i + j * M] = BC[i + j * M];
                            
                            
                        } else if (BC[i + j * M] < 0) {
                           // tic=omp_get_wtime();
                            cartto2CG(i,j,N+1,K,DIRkLIN);
                            //toc=omp_get_wtime()-tic;
                            //carttime+=toc;
                            dir=DIRkLIN[0];
                            k=DIRkLIN[1];
                            lin=DIRkLIN[2];
                            
                           
                            
                          tempwindx = 0;
                          tempwindy = 0;
                          if (j == z) {
                              //tic=omp_get_wtime();
                              cartto2CG(i,j+K,N+1,K,DIRkLINnby2);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              dirnby2=DIRkLINnby2[0];
                              knby2=DIRkLINnby2[1];
                              linnby2=DIRkLINnby2[2];
                              m1 = Ucoarse[dirnby2][knby2][linnby2];
                              tempwindy = 7;
                              hx = H;
                          } else if (j == M - 1 - (K - z) && z != 0) {
                              // tic=omp_get_wtime();
                              cartto2CG(i,j-K,N+1,K,DIRkLINnby1);
                              //toc=omp_get_wtime()-tic;
                             // carttime+=toc;
                              dirnby1=DIRkLINnby1[0];
                              knby1=DIRkLINnby1[1];
                              linnby1=DIRkLINnby1[2];
                              
                              m1 =Ucoarse[dirnby1][knby1][linnby1];
                              hx = H;
                              tempwindy = 5;
                          } else if (j == M - 1) {
                             // tic=omp_get_wtime();
                              cartto2CG(i,j-K,N+1,K,DIRkLINnby1);
                             // toc=omp_get_wtime()-tic;
                             // carttime+=toc;
                              dirnby1=DIRkLINnby1[0];
                              knby1=DIRkLINnby1[1];
                              linnby1=DIRkLINnby1[2];
                              
                              m1 =Ucoarse[dirnby1][knby1][linnby1];
                              
                             
                              hx = H;
                              tempwindy = 5;
                          } else {
                             // tic=omp_get_wtime();
                              cartto2CG(i,j+K,N+1,K,DIRkLINnby2);
                             // toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              dirnby2=DIRkLINnby2[0];
                              knby2=DIRkLINnby2[1];
                              linnby2=DIRkLINnby2[2];
                              // tic=omp_get_wtime();
                              cartto2CG(i,j-K,N+1,K,DIRkLINnby1);
                              //toc=omp_get_wtime()-tic;
                             // carttime+=toc;
                              dirnby1=DIRkLINnby1[0];
                              knby1=DIRkLINnby1[1];
                              linnby1=DIRkLINnby1[2];
                              m1 = min(Ucoarse[dirnby1][knby1][linnby1], Ucoarse[dirnby2][knby2][linnby2]);
                              if (Ucoarse[dirnby1][knby1][linnby1] <Ucoarse[dirnby2][knby2][linnby2]){
                                  tempwindy = 5;
                              } else {
                                  tempwindy = 7;
                              }
                              hx = H;
                          }

                          if (i == p) {
                               //tic=omp_get_wtime();
                              cartto2CG(i+K,j,N+1,K,DIRkLINnbx2);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              dirnbx2=DIRkLINnbx2[0];
                              knbx2=DIRkLINnbx2[1];
                              linnbx2=DIRkLINnbx2[2];
                              
                              m2 =Ucoarse[dirnbx2][knbx2][linnbx2];
                              hy = H;
                              tempwindx = 8;
                          } else if (i == M - 1 - (K - p) && p != 0) {
                               //tic=omp_get_wtime();
                              cartto2CG(i-K,j,N+1,K,DIRkLINnbx1);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              dirnbx1=DIRkLINnbx1[0];
                              knbx1=DIRkLINnbx1[1];
                              linnbx1=DIRkLINnbx1[2];
                              
                              m2 =Ucoarse[dirnbx1][knbx1][linnbx1];
                              
                              tempwindx = 6;
                              hy = H;
                          } else if (i == M - 1) {
                               //tic=omp_get_wtime();
                              cartto2CG(i-K,j,N+1,K,DIRkLINnbx1);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              dirnbx1=DIRkLINnbx1[0];
                              knbx1=DIRkLINnbx1[1];
                              linnbx1=DIRkLINnbx1[2];
                              
                              m2 =Ucoarse[dirnbx1][knbx1][linnbx1];
                              tempwindx = 6;
                              hy = H;
                          } else {
                              //tic=omp_get_wtime();
                              cartto2CG(i-K,j,N+1,K,DIRkLINnbx1);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              dirnbx1=DIRkLINnbx1[0];
                              knbx1=DIRkLINnbx1[1];
                              linnbx1=DIRkLINnbx1[2];
                               //tic=omp_get_wtime();
                              cartto2CG(i+K,j,N+1,K,DIRkLINnbx2);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              dirnbx2=DIRkLINnbx2[0];
                              knbx2=DIRkLINnbx2[1];
                              linnbx2=DIRkLINnbx2[2];
                              m2 =min(Ucoarse[dirnbx1][knbx1][linnbx1],Ucoarse[dirnbx2][knbx2][linnbx2]);
                              hy = H;
                              if (Ucoarse[dirnbx1][knbx1][linnbx1] < Ucoarse[dirnbx2][knbx2][linnbx2]){
                                  tempwindx = 6;
                              } else {
                                  tempwindx = 8;
                              }
                          }

                            if (m1 < 500 || m2 < 500) {
                                F = R[i + j * M];


                                discr = 2 * pow(F, 2) * pow(H, 2) - pow(m1 - m2, 2);

                                if (abs(m1 - m2) < F * H) {
                                    unew = 0.5 * (m1 + m2 + sqrt(discr));
                                    old =  Ucoarse[dir][k][lin];
                                    if (old < unew) {
                                        currError=0;
                                        Ucoarse[dir][k][lin] = old;
                                        maxError=max(currError,maxError);
                                    } else {

                                        if (tempwindy == 0) {
                                            wind[dir][k][lin] = tempwindx;
                                        } else if (tempwindx == 0) {
                                            wind[dir][k][lin] = tempwindy;
                                        } else {
                                            if (tempwindx == 6 && tempwindy == 5) {
                                                wind[dir][k][lin] = 2;
                                            } else if (tempwindx == 6 && tempwindy == 7) {
                                                wind[dir][k][lin] = 3;
                                            } else if (tempwindx == 8 && tempwindy == 5) {
                                                wind[dir][k][lin]= 1;
                                            } else if (tempwindx == 8 && tempwindy == 7) {
                                                wind[dir][k][lin]= 4;
                                            }
                                        }
                                        currError=abs(unew-old);
                                        Ucoarse[dir][k][lin]= unew;
                                         maxError=max(currError,maxError);
                                    }

                                } else {
                                    unew = min(m1, m2) + F * H;
                                    old =  Ucoarse[dir][k][lin];
                                    currError=0;
                                    Ucoarse[dir][k][lin] = old;
                                    maxError=max(currError,maxError);
                                    if (old < unew) {
                                        currError=0;
                                        Ucoarse[dir][k][lin] = old;
                                        maxError=max(currError,maxError);
                                    } else {
                                        if (m1 < m2) {
                                            wind[dir][k][lin] = tempwindy;
                                        } else {
                                            wind[dir][k][lin] = tempwindx;
                                        }

                                        currError=abs(unew-old);
                                        Ucoarse[dir][k][lin]= unew;
                                        maxError=max(currError,maxError);
                                    }


                                }
                              }
                            }
                          }


                      }
                    }
                }
        if (maxError<1e-12){
            totalsweeps[0]=iter;
            return;
        }
    }

    totalsweeps[0]=sweepiter;
   
    return;
}




//#endif //EIKONALPROJECT_COARSEINITIAL_H
