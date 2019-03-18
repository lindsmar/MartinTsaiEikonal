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

inline void cartto2CG(int i, int j, int k, int N, int K,int *DHBlin){
    
        DHBlin[0]=i%K;
        DHBlin[1]=j%K;
        DHBlin[2]=k%K;
        DHBlin[3]=i/K+j/K*N+k/K*N*N;

    return;
}


void coarseinitial(double *Ucoarse[DIM1][DIM2][DIM3], double *R, double *BC,double *obs, int N,int K, int p, int z, int r, double *wind[DIM1][DIM2][DIM3], int sweepiter, int *totalsweeps)
{
    int M;
    double H,h,old,m1,m2,m3,M1,M2,M3,discr,hx,hy,unew,F;
    double m11,m12,m21,m22,m31,m32;
    int tempwindx, tempwindy,tempwindz,istart1,istart2,jstart1,jstart2,kstart1,kstart2;
    double currError;

    H=(double)1/(double) N;
    h=(double)1/(double)(N*K);
    M=N*K+1;
    
    int s1,s2,s3,i,j,k;

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
    
    if (r==0){
        kstart1=M-1;
        kstart2=0;
    }
    else{
        kstart1=M-1-(K-r);
        kstart2=r;
    }
  
    double carttime=0;
    double tic,toc;
    
    
    int DHBlin[4];
    int Hs,Bs,Ds,lin;
    int DHBlinnb[4];
    int Hsnb,Bsnb,Dsnb,linnb;
   
    for (int iter=0; iter<sweepiter; iter++) {
        double maxError=0;
        for (s1 = 1; s1 >= -1; s1 -= 2){
            for (s2 = 1; s2 >= -1; s2 -= 2){
                 for (s3 = 1; s3 >= -1; s3 -= 2){
                for (i = (s1 < 0 ? istart1 : istart2); (s1 < 0 ? i >= 0 : i <= M - 1); i += K * s1){
                    for (j = (s2 < 0 ? jstart1 : jstart2); (s2 < 0 ? j >= 0 : j <= M - 1); j += K * s2) {
                         for (k = (s3 < 0 ? kstart1 : kstart2); (s3 < 0 ? k >= 0 : k <= M - 1); k += K * s3) {
                        if (obs[i + j * M + k*M*M] > 1) {
                            //don't do anything; point is in obstacle
                        } else if (BC[i + j * M + k*M*M] >= 0) {
                            //tic=omp_get_wtime();
                            cartto2CG(i,j,k,N+1,K,DHBlin);
                            //toc=omp_get_wtime()-tic;
                            //toc;
                            Ds=DHBlin[0];
                            Hs=DHBlin[1];
                            Bs=DHBlin[2];
                            lin=DHBlin[3];
                            Ucoarse[Ds][Hs][Bs][lin] = BC[i + j * M + k*M*M];
                            
                            
                        } else if (BC[i + j * M + k*M*M] < 0) {
                           // tic=omp_get_wtime();
                            cartto2CG(i,j,k,N+1,K,DHBlin);
                            //toc=omp_get_wtime()-tic;
                            //carttime+=toc;
                            Ds=DHBlin[0];
                            Hs=DHBlin[1];
                            Bs=DHBlin[2];
                            lin=DHBlin[3];
                            
                           
                            
                          tempwindx = 0;
                          tempwindy = 0;
                           tempwindz = 0;
                          if (j == z) {
                              //tic=omp_get_wtime();
                              cartto2CG(i,j+K,k,N+1,K,DHBlinnb);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              Dsnb=DHBlinnb[0];
                              Hsnb=DHBlinnb[1];
                              Bsnb=DHBlinnb[2];
                              linnb=DHBlinnb[3];
                              m1 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                              tempwindy = 7;
                              
                          } else if (j == M - 1 - (K - z) && z != 0) {
                              // tic=omp_get_wtime();
                              cartto2CG(i,j-K,k,N+1,K,DHBlinnb);
                              //toc=omp_get_wtime()-tic;
                             // carttime+=toc;
                              Dsnb=DHBlinnb[0];
                              Hsnb=DHBlinnb[1];
                              Bsnb=DHBlinnb[2];
                              linnb=DHBlinnb[3];
                              m1 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                              
                              tempwindy = 5;
                          } else if (j == M - 1) {
                             // tic=omp_get_wtime();
                              cartto2CG(i,j-K,k,N+1,K,DHBlinnb);
                             // toc=omp_get_wtime()-tic;
                             // carttime+=toc;
                              Dsnb=DHBlinnb[0];
                              Hsnb=DHBlinnb[1];
                              Bsnb=DHBlinnb[2];
                              linnb=DHBlinnb[3];
                              m1 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                              
                             
                              
                              tempwindy = 5;
                          } else {
                             // tic=omp_get_wtime();
                              cartto2CG(i,j+K,k,N+1,K,DHBlinnb);
                             // toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              Dsnb=DHBlinnb[0];
                              Hsnb=DHBlinnb[1];
                              Bsnb=DHBlinnb[2];
                              linnb=DHBlinnb[3];
                              m11 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                              // tic=omp_get_wtime();
                              cartto2CG(i,j-K,k,N+1,K,DHBlinnb);
                              //toc=omp_get_wtime()-tic;
                             // carttime+=toc;
                              Dsnb=DHBlinnb[0];
                              Hsnb=DHBlinnb[1];
                              Bsnb=DHBlinnb[2];
                              linnb=DHBlinnb[3];
                              m12 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                              m1 = min(m11,m12);
                              
                              
                          }

                          if (i == p) {
                               //tic=omp_get_wtime();
                              cartto2CG(i+K,j,k,N+1,K,DHBlinnb);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              Dsnb=DHBlinnb[0];
                              Hsnb=DHBlinnb[1];
                              Bsnb=DHBlinnb[2];
                              linnb=DHBlinnb[3];
                              m2 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                              
                             
                              tempwindx = 8;
                          } else if (i == M - 1 - (K - p) && p != 0) {
                               //tic=omp_get_wtime();
                              cartto2CG(i-K,j,k,N+1,K,DHBlinnb);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              Dsnb=DHBlinnb[0];
                              Hsnb=DHBlinnb[1];
                              Bsnb=DHBlinnb[2];
                              linnb=DHBlinnb[3];
                              m2 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                              
                              tempwindx = 6;
                              
                          } else if (i == M - 1) {
                               //tic=omp_get_wtime();
                              cartto2CG(i-K,j,k,N+1,K,DHBlinnb);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              Dsnb=DHBlinnb[0];
                              Hsnb=DHBlinnb[1];
                              Bsnb=DHBlinnb[2];
                              linnb=DHBlinnb[3];
                              m2 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                              tempwindx = 6;
                             
                          } else {
                              //tic=omp_get_wtime();
                              cartto2CG(i-K,j,k,N+1,K,DHBlinnb);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              Dsnb=DHBlinnb[0];
                              Hsnb=DHBlinnb[1];
                              Bsnb=DHBlinnb[2];
                              linnb=DHBlinnb[3];
                              m21 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                               //tic=omp_get_wtime();
                              cartto2CG(i+K,j,k,N+1,K,DHBlinnb);
                              //toc=omp_get_wtime()-tic;
                              //carttime+=toc;
                              Dsnb=DHBlinnb[0];
                              Hsnb=DHBlinnb[1];
                              Bsnb=DHBlinnb[2];
                              linnb=DHBlinnb[3];
                              m22 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                              m2 =min(m21,m22);
                              
                             
                          }
                            
                            if (k == r) {
                                //tic=omp_get_wtime();
                                cartto2CG(i,j,k+K,N+1,K,DHBlinnb);
                                //toc=omp_get_wtime()-tic;
                                //carttime+=toc;
                                Dsnb=DHBlinnb[0];
                                Hsnb=DHBlinnb[1];
                                Bsnb=DHBlinnb[2];
                                linnb=DHBlinnb[3];
                                m3= Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                
                                tempwindx = 8;
                            } else if (k == M - 1 - (K - r) && r != 0) {
                                //tic=omp_get_wtime();
                                cartto2CG(i,j,k-K,N+1,K,DHBlinnb);
                                //toc=omp_get_wtime()-tic;
                                //carttime+=toc;
                                Dsnb=DHBlinnb[0];
                                Hsnb=DHBlinnb[1];
                                Bsnb=DHBlinnb[2];
                                linnb=DHBlinnb[3];
                                m3= Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                
                                tempwindx = 6;
                                
                            } else if (k == M - 1) {
                                //tic=omp_get_wtime();
                                cartto2CG(i,j,k-K,N+1,K,DHBlinnb);
                                //toc=omp_get_wtime()-tic;
                                //carttime+=toc;
                                Dsnb=DHBlinnb[0];
                                Hsnb=DHBlinnb[1];
                                Bsnb=DHBlinnb[2];
                                linnb=DHBlinnb[3];
                                m3= Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                tempwindx = 6;
                                
                            } else {
                                //tic=omp_get_wtime();
                                cartto2CG(i,j,k-K,N+1,K,DHBlinnb);
                                //toc=omp_get_wtime()-tic;
                                //carttime+=toc;
                                Dsnb=DHBlinnb[0];
                                Hsnb=DHBlinnb[1];
                                Bsnb=DHBlinnb[2];
                                linnb=DHBlinnb[3];
                                m31= Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                //tic=omp_get_wtime();
                                cartto2CG(i,j,k+K,N+1,K,DHBlinnb);
                                //toc=omp_get_wtime()-tic;
                                //carttime+=toc;
                                Dsnb=DHBlinnb[0];
                                Hsnb=DHBlinnb[1];
                                Bsnb=DHBlinnb[2];
                                linnb=DHBlinnb[3];
                                m32= Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                m3 =min(m31,m32);
                               
                               
                            }

                            
                            
                            if (m1 <= m2){
                                if (m1 <= m3){
                                    M1=m1;
                                    if (m2<=m3)
                                    {
                                        M2=m2;
                                        M3=m3;
                                    }
                                    else{
                                        M2=m3;
                                        M3=m2;
                                    }
                                }
                                else{
                                    M1=m3;
                                    M2=m1;
                                    M3=m2;
                                }
                            }
                            else{
                                if (m2<=m3){
                                    M1=m2;
                                    if (m1<=m3){
                                        M2=m1;
                                        M3=m3;
                                    }
                                    else{
                                        M2=m3;
                                        M3=m1;
                                    }
                                }
                                else{
                                    M1=m3;
                                    M2=m2;
                                    M3=m1;
                                }
                            }
                            
                            unew=M1+R[i + j * M + k*M*M]*H;
                            
                            if (unew <= M2){
                                old=Ucoarse[Ds][Hs][Bs][lin];
                                if (old < unew) {
                                    currError=0;
                                    Ucoarse[Ds][Hs][Bs][lin] = old;
                                    maxError=max(currError,maxError);
                                } else{
                                    /*if (tempwindy == 0) {
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
                                    }*/
                                    
                                    currError=abs(unew-old);
                                    Ucoarse[Ds][Hs][Bs][lin]= unew;
                                    maxError=max(currError,maxError);
                                    
                                }
                            }
                            else{
                                discr=2*pow(R[i+j*N+k*N*N],2)*pow(H,2)-pow(M1-M2,2);
                                unew=0.5*(M1+M2+sqrt(discr));
                                if (unew <= M3){
                                    old=Ucoarse[Ds][Hs][Bs][lin];
                                    if (old < unew) {
                                        currError=0;
                                        Ucoarse[Ds][Hs][Bs][lin] = old;
                                        maxError=max(currError,maxError);
                                    } else{
                                        /*if (tempwindy == 0) {
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
                                         }*/
                                        
                                        currError=abs(unew-old);
                                        Ucoarse[Ds][Hs][Bs][lin]= unew;
                                        maxError=max(currError,maxError);
                                        
                                    }
                                }
                                else{
                                    discr=4*pow(M1+M2+M3,2)-12*(M1*M1+M2*M2+M3*M3-pow(R[i+j*N+k*N*N],2)*pow(H,2));
                                    unew=(2*(M1+M2+M3)+sqrt(discr))/6;
                                    if (unew<1000){
                                        old=Ucoarse[Ds][Hs][Bs][lin];
                                        if (old < unew) {
                                            currError=0;
                                            Ucoarse[Ds][Hs][Bs][lin] = old;
                                            maxError=max(currError,maxError);
                                        } else{
                                            /*if (tempwindy == 0) {
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
                                             }*/
                                            
                                            currError=abs(unew-old);
                                            Ucoarse[Ds][Hs][Bs][lin]= unew;
                                            maxError=max(currError,maxError);
                                            
                                        }
                                    }
                                }
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
