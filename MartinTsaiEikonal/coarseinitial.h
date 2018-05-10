//
// Created by Lindsay Martin on 10/29/16.
//

#ifndef EIKONALPROJECT_COARSEINITIAL_H
#define EIKONALPROJECT_COARSEINITIAL_H
//
//  coarseupdate.h
//  Eikonal Project 3
//
//  Created by Lindsay Martin on 9/29/16.
//  Copyright © 2016 Lindsay Martin. All rights reserved.
//





#include <algorithm>
#include <cmath>
using namespace std;


double * coarseinitial(double *Ucoarse,double * R, double *BC,double *obs, int N,int K, int p, int z, double *wind, int sweepiter);

double * coarseinitial(double *Ucoarse, double *R, double *BC,double *obs, int N,int K, int p, int z, double *wind, int sweepiter)
{
    int M;
    double H,h,old,m1,m2,discr,hx,hy,unew,F;
    int tempwindx, tempwindy,istart1,istart2,jstart1,jstart2;;


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

    for (int iter=0; iter<sweepiter; iter++) {
        for (s1 = 1; s1 >= -1; s1 -= 2){
            for (s2 = 1; s2 >= -1; s2 -= 2){

                for (i = (s1 < 0 ? istart1 : istart2); (s1 < 0 ? i >= 0 : i <= M - 1); i += K * s1){
                    for (j = (s2 < 0 ? jstart1 : jstart2); (s2 < 0 ? j >= 0 : j <= M - 1); j += K * s2) {
                        if (obs[i + j * M] > 1) {
                            //don't do anything; point is in obstacle
                        } else if (BC[i + j * M] >= 0) {
                            Ucoarse[i + j * M] = BC[i + j * M];
                        } else if (BC[i + j * M] < 0) {
                          tempwindx = 0;
                          tempwindy = 0;
                          if (j == z) {
                              m1 = Ucoarse[i + (j + K) * M];
                              tempwindy = 7;
                              hx = H;
                          } else if (j == M - 1 - (K - z) && z != 0) {
                              m1 = Ucoarse[i + (j - K) * M];
                              hx = H;
                              tempwindy = 5;
                          } else if (j == M - 1) {
                              m1 = Ucoarse[i + (j - K) * M];
                              hx = H;
                              tempwindy = 5;
                          } else {
                              m1 = min(Ucoarse[i + (j - K) * M], Ucoarse[i + (j + K) * M]);
                              if (Ucoarse[i + (j - K) * M] < Ucoarse[i + (j + K) * M]){
                                  tempwindy = 5;
                              } else {
                                  tempwindy = 7;
                              }
                              hx = H;
                          }

                          if (i == p) {
                              m2 = Ucoarse[(i + K) + j * M];
                              hy = H;
                              tempwindx = 8;
                          } else if (i == M - 1 - (K - p) && p != 0) {
                              m2 = Ucoarse[(i - K) + j * M];
                              tempwindx = 6;
                              hy = H;
                          } else if (i == M - 1) {
                              m2 = Ucoarse[(i - K) + j * M];
                              tempwindx = 6;
                              hy = H;
                          } else {
                              m2 =min(Ucoarse[(i - K) + j * M],Ucoarse[(i + K) + j * M]);
                              hy = H;
                              if (Ucoarse[(i - K) + j * M] < Ucoarse[(i + K) + j * M]){
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
                                    old = Ucoarse[i + j * M];
                                    if (old < unew) {
                                        Ucoarse[i + j * M] = old;
                                    } else {

                                        if (tempwindy == 0) {
                                            wind[i + j * M] = tempwindx;
                                        } else if (tempwindx == 0) {
                                            wind[i + j * M] = tempwindy;
                                        } else {
                                            if (tempwindx == 6 && tempwindy == 5) {
                                                wind[i + j * M] = 2;
                                            } else if (tempwindx == 6 && tempwindy == 7) {
                                                wind[i + j * M] = 3;
                                            } else if (tempwindx == 8 && tempwindy == 5) {
                                                wind[i + j * M] = 1;
                                            } else if (tempwindx == 8 && tempwindy == 7) {
                                                wind[i + j * M] = 4;
                                            }
                                        }

                                        Ucoarse[i + j * M] = unew;
                                    }

                                } else {
                                    unew = min(m1, m2) + F * H;
                                    old = Ucoarse[i + j * M];
                                    if (old < unew) {
                                        Ucoarse[i + j * M] = old;
                                    } else {
                                        if (m1 < m2) {
                                            wind[i + j * M] = tempwindy;
                                        } else {
                                            wind[i + j * M] = tempwindx;
                                        }

                                        Ucoarse[i + j * M] = unew;
                                    }


                                }
                              }
                            }
                          }


                      }
                    }
                }
    }





    return Ucoarse;
}




#endif //EIKONALPROJECT_COARSEINITIAL_H
