//
//

#ifndef EIKONALPROJECT_COARSEINITIAL_ani_H
#define EIKONALPROJECT_COARSEINITIAL_ani_H
//
//  coarseinitial_ani.h
//  Eikonal Project 3
//
//  Copyright © 2016 Lindsay Martin. All rights reserved.
//





#include <algorithm>
#include <cmath>
using namespace std;


double * coarseinitial_ani(double *Ucoarse,double * cbar, double *BC,double *obs, int N,int K, int p, int z, double *wind, int sweepiter);

double * coarseinitial_ani(double *Ucoarse, double *cbar, double *BC,double *obs, int N,int K, int p, int z, double *wind, int sweepiter)
{
    int M;
    double H,h,old,unew;
    int istart1,istart2,jstart1,jstart2;;
    double neighbors[16];
    int index;
    double tempmin,tempwind;

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

                          old=Ucoarse[i+j*M];
                          if (j == jstart2 && i==istart2){
                        neighbors[0]=1000;
                        neighbors[1]=1000;
                        neighbors[2]=H/cbar[0]+Ucoarse[i+(j+K)*M];
                        neighbors[3]=sqrt(2)*H/cbar[7]+Ucoarse[(i+K)+(j+K)*M];
                        neighbors[4]=H/cbar[6]+Ucoarse[(i+K)+j*M];
                        neighbors[5]=1000;
                        neighbors[6]=1000;
                        neighbors[7]=1000;
                        neighbors[8]=6;
                        neighbors[9]=3;
                        neighbors[10]=7;
                        neighbors[11]=4;
                        neighbors[12]=8;
                        neighbors[13]=1;
                        neighbors[14]=5;
                        neighbors[15]=2;
                      }
                      else if (j==jstart2 && i!=istart2 && i!=istart1){
                        neighbors[0]=H/cbar[2]+Ucoarse[(i-K)+j*M];
                        neighbors[1]=sqrt(2)*H/cbar[1]+Ucoarse[(i-K)+(j+K)*M];
                        neighbors[2]=H/cbar[0]+Ucoarse[i+(j+K)*M];
                        neighbors[3]=sqrt(2)*H/cbar[7]+Ucoarse[(i+K)+(j+K)*M];
                        neighbors[4]=H/cbar[6]+Ucoarse[(i+K)+j*M];
                        neighbors[5]=1000;
                        neighbors[6]=1000;
                        neighbors[7]=1000;
                        neighbors[8]=6;
                        neighbors[9]=3;
                        neighbors[10]=7;
                        neighbors[11]=4;
                        neighbors[12]=8;
                        neighbors[13]=1;
                        neighbors[14]=5;
                        neighbors[15]=2;
                      }
                      else if (j==jstart2 && i==istart1){
                        neighbors[0]=H/cbar[2]+Ucoarse[(i-K)+j*M];
                        neighbors[1]=sqrt(2)*H/cbar[1]+Ucoarse[(i-K)+(j+K)*M];
                        neighbors[2]=H/cbar[0]+Ucoarse[i+(j+K)*M];
                        neighbors[3]=1000;
                        neighbors[4]=1000;
                        neighbors[5]=1000;
                        neighbors[6]=1000;
                        neighbors[7]=1000;
                        neighbors[8]=6;
                        neighbors[9]=3;
                        neighbors[10]=7;
                        neighbors[11]=4;
                        neighbors[12]=8;
                        neighbors[13]=1;
                        neighbors[14]=5;
                        neighbors[15]=2;
                      }
                      else if (j!=jstart2 && j!=jstart1 && i==istart1){
                        neighbors[0]=H/cbar[2]+Ucoarse[(i-K)+j*M];
                        neighbors[1]=sqrt(2)*H/cbar[1]+Ucoarse[(i-K)+(j+K)*M];
                        neighbors[2]=H/cbar[0]+Ucoarse[i+(j+K)*M];
                        neighbors[3]=1000;
                        neighbors[4]=1000;
                        neighbors[5]=1000;
                        neighbors[6]=H/cbar[4]+Ucoarse[i+(j-K)*M];
                        neighbors[7]=sqrt(2)*H/cbar[3]+Ucoarse[(i-K)+(j-K)*M];
                        neighbors[8]=6;
                        neighbors[9]=3;
                        neighbors[10]=7;
                        neighbors[11]=4;
                        neighbors[12]=8;
                        neighbors[13]=1;
                        neighbors[14]=5;
                        neighbors[15]=2;
                      }
                      else if (j==jstart1 && i==istart1){
                        neighbors[0]=H/cbar[2]+Ucoarse[(i-K)+j*M];
                        neighbors[1]=1000;
                        neighbors[2]=1000;
                        neighbors[3]=1000;
                        neighbors[4]=1000;
                        neighbors[5]=1000;
                        neighbors[6]=H/cbar[4]+Ucoarse[i+(j-K)*M];
                        neighbors[7]=sqrt(2)*H/cbar[3]+Ucoarse[(i-K)+(j-K)*M];
                        neighbors[8]=6;
                        neighbors[9]=3;
                        neighbors[10]=7;
                        neighbors[11]=4;
                        neighbors[12]=8;
                        neighbors[13]=1;
                        neighbors[14]=5;
                        neighbors[15]=2;
                      }
                      else if (j==jstart1 && i!=istart2 && i!=istart1){
                        neighbors[0]=H/cbar[2]+Ucoarse[(i-K)+j*M];
                        neighbors[1]=1000;
                        neighbors[2]=1000;
                        neighbors[3]=1000;
                        neighbors[4]=H/cbar[6]+Ucoarse[(i+K)+j*M];
                        neighbors[5]=sqrt(2)*H/cbar[5]+Ucoarse[(i+K)+(j-K)*M];
                        neighbors[6]=H/cbar[4]+Ucoarse[i+(j-K)*M];
                        neighbors[7]=sqrt(2)*H/cbar[3]+Ucoarse[(i-K)+(j-K)*M];
                        neighbors[8]=6;
                        neighbors[9]=3;
                        neighbors[10]=7;
                        neighbors[11]=4;
                        neighbors[12]=8;
                        neighbors[13]=1;
                        neighbors[14]=5;
                        neighbors[15]=2;
                      }
                      else if (j==jstart1 && i==istart2){
                        neighbors[0]=1000;
                        neighbors[1]=1000;
                        neighbors[2]=1000;
                        neighbors[3]=1000;
                        neighbors[4]=H/cbar[6]+Ucoarse[(i+K)+j*M];
                        neighbors[5]=sqrt(2)*H/cbar[5]+Ucoarse[(i+K)+(j-K)*M];
                        neighbors[6]=H/cbar[4]+Ucoarse[i+(j-K)*M];
                        neighbors[7]=1000;
                        neighbors[8]=6;
                        neighbors[9]=3;
                        neighbors[10]=7;
                        neighbors[11]=4;
                        neighbors[12]=8;
                        neighbors[13]=1;
                        neighbors[14]=5;
                        neighbors[15]=2;
                      }
                      else if (j!=jstart1 && i==istart2 && j!=jstart2){
                        neighbors[0]=1000;
                        neighbors[1]=1000;
                        neighbors[2]=H/cbar[0]+Ucoarse[i+(j+K)*M];
                        neighbors[3]=sqrt(2)*H/cbar[7]+Ucoarse[(i+K)+(j+K)*M];
                        neighbors[4]=H/cbar[6]+Ucoarse[(i+K)+j*M];
                        neighbors[5]=sqrt(2)*H/cbar[5]+Ucoarse[(i+K)+(j-K)*M];
                        neighbors[6]=H/cbar[4]+Ucoarse[i+(j-K)*M];
                        neighbors[7]=1000;
                        neighbors[8]=6;
                        neighbors[9]=3;
                        neighbors[10]=7;
                        neighbors[11]=4;
                        neighbors[12]=8;
                        neighbors[13]=1;
                        neighbors[14]=5;
                        neighbors[15]=2;
                      }
                      else{
                        neighbors[0]=H/cbar[2]+Ucoarse[(i-K)+j*M];
                        neighbors[1]=sqrt(2)*H/cbar[1]+Ucoarse[(i-K)+(j+K)*M];
                        neighbors[2]=H/cbar[0]+Ucoarse[i+(j+K)*M];
                        neighbors[3]=sqrt(2)*H/cbar[7]+Ucoarse[(i+K)+(j+K)*M];
                        neighbors[4]=H/cbar[6]+Ucoarse[(i+K)+j*M];
                        neighbors[5]=sqrt(2)*H/cbar[5]+Ucoarse[(i+K)+(j-K)*M];
                        neighbors[6]=H/cbar[4]+Ucoarse[i+(j-K)*M];
                        neighbors[7]=sqrt(2)*H/cbar[3]+Ucoarse[(i-K)+(j-K)*M];
                        neighbors[8]=6;
                        neighbors[9]=3;
                        neighbors[10]=7;
                        neighbors[11]=4;
                        neighbors[12]=8;
                        neighbors[13]=1;
                        neighbors[14]=5;
                        neighbors[15]=2;
                      }
                      index=0;
                      for (int i=1;i<=7;i++){
                          if (neighbors[i]<neighbors[index]){
                          index=i;
                          tempmin=neighbors[index];
                        }
                        else{
                          tempmin=neighbors[index];
                        }
                      }
                      old=Ucoarse[i+j*M];
                      if (tempmin<old){
                        Ucoarse[i+j*M]=tempmin;
                        wind[i+j*M]=neighbors[index+8];
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
