//
// Created by Lindsay Martin on 09/06/17.
//

#ifndef EIKONALPROJECT_COARSEUPDATE2_ani_H
#define EIKONALPROJECT_COARSEUPDATE2_ani_H
//
// Created by Lindsay Martin on 10/23/16.
//

#include "consistent.h"
#include "smooth.h"

double * coarseupdate_ani2(double *Ucoarse,double *Ucoarse1, double *Ufine1, double *Ufine2,  double *UK1, double *UK2, double *UK3, double  * cbar, double *BC,double *obs, int N,int K,int p, int z, double *wind, double *fwind,  int sweepiter, double epsilon, double delta, double x0, int wt1, int wt2, int wt3);

double * coarseupdate_ani2(double *Ucoarse,double *Ucoarse1, double *Ufine1, double *Ufine2,  double *UK1, double *UK2, double *UK3, double *cbar, double *BC,double *obs, int N,int K, int p, int z, double *wind, double *fwind,  int sweepiter, double epsilon, double delta, double x0, int wt1, int wt2, int wt3)
{
    int M;
    double H,h,unew, utemp, uknew;
    int istart1,istart2,jstart1,jstart2;
    double prevW, newwind,prevw, tempmin;
    double neighbors[16];
    bool consis;
    int index;
    H=(double)1/(double) N;
    h=(double)1/(double)(N*K);
    M=N*K+1;


    double T, wtdsum, lower;
    double *marked, *UKtemp1, *tempwind;

    marked= (double *) mxCalloc(M*M, sizeof(double));
    UKtemp1= (double *) mxCalloc(M*M, sizeof(double));
    tempwind= (double *) mxCalloc(M*M, sizeof(double));
    for (int c=0; c<M*M; c++){
        marked[c]=0;
        tempwind[c]=0;
    }

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
       for (s2 = -1; s2 <= 1; s2 += 2){
        for (s1 = -1; s1 <= 1; s1 += 2){
          //for (s2 = 1; s2 >= -1; s2 -= 2)
          // for (s1 = 1; s1 >= -1; s1 -= 2)
        //s2=1;
        //s1=1;
                for (i = (s1 < 0 ? istart1 : istart2); (s1 < 0 ? i >= 0 : i <= M - 1); i += K * s1){
                    for (j = (s2 < 0 ? jstart1 : jstart2); (s2 < 0 ? j >= 0 : j <= M - 1); j += K * s2) {
                        if (obs[i + j * M] > 1) {
                            //don't do anything; point is in obstacle
                        } else if (BC[i + j * M] >= 0) {
                            Ucoarse[i + j * M] = BC[i + j * M];
                        }else if (BC[i + j * M] < 0){
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
                        prevW=wind[i + j * M];
                        prevw=fwind[i + j * M];
                        consis=consistent(prevW,prevw,p,z);
                        newwind=neighbors[index+8];
                        if (marked[i+j*M]==0){
                          if (consis){
                            if (prevW==newwind){
                              unew=tempmin;
                              if (abs(unew-UK1[i+j*M])<1e-15){// || abs(Ufine[i+j*M]-UK[i+j*M])<1e-13)
                                utemp = Ufine1[i + j * M];
                              }
                              else{

                                wtdsum=(wt1*(unew-UK1[i+j*M])+wt2*(UK1[i+j*M]-UK2[i+j*M])+wt3*(UK2[i+j*M]-UK3[i+j*M]))/(wt1+wt2+wt3);
                              T=(Ufine1[i+j*M]-Ufine2[i+j*M])/wtdsum;
                              //T=(ORIG[i+j*M]-Ufine[i+j*M])/(unew-UK[i+j*M]);
                              //T=T-.01;
                              lower=(Ucoarse1[i+j*M]-Ufine1[i+j*M])/(unew-UK1[i+j*M]);
                              if (T<0){
                                T=0;//lower+1;
                              }
                              else{
                                T=smooth(T,epsilon,delta,x0);
                              }
                              utemp = T * unew + Ufine1[i + j * M] - T * UK1[i + j * M];
                              }
                              uknew=unew;
                              //marked[i+j*M]=1;

                              Ucoarse[i + j * M] = utemp;
                              UKtemp1[i + j * M] = unew;
                              tempwind[i + j * M] = fwind[i+j*M];
                              continue;
                            }
                            else{
                              utemp=Ufine1[i + j * M];
                              uknew=utemp;
                              newwind=fwind[i+j*M];
                            }
                          }
                          else{

                            utemp=Ufine1[i + j * M];
                            uknew=utemp;
                            newwind=fwind[i+j*M];
                          }
                          //old = Ucoarse[i + j * M];
                          //if (old <= utemp) {
                          //    Ucoarse[i + j * M] = old;
                        //  } else {
                            Ucoarse[i + j * M] = utemp;
                            UKtemp1[i + j * M] = unew;

                            tempwind[i + j * M] = newwind;
                        //  }


                        }
                        else{
                          //do nothing
                        }
                      }
                    }
                  }
                }
              }
            }




    for (i = istart2; i <= M - 1; i += K)
        for (j = jstart2;  j <= M - 1; j += K) {
          UK3[i+j*M]=UK2[i+j*M];
          UK2[i+j*M]=UK1[i+j*M];
          UK1[i+j*M]=UKtemp1[i+j*M];
            wind[i+j*M]=tempwind[i+j*M];
    }

mxFree(marked);
mxFree(UKtemp1);
mxFree(tempwind);
    return Ucoarse;
}


#endif //EIKONALPROJECT_COARSEUPDATE2_H
