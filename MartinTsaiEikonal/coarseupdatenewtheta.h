//
// Created by Lindsay Martin on 09/06/17.
//

#ifndef EIKONALPROJECT_COARSEUPDATENEWTHETA_H
#define EIKONALPROJECT_COARSEUPDATENEWTHETA_H
//
// Created by Lindsay Martin on 10/23/16.
//

#include "consistent.h"
#include "smooth.h"

double * coarseupdatenewtheta(double *Ucoarse, double *Ucoarse1, double *Ufine1, double *Ufine2,  double *UK1, double *UK2, double *UK3,  double * R, double *BC,double *obs, int N,int K, int p, int z, double *wind, double *fwind, int sweepiter, double epsilon, double delta, double x0, double *theta, int wt1, int wt2, int wt3);

double * coarseupdatenewtheta(double *Ucoarse, double *Ucoarse1, double *Ufine1, double *Ufine2,  double *UK1, double *UK2, double *UK3,  double * R, double *BC,double *obs, int N,int K, int p, int z, double *wind, double *fwind, int sweepiter, double epsilon, double  delta, double x0, double *theta, int wt1, int wt2, int wt3)
{
    int M;
    double H,h,old,m1,m2,discr,hx,hy,unew,F, utemp, uknew;
    int tempwindx, tempwindy, istart1,istart2,jstart1,jstart2;
    double prevW, newwind,prevw;

    bool consis;

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
                                m1 = min(Ucoarse[i + (j - K) * M],Ucoarse[i + (j + K) * M]);
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
                                m2 =min(Ucoarse[(i - K) + j * M], Ucoarse[(i + K) + j * M]);
                                hy = H;
                                if (Ucoarse[(i - K) + j * M] < Ucoarse[(i + K) + j * M]){
                                    tempwindx = 6;
                                } else {
                                    tempwindx = 8;
                                }
                            }




                            if (m1 < 500 || m2 < 500) {
                                F = R[i + j * M];
                                prevW=wind[i + j * M];
                                prevw=fwind[i + j * M];
                                consis=consistent(prevW,prevw,p,z);


                                discr = 2 * pow(F, 2) * pow(H, 2) - pow(m1 - m2, 2);

                                if (abs(m1 - m2) < F * H) {

                                    if (tempwindx == 6 && tempwindy == 5) {
                                        newwind = 2;
                                    } else if (tempwindx == 6 && tempwindy == 7) {
                                        newwind = 3;
                                    } else if (tempwindx == 8 && tempwindy == 5) {
                                        newwind = 1;
                                    } else if (tempwindx == 8 && tempwindy == 7) {
                                        newwind = 4;
                                    }
                                    /*unew = 0.5 * (m1 + m2 + sqrt(discr));
                                    if (abs(unew-UK[i+j*M])<1e-15){
                                      utemp = Ufine[i + j * M];
                                      Ucoarse[i + j * M] = utemp;
                                       UKtemp[i + j * M] = unew;
                                       tempwind[i + j * M] = newwind;
                                      continue;
                                    }*/



                                    if (marked[i+j*M]==0){
                                      if (consis){
                                        if (prevW==newwind){
                                          unew = 0.5 * (m1 + m2 + sqrt(discr));
                                          if (abs(unew-UK1[i+j*M])<1e-15){ //|| abs(Ufine1[i+j*M]-UK1[i+j*M])<1e-15)
                                            utemp = Ufine1[i + j * M];
                                          }
                                          else{
                                              wtdsum=(wt1*(unew-UK1[i+j*M])+wt2*(UK1[i+j*M]-UK2[i+j*M])+wt3*(UK2[i+j*M]-UK3[i+j*M]))/(wt1+wt2+wt3);
                                            T=(Ufine1[i+j*M]-Ufine2[i+j*M])/wtdsum;
                                            //T=(ORIG[i+j*M]-Ufine[i+j*M])/(unew-UK[i+j*M]);
                                            //T=T-.01;
                                            lower=(Ucoarse1[i+j*M]-Ufine1[i+j*M])/(unew-UK1[i+j*M]);
                                            if (T<0){
                                              T=0;
                                            }
                                            else{
                                              T=smooth(T,epsilon,delta,x0);
                                            }
                                            theta[i+j*M]=T;
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
                                          //utemp =0.5 * (m1 + m2 + sqrt(discr));
                                          //uknew = utemp;
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
                                  //    old = Ucoarse[i + j * M];
                                  //    if (old <= utemp) {
                                  //        Ucoarse[i + j * M] = old;
                                  //    } else {
                                          Ucoarse[i + j * M] = utemp;
                                          UKtemp1[i + j * M] = unew;

                                          tempwind[i + j * M] = newwind;
                                        //}


                                    }
                                    else{
                                      //do nothing
                                    }


                                } else {
                                    if (m1 < m2) {
                                        newwind = tempwindy;
                                    } else {
                                        newwind = tempwindx;
                                    }
                                    /*unew = min(m1, m2) + F * H;
                                    if (abs(unew-UK[i+j*M])<1e-15){
                                      utemp = Ufine[i + j * M];
                                      Ucoarse[i + j * M] = utemp;
                                       UKtemp[i + j * M] = unew;
                                       tempwind[i + j * M] = newwind;
                                      continue;
                                    }*/

                                    if (marked[i+j*M]==0){
                                      if (consis){
                                        if (prevW==newwind){
                                          unew = min(m1, m2) + F * H;
                                          if (abs(unew-UK1[i+j*M])<1e-15){// || abs(Ufine1[i+j*M]-UK1[i+j*M])<1e-15)
                                            utemp = Ufine1[i + j * M];
                                          }
                                          else{
                                            wtdsum=(wt1*(unew-UK1[i+j*M])+wt2*(UK1[i+j*M]-UK2[i+j*M])+wt3*(UK2[i+j*M]-UK3[i+j*M]))/(wt1+wt2+wt3);
                                            T=(Ufine1[i+j*M]-Ufine2[i+j*M])/wtdsum;
                                            //T=(ORIG[i+j*M]-Ufine[i+j*M])/(unew-UK[i+j*M]);
                                            //T=T-.01;
                                            lower=(Ucoarse1[i+j*M]-Ufine1[i+j*M])/(unew-UK1[i+j*M]);
                                            if (T<0){
                                              T=0;
                                            }
                                            else{
                                              T=smooth(T,epsilon,delta,x0);
                                            }
                                            theta[i+j*M]=T;
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
                                          //utemp =0.5 * (m1 + m2 + sqrt(discr));
                                          //uknew = utemp;
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
                                    //  old = Ucoarse[i + j * M];
                                    //  if (old <= utemp) {
                                    //      Ucoarse[i + j * M] = old;
                                  //    } else {
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
                }


    }

    for (i = istart2; i <= M - 1; i += K){
        for (j = jstart2;  j <= M - 1; j += K) {
          UK3[i+j*M]=UK2[i+j*M];
          UK2[i+j*M]=UK1[i+j*M];
          UK1[i+j*M]=UKtemp1[i+j*M];
            wind[i+j*M]=tempwind[i+j*M];
          }
      }



mxFree(marked);
mxFree(UKtemp1);
mxFree(tempwind);
    return Ucoarse;
}


#endif //EIKONALPROJECT_COARSEUPDATENEWTHETA_H
