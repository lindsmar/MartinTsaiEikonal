//
//  fineupdate2.h
//  Eikonal Project 5
//
//  Created by Lindsay Martin on 6/22/16.
//  Copyright © 2016 Lindsay Martin. All rights reserved.
//

#ifndef fineupdate2_h
#define fineupdate2_h
#include "whichvec.h"
#include <algorithm>
#include <cmath>
using namespace std;

void fineupdate2(double* FINE, double *Ucoarse, double *f, double *BC,double *obs, int N,int K, double *wind, double *fwind);
void fineupdate2(double * FINE, double *Ucoarse, double *f,double *BC,double *obs, int N,int K, double *wind, double *fwind)
{
    int M;
    int MM;
    double h;

    h=(double)1/(double) (N*K);
    M=N*K+1;
    MM=M+N-1;
    int s1;
    int s2;
    int cont, tempwindx,tempwindy, windtemp;


    //double *FINE, *FINE2, *exwind, *cellbdry;
    double *fine, *WIND, *BDRY, *UFM;

    //fine2=(double *) mxCalloc(M*M, sizeof(double));
    fine=(double *) mxCalloc(M*M, sizeof(double));
    //FINE=mxCreateDoubleMatrix(M,M,mxREAL);
    //fine=mxGetPr(FINE);
    //FINE2=mxCreateDoubleMatrix(M,M,mxREAL);
    //fine2=mxGetPr(FINE2);
    UFM=(double *) mxCalloc(MM*MM, sizeof(double));
    WIND=(double *) mxCalloc(MM*MM, sizeof(double));
    BDRY=(double *) mxCalloc(MM*MM, sizeof(double));
    int c;
    for (c=0; c<(M+N-1)*(M+N-1); c++){
        UFM[c]=1000;
    }
    /*exwind=mxCreateDoubleMatrix(MM,MM,mxREAL);
     WIND=mxGetPr(exwind);

     cellbdry=mxCreateDoubleMatrix(MM,MM,mxREAL);
     BDRY=mxGetPr(cellbdry);*/

    for (int z=0; z<M*M; z++){
        fine[z]=1000;
        //   FINE[z]=Ucoarse[z];
        fwind[z]=wind[z];
    }

    int z;
    for (z=0; z<(M+N-1)*(M+N-1); z++){
        WIND[z]=1000;
    }

    for (z=0; z<(M+N-1)*(M+N-1); z++){
        BDRY[z]=-1;
    }

    double unew,m1,m2,discr,old;
    int p,m,i,j,v;

    //Set Bdry condition for each subdomain
    for (m=0; m<=N-1; m++){
        for (p=0; p<=N-1; p++){
            for (i=p*K+p; i<=(p+1)*K+p; i++){
                for (j=m*K+m; j<=(m+1)*K+m; j++)
                {
                    if (obs[(i-p)+(j-m)*M]>1)
                    {
                        //don't do anything; point is in obstacle
                    }
                    else if (BC[(i-p)+(j-m)*M]>=0)
                    {
                        UFM[i+j*MM]=BC[(i-p)+(j-m)*M];
                        BDRY[i+j*MM]=BC[(i-p)+(j-m)*M];
                    }
                    else if (wind[(i-p)+(j-m)*M]!=0){

                        if (i==p*K+p && j==m*K+m){
                            if (wind[(i-p)+(j-m)*M]==5  ||  wind[(i-p)+(j-m)*M]==6 || wind[(i-p)+(j-m)*M]==3 || wind[(i-p)+(j-m)*M]==1 || wind[(i-p)+(j-m)*M]==2){
                                UFM[i+j*MM]=Ucoarse[(i-p)+(j-m)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[(i-p)+(j-m)*M];
                            }
                            else{
                                UFM[i+j*MM]=1000;
                            }
                        }
                        else if (i==p*K+p && j==(m+1)*K+m){
                            if (wind[(i-p)+(j-m)*M]==6  ||  wind[(i-p)+(j-m)*M]==7 || wind[(i-p)+(j-m)*M]==4 || wind[(i-p)+(j-m)*M]==2 || wind[(i-p)+(j-m)*M]==3){
                                UFM[i+j*MM]=Ucoarse[(i-p)+(j-m)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[(i-p)+(j-m)*M];
                            }
                            else{
                                UFM[i+j*MM]=1000;
                            }
                        }
                        else if (i==(p+1)*K+p && j==m*K+m){
                            if (wind[(i-p)+(j-m)*M]==8  ||  wind[(i-p)+(j-m)*M]==5 || wind[(i-p)+(j-m)*M]==4 || wind[(i-p)+(j-m)*M]==2|| wind[(i-p)+(j-m)*M]==1){
                                UFM[i+j*MM]=Ucoarse[(i-p)+(j-m)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[(i-p)+(j-m)*M];
                            }
                            else{
                                UFM[i+j*MM]=1000;
                            }
                        }
                        else if (i==(p+1)*K+p && j==(m+1)*K+m){
                            if (wind[(i-p)+(j-m)*M]==8  ||  wind[(i-p)+(j-m)*M]==7 || wind[(i-p)+(j-m)*M]==3 || wind[(i-p)+(j-m)*M]==1 || wind[(i-p)+(j-m)*M]==4){
                                UFM[i+j*MM]=Ucoarse[(i-p)+(j-m)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[(i-p)+(j-m)*M];
                            }
                            else{
                                UFM[i+j*MM]=1000;
                            }

                        }
                        else if (i==p*K+p){
                            if (wind[p*K+(j-m)*M]==2 || wind[p*K+(j-m)*M]==3 || wind[p*K+(j-m)*M]==5 || wind[p*K+(j-m)*M]==7|| wind[p*K+(j-m)*M]==6){
                                UFM[p*K+p+j*MM]=Ucoarse[p*K+(j-m)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[p*K+(j-m)*M];
                            }
                            else{
                                UFM[p*K+p+j*MM]=1000;
                            }
                        }
                        else if (i==(p+1)*K+p){
                            if (wind[(p+1)*K+(j-m)*M]==1 || wind[(p+1)*K+(j-m)*M]==4 || wind[(p+1)*K+(j-m)*M]==5 || wind[(p+1)*K+(j-m)*M]==7 || wind[(p+1)*K+(j-m)*M]==8){
                                UFM[(p+1)*K+p+(j)*MM]=Ucoarse[(p+1)*K+(j-m)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[(p+1)*K+(j-m)*M];
                            }
                            else{
                                UFM[(p+1)*K+p+(j)*MM]=1000;
                            }
                        }
                        else if (j==m*K+m){
                            if (wind[i-p+(m*K)*M]==1 || wind[i-p+(m*K)*M]==2 || wind[i-p+(m*K)*M]==6 || wind[i-p+(m*K)*M]==8 || wind[i-p+(m*K)*M]==5){
                                UFM[i+(m*K+m)*MM]=Ucoarse[i-p+(m*K)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[i-p+(m*K)*M];
                            }
                            else{
                                UFM[i+(m*K+m)*MM]=1000;
                            }
                        }
                        else if (j==(m+1)*K+m){
                            if (wind[i-p+((m+1)*K)*M]==3 || wind[i-p+((m+1)*K)*M]==4 || wind[i-p+((m+1)*K)*M]==6 || wind[i-p+((m+1)*K)*M]==8 || wind[i-p+((m+1)*K)*M]==7){
                                UFM[i+((m+1)*K+m)*MM]=Ucoarse[i-p+((m+1)*K)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[i-p+((m+1)*K)*M];
                            }
                            else{
                                UFM[i+((m+1)*K+m)*MM]=1000;
                            }

                        }

                    }

                }
            }
        }
    }





    //fast sweep
    //s1=-1;
    //s2=1;
    for (m=0; m<=N-1; m++){ //should be parallel for loop
        for (p=0; p<=N-1; p++){ //should be parallel for loop
            for (int iter=0; iter<4;iter++){
               for (s1=1; s1>=-1; s1-=2){
                 for (s2=1; s2>=-1; s2-=2){
                 //s1=-1;
                 //s2=1;
                        for (i=(s1<0 ? (p+1)*K+p : p*K+p); (s1<0 ? i>=p*K+p : i<= (p+1)*K+p); i+=s1){
                            for (j=(s2<0 ? (m+1)*K+m : m*K+m); (s2<0 ? j>=m*K+m : j<= (m+1)*K+m); j+=s2)
                            {
                                if (obs[(i-p)+(j-m)*M]>1)
                                {
                                    //don't do anything; point is in obstacle
                                }
                                else if (BC[(i-p)+(j-m)*M]>=0)
                                {
                                    UFM[i+j*MM]=BC[(i-p)+(j-m)*M];
                                }
                                //else if (BDRY[i+j*MM]>=0){
                                  // UFM[i+j*MM]=BDRY[i+j*MM];
                                //}
                                /*else if (i==(p+1)*K+p && j==(m+1)*K+m && wind[(i-p)+(j-m)*M]==4){
                                  UFM[i+j*MM]=1000;
                                }
                                else if (i==(p+1)*K+p && j==m*K+m && wind[(i-p)+(j-m)*M]==1){
                                  UFM[i+j*MM]=1000;
                                }
                                else if (i==p*K+p && j==(m+1)*K+m && wind[(i-p)+(j-m)*M]==3){
                                  UFM[i+j*MM]=1000;
                                }
                                else if (i==p*K+p && j==m*K+m && wind[(i-p)+(j-m)*M]==2){
                                  UFM[i+j*MM]=1000;
                                }*/
                                else if (BC[(i-p)+(j-m)*M]<0)
                                {

                                    if (j==m*K+m){
                                        m1=(s2<0 ? UFM[i+(j+1)*MM] : 1000);//UFM[i+(j+1)*MM];//
                                        tempwindy = (s2<0 ? 7 : 5);;
                                        //tempwindy=7;
                                    }
                                    else if (j==(m+1)*K+m){
                                        m1=(s2<0 ? 1000 : UFM[i+(j-1)*MM]);//UFM[i+(j-1)*MM];//
                                        tempwindy = (s2<0 ? 7 : 5);
                                        //tempwindy=5;
                                    }
                                    else{
                                        m1=(s2>0 ? UFM[i+(j-1)*MM] : UFM[i+(j+1)*MM]);//min(UFM[i+(j-1)*MM],UFM[i+(j+1)*MM]);//
                                        if (s2>0){//(UFM[i+(j-1)*MM]<UFM[i+(j+1)*MM])//
                                            tempwindy=5;
                                        }
                                        else{
                                            tempwindy=7;
                                        }
                                    }
                                    if (i==p*K+p){
                                        m2=(s1<0 ? UFM[(i+1)+j*MM] : 1000);//UFM[(i+1)+j*MM];//
                                        tempwindx = (s1<0 ? 8 : 6);
                                        //tempwindx=8;
                                    }
                                    else if (i==(p+1)*K+p){
                                        m2=(s1<0 ? 1000 : UFM[(i-1)+j*MM]);//UFM[(i-1)+j*MM];// ;
                                        tempwindx = (s1<0 ? 8 : 6);
                                        //tempwindx=6;
                                    }
                                    else{
                                        m2=(s1>0 ? UFM[(i-1)+j*MM] : UFM[(i+1)+j*MM]);//min(UFM[(i-1)+j*MM], UFM[(i+1)+j*MM]);//(s1>0 ? UFM[(i-1)+j*MM] : UFM[(i+1)+j*MM]);
                                        if (s1>0){//(UFM[(i-1)+j*MM]<UFM[(i+1)+j*MM])
                                            tempwindx=6;
                                        }
                                        else{
                                            tempwindx=8;
                                        }
                                    }

                                    //discr=2*pow(f[(i-p)+(j-m)*M],2)*pow(h,2)-pow(m1-m2,2);
                                    discr=2*pow(f[(i-p)+(j-m)*M],2)*pow(h,2)-pow(m1-m2,2);
                                    if (abs(m1-m2)<f[(i-p)+(j-m)*M]*h)
                                    {
                                        unew=0.5*(m1+m2+sqrt(discr));

                                        if (tempwindx==6 && tempwindy==5){
                                            windtemp=2;
                                        }
                                        else if (tempwindx==6 && tempwindy==7){
                                            windtemp=3;
                                        }
                                        else if (tempwindx==8 && tempwindy==5){
                                            windtemp=1;
                                        }
                                        else if (tempwindx==8 && tempwindy==7){
                                            windtemp=4;
                                        }



                                        old=UFM[i+j*MM];
                                        if (old<unew){
                                            UFM[i+j*MM]=old;

                                        }
                                        else{

                                            WIND[i+j*MM]=windtemp;

                                            UFM[i+j*MM]=unew;
                                        }
                                    }
                                    else
                                    {
                                        unew=min(m1,m2)+f[(i-p)+(j-m)*M]*h;

                                        if (m1<m2){
                                            windtemp=tempwindy;
                                        }
                                        else{
                                            windtemp=tempwindx;
                                        }



                                        /*if (i==p*K+p){
                                          if (wind[(i-p)+(j-m)*M]==2 && windtemp==5){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==2 && windtemp==7){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==3 && windtemp==5){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==3 && windtemp==7){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==3 && windtemp==6){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==3 && windtemp==8){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==2 && windtemp==6){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==2 && windtemp==8){
                                            continue;
                                          }
                                        }
                                        else if (i==(p+1)*K+p){
                                          if (wind[(i-p)+(j-m)*M]==1 && windtemp==5){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==1 && windtemp==7){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==4 && windtemp==5){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==4 && windtemp==7){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==4 && windtemp==6){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==4 && windtemp==8){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==1 && windtemp==6){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==1 && windtemp==8){
                                            continue;
                                          }
                                        }

                                        if (j==m*K+m){
                                          if (wind[(i-p)+(j-m)*M]==1 && windtemp==8){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==1 && windtemp==6){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==1 && windtemp==7){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==1 && windtemp==5){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==2 && windtemp==8){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==2 && windtemp==7){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==2 && windtemp==6){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==2 && windtemp==5){
                                            continue;
                                          }
                                        }
                                        else if (j==(m+1)*K+m){
                                          if (wind[(i-p)+(j-m)*M]==4 && windtemp==8){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==4 && windtemp==6){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==4 && windtemp==7){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==4 && windtemp==5){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==3 && windtemp==8){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==3 && windtemp==7){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==3 && windtemp==6){
                                            continue;
                                          }
                                          else if (wind[(i-p)+(j-m)*M]==3 && windtemp==5){
                                            continue;
                                          }
                                        }*/
                                        old=UFM[i+j*MM];
                                        if (old<unew){
                                            UFM[i+j*MM]=old;
                                        }
                                        else{
                                            WIND[i+j*MM]=windtemp;
                                            UFM[i+j*MM]=unew;
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





    //now take minimums along shared boundaries

double temp1,temp2,temp3,temp4,temp;
double * vec1, * vec2, *vec3, *vec4, *vec5;
vec1=(double *) mxCalloc(2, sizeof(double));
vec2=(double *) mxCalloc(2, sizeof(double));
vec3=(double *) mxCalloc(2, sizeof(double));
vec4=(double *) mxCalloc(2, sizeof(double));
vec5=(double *) mxCalloc(2, sizeof(double));
//nonshifted
    for (i=0; i<=M-1; i+=K){
      for (j=0; j<=M-1; j+=K){
        if (BC[i+j*M]>=0){
          FINE[i+j*M]=BC[i+j*M];
        }
        else{
          if (i==0){
            if (j==0){
              fine[0]=UFM[0];
              FINE[0]=fine[0];
              fwind[0]=WIND[0];
            }
            else if (j==M-1) {
              fine[i+j*M]=UFM[(MM-1)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[(MM-1)*MM];
            }
            else{
              temp1=UFM[(j+j/K-1)*MM];
              temp2=UFM[(j+j/K)*MM];
              whichvec(WIND[(j+j/K-1)*MM],vec1);
              whichvec(WIND[(j+j/K)*MM],vec2);
              whichvec(wind[i+j*M],vec3);

              if (vec1[0]>=0 && vec2[0]>=0 && vec3[0]>=0){
                fine[i+j*M]=UFM[(j+j/K-1)*MM];
                FINE[i+j*M]=fine[i+j*M];
                fwind[i+j*M]=WIND[(j+j/K-1)*MM];
              }
              else if (-1*vec1[0]>=0 && -1*vec2[0]>=0 && -1*vec3[0]>=0){
                fine[i+j*M]=UFM[(j+j/K)*MM];
                FINE[i+j*M]=fine[i+j*M];
                fwind[i+j*M]=WIND[(j+j/K)*MM];
              }
              else{
                if (temp1<=temp2){
                  fine[i+j*M]=UFM[(j+j/K-1)*MM];
                  FINE[i+j*M]=fine[i+j*M];
                  fwind[i+j*M]=WIND[(j+j/K-1)*MM];
                }
                else{
                  fine[i+j*M]=UFM[(j+j/K)*MM];
                  FINE[i+j*M]=fine[i+j*M];
                  fwind[i+j*M]=WIND[(j+j/K)*MM];
                }
              }
            }
          }
          else if(i==M-1){
            if (j==0){
              fine[i+j*M]=UFM[MM-1];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[MM-1];
            }
            else if (j==M-1) {
              fine[i+j*M]=UFM[MM-1+(MM-1)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[MM-1+(MM-1)*MM];

            }
            else{
              temp1=UFM[(j+j/K-1)*MM];
              temp2=UFM[(j+j/K)*MM];
              whichvec(WIND[(j+j/K-1)*MM],vec1);
              whichvec(WIND[(j+j/K)*MM],vec2);
              whichvec(wind[i+j*M],vec3);

              if (vec1[0]>=0 && vec2[0]>=0 && vec3[0]>=0){
                fine[i+j*M]=UFM[(j+j/K-1)*MM];
                FINE[i+j*M]=fine[i+j*M];
                fwind[i+j*M]=WIND[(j+j/K-1)*MM];
              }
              else if (-1*vec1[0]>=0 && -1*vec2[0]>=0 && -1*vec3[0]>=0){
                fine[i+j*M]=UFM[(j+j/K)*MM];
                FINE[i+j*M]=fine[i+j*M];
                fwind[i+j*M]=WIND[(j+j/K)*MM];
              }
              else{
                if (temp1<=temp2){
                  fine[i+j*M]=UFM[(j+j/K-1)*MM];
                  FINE[i+j*M]=fine[i+j*M];
                  fwind[i+j*M]=WIND[(j+j/K-1)*MM];
                }
                else{
                  fine[i+j*M]=UFM[(j+j/K)*MM];
                  FINE[i+j*M]=fine[i+j*M];
                  fwind[i+j*M]=WIND[(j+j/K)*MM];
                }
              }
            }
          }
          else if (j==0){
            temp1=UFM[i+i/K-1];
            temp2=UFM[i+i/K];
            whichvec(WIND[i+i/K-1],vec1);
            whichvec(WIND[i+i/K],vec2);
            whichvec(wind[i+j*M],vec3);

            if (vec1[1]>=0 && vec2[1]>=0 && vec3[1]>=0){
              fine[i+j*M]=UFM[i+i/K-1];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[i+i/K-1];
            }
            else if (-1*vec1[1]>=0 && -1*vec2[1]>=0 && -1*vec3[1]>=0){
              fine[i+j*M]=UFM[i+i/K];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[i+i/K];
            }
            else{
              if (temp1<=temp2){
                fine[i+j*M]=UFM[i+i/K-1];
                FINE[i+j*M]=fine[i+j*M];
                fwind[i+j*M]=WIND[i+i/K-1];
              }
              else{
                fine[i+j*M]=UFM[i+i/K];
                FINE[i+j*M]=fine[i+j*M];
                fwind[i+j*M]=WIND[i+i/K];
              }
            }
          }
          else if (j==M-1){
            temp1=UFM[i+i/K-1+(MM-1)*MM];
            temp2=UFM[i+i/K+(MM-1)*MM];
            whichvec(WIND[i+i/K-1+(MM-1)*MM],vec1);
            whichvec(WIND[i+i/K+(MM-1)*MM],vec2);
            whichvec(wind[i+j*M],vec3);

            if (vec1[1]>=0 && vec2[1]>=0 && vec3[1]>=0){
              fine[i+j*M]=UFM[i+i/K-1+(MM-1)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[i+i/K-1+(MM-1)*MM];
            }
            else if (-1*vec1[1]>=0 && -1*vec2[1]>=0 && -1*vec3[1]>=0){
              fine[i+j*M]=UFM[i+i/K+(MM-1)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[i+i/K+(MM-1)*MM];
            }
            else{
              if (temp1<=temp2){
                fine[i+j*M]=UFM[i+i/K-1+(MM-1)*MM];
                FINE[i+j*M]=fine[i+j*M];
                fwind[i+j*M]=WIND[i+i/K-1+(MM-1)*MM];
              }
              else{
                fine[i+j*M]=UFM[i+i/K+(MM-1)*MM];
                FINE[i+j*M]=fine[i+j*M];
                fwind[i+j*M]=WIND[i+i/K+(MM-1)*MM];
              }
            }
          }
        else{
          temp1=UFM[i+i/K-1+(j+j/K-1)*MM];
          temp2=UFM[i+i/K-1+(j+j/K)*MM];
          temp3=UFM[i+i/K+(j+j/K-1)*MM];
          temp4=UFM[i+i/K+(j+j/K)*MM];
          whichvec(WIND[i+i/K-1+(j+j/K-1)*MM],vec1);
          whichvec(WIND[i+i/K-1+(j+j/K)*MM],vec2);
          whichvec(WIND[i+i/K+(j+j/K-1)*MM],vec3);
          whichvec(WIND[i+i/K+(j+j/K)*MM],vec4);
          whichvec(wind[i+j*M],vec5);
          if ((vec1[0]-vec1[1]>0) && (vec2[0]-vec2[1]>0) && (vec3[0]-vec3[1]>0) && (vec4[0]-vec4[1]>0) && (vec5[0]-vec5[1]>0)){
            fine[i+j*M]=UFM[i+i/K-1+(j+j/K-1)*MM];
            FINE[i+j*M]=fine[i+j*M];
            fwind[i+j*M]=WIND[i+i/K-1+(j+j/K-1)*MM];
          }
          else if ((-1*vec1[0]-vec1[1]>0) && (-1*vec2[0]-vec2[1]>0) && (-1*vec3[0]-vec3[1]>0) && (-1*vec4[0]-vec4[1]>0) && (-1*vec5[0]-vec5[1]>0)){
            fine[i+j*M]=UFM[i+i/K-1+(j+j/K)*MM];
            FINE[i+j*M]=fine[i+j*M];
            fwind[i+j*M]=WIND[i+i/K-1+(j+j/K)*MM];
          }
          else if ((vec1[0]+vec1[1]>0) && (vec2[0]+vec2[1]>0) && (vec3[0]+vec3[1]>0) && (vec4[0]+vec4[1]>0) && (vec5[0]+vec5[1]>0)){
            fine[i+j*M]=UFM[i+i/K+(j+j/K-1)*MM];
            FINE[i+j*M]=fine[i+j*M];
            fwind[i+j*M]=WIND[i+i/K+(j+j/K-1)*MM];
          }
          else if ((-1*vec1[0]+vec1[1]>0) && (-1*vec2[0]+vec2[1]>0) && (-1*vec3[0]+vec3[1]>0) && (-1*vec4[0]+vec4[1]>0) && (-1*vec5[0]+vec5[1]>0)){
            fine[i+j*M]=UFM[i+i/K+(j+j/K)*MM];
            FINE[i+j*M]=fine[i+j*M];
            fwind[i+j*M]=WIND[i+i/K+(j+j/K)*MM];
          }
          else{
            temp=min({temp1,temp2,temp3,temp4});
            if (temp==temp1){
              fine[i+j*M]=UFM[i+i/K-1+(j+j/K-1)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[i+i/K-1+(j+j/K-1)*MM];
            }
            else if (temp==temp2){
              fine[i+j*M]=UFM[i+i/K-1+(j+j/K)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[i+i/K-1+(j+j/K)*MM];
            }
            else if (temp==temp3){
              fine[i+j*M]=UFM[i+i/K+(j+j/K-1)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[i+i/K+(j+j/K-1)*MM];
            }
            else{
              fine[i+j*M]=UFM[i+i/K+(j+j/K)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[i+i/K+(j+j/K)*MM];
            }
          }
        }
      }
    }
  }
//horizontal shifts
  for (p=1;p<K;p++){
    for (i=0;i<=M-1;i+=K){
      for (j=p;j<=M-1;j+=K){
        if (BC[i+j*M]>=0){
          FINE[i+j*M]=BC[i+j*M];
        }
        else{
          if (i==0){
            fine[i+j*M]=UFM[(j+(j-p)/K)*MM];
            FINE[i+j*M]=fine[i+j*M];
            fwind[i+j*M]=WIND[(j+(j-p)/K)*MM];
          }
          else if (i==M-1){
            fine[i+j*M]=UFM[MM-1+(j+(j-p)/K)*MM];
            FINE[i+j*M]=fine[i+j*M];
            fwind[i+j*M]=WIND[MM-1+(j+(j-p)/K)*MM];
          }
          else{
            temp1=UFM[i+i/K-1+(j+(j-p)/K)*MM];
            temp2=UFM[i+i/K+(j+(j-p)/K)*MM];
            whichvec(WIND[i+i/K-1+(j+(j-p)/K)*MM],vec1);
            whichvec(WIND[i+i/K+(j+(j-p)/K)*MM],vec2);
            whichvec(wind[i+j*M],vec3);

            if (vec1[1]>=0 && vec2[1]>=0 && vec3[1]>=0){
              fine[i+j*M]=UFM[i+i/K-1+(j+(j-p)/K)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[i+i/K-1+(j+(j-p)/K)*MM];
            }
            else if (-1*vec1[1]>=0 && -1*vec2[1]>=0 && -1*vec3[1]>=0){
              fine[i+j*M]=UFM[i+i/K+(j+(j-p)/K)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[i+i/K+(j+(j-p)/K)*MM];
            }
            else{
              if (temp1<=temp2){
                fine[i+j*M]=UFM[i+i/K-1+(j+(j-p)/K)*MM];
                FINE[i+j*M]=fine[i+j*M];
                fwind[i+j*M]=WIND[i+i/K-1+(j+(j-p)/K)*MM];
              }
              else{
                fine[i+j*M]=UFM[i+i/K+(j+(j-p)/K)*MM];
                FINE[i+j*M]=fine[i+j*M];
                fwind[i+j*M]=WIND[i+i/K+(j+(j-p)/K)*MM];
              }
            }
          }
        }
      }
    }
  }

//vertical shifts
for (p=1;p<K;p++){
  for (i=p;i<=M-1;i+=K){
    for (j=0;j<=M-1;j+=K){
      if (BC[i+j*M]>=0){
        FINE[i+j*M]=BC[i+j*M];
      }
      else{
        if (j==0){
          fine[i+j*M]=UFM[i+(i-p)/K];
          FINE[i+j*M]=fine[i+j*M];
          fwind[i+j*M]=WIND[i+(i-p)/K];
        }
        else if (j==M-1){
          fine[i+j*M]=UFM[i+(i-p)/K+(MM-1)*MM];
          FINE[i+j*M]=fine[i+j*M];
          fwind[i+j*M]=WIND[i+(i-p)/K+(MM-1)*MM];
        }
        else{
          temp1=UFM[(i+(i-p)/K)+(j+j/K-1)*MM];
          temp2=UFM[(i+(i-p)/K)+(j+j/K)*MM];
          whichvec(WIND[(i+(i-p)/K)+(j+j/K-1)*MM],vec1);
          whichvec(WIND[(i+(i-p)/K)+(j+j/K)*MM],vec2);
          whichvec(wind[i+j*M],vec3);

          if (vec1[0]>=0 && vec2[0]>=0 && vec3[0]>=0){
            fine[i+j*M]=UFM[(i+(i-p)/K)+(j+j/K-1)*MM];
            FINE[i+j*M]=fine[i+j*M];
            fwind[i+j*M]=WIND[(i+(i-p)/K)+(j+j/K-1)*MM];
          }
          else if (-1*vec1[0]>=0 && -1*vec2[0]>=0 && -1*vec3[0]>=0){
            fine[i+j*M]=UFM[(i+(i-p)/K)+(j+j/K)*MM];
            FINE[i+j*M]=fine[i+j*M];
            fwind[i+j*M]=WIND[(i+(i-p)/K)+(j+j/K)*MM];
          }
          else{
            if (temp1<=temp2){
              fine[i+j*M]=UFM[(i+(i-p)/K)+(j+j/K-1)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[(i+(i-p)/K)+(j+j/K-1)*MM];
            }
            else{
              fine[i+j*M]=UFM[(i+(i-p)/K)+(j+j/K)*MM];
              FINE[i+j*M]=fine[i+j*M];
              fwind[i+j*M]=WIND[(i+(i-p)/K)+(j+j/K)*MM];
            }
          }
        }
      }
    }
  }
}

//fill in subdomains

for (m=0;m<=N-1;m++){
  for (p=0;p<=N-1;p++){
    for (i=p*K+p+1; i<(p+1)*K+p; i++){
      for (j=m*K+m+1; j<(m+1)*K+m; j++){
        fine[(i-p)+(j-m)*M]=UFM[i+j*MM];
        FINE[(i-p)+(j-m)*M]=fine[(i-p)+(j-m)*M];
        fwind[(i-p)+(j-m)*M]=WIND[i+j*MM];
      }
    }
  }
}









    mxFree(WIND);
    mxFree(BDRY);
    mxFree(fine);
    mxFree(UFM);
}

#endif /* fineupdate_h */
