//
//  fineupdate.h
//  Eikonal Project 5
//
//  Created by Lindsay Martin on 6/22/16.
//  Copyright © 2016 Lindsay Martin. All rights reserved.
//

#ifndef fineupdate_h
#define fineupdate_h
#include <algorithm>
#include <cmath>
using namespace std;

void fineupdate(double* FINE, double *Ucoarse, double *f, double *BC,double *obs, int N,int K, double *wind, double *fwind);
void fineupdate(double * FINE, double *Ucoarse, double *f,double *BC,double *obs, int N,int K, double *wind, double *fwind)
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
                            if (wind[p*K+(j-m)*M]==2 || wind[p*K+(j-m)*M]==3 || wind[p*K+(j-m)*M]==6){// || wind[p*K+(j-m)*M]==5 || wind[p*K+(j-m)*M]==7 ){
                                UFM[p*K+p+j*MM]=Ucoarse[p*K+(j-m)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[p*K+(j-m)*M];
                            }
                            else{
                                UFM[p*K+p+j*MM]=1000;
                            }
                        }
                        else if (i==(p+1)*K+p){
                            if (wind[(p+1)*K+(j-m)*M]==1 || wind[(p+1)*K+(j-m)*M]==4 || wind[(p+1)*K+(j-m)*M]==8){// || wind[(p+1)*K+(j-m)*M]==7 || wind[(p+1)*K+(j-m)*M]==5){
                                UFM[(p+1)*K+p+(j)*MM]=Ucoarse[(p+1)*K+(j-m)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[(p+1)*K+(j-m)*M];
                            }
                            else{
                                UFM[(p+1)*K+p+(j)*MM]=1000;
                            }
                        }
                        else if (j==m*K+m){
                            if (wind[i-p+(m*K)*M]==1 || wind[i-p+(m*K)*M]==2 || wind[i-p+(m*K)*M]==5){// || wind[i-p+(m*K)*M]==8 || wind[i-p+(m*K)*M]==6){
                                UFM[i+(m*K+m)*MM]=Ucoarse[i-p+(m*K)*M];
                                WIND[i+j*MM]=wind[(i-p)+(j-m)*M];
                                BDRY[i+j*MM]=Ucoarse[i-p+(m*K)*M];
                            }
                            else{
                                UFM[i+(m*K+m)*MM]=1000;
                            }
                        }
                        else if (j==(m+1)*K+m){
                            if (wind[i-p+((m+1)*K)*M]==3 || wind[i-p+((m+1)*K)*M]==4 || wind[i-p+((m+1)*K)*M]==7){// || wind[i-p+((m+1)*K)*M]==8 || wind[i-p+((m+1)*K)*M]==6){
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

/*for (m=0; m<=N-1; m++){ //should be parallel for loop
  for (p=0; p<=N-1; p++){ //should be parallel for loop
    for (s1=1; s1>=-1; s1-=2){
      for (s2=1; s2>=-1; s2-=2){
         for (i=(s1<0 ? (p+1)*K+p : p*K+p); (s1<0 ? i>=p*K+p : i<= (p+1)*K+p); i+=s1){
            for (j=(s2<0 ? (m+1)*K+m : m*K+m); (s2<0 ? j>=m*K+m : j<= (m+1)*K+m); j+=s2){
              if (i==p*K+p && j==m*K+m){
                if(BDRY[i+j*MM]<500){
                  if (wind[(i-p)+(j-m)*M]==1){
                    if (UFM[i+j*MM]<UFM[i+1+j*MM]){
                      UFM[i+j*MM]=UFM[i+1+j*MM];
                    }
                  }
                  else if (wind[(i-p)+(j-m)*M]==3){
                    if (UFM[i+j*MM]<UFM[i+(j+1)*MM]){
                      UFM[i+j*MM]=UFM[i+(j+1)*MM];
                    }
                  }
                }
              }
              else if (i==p*K+p && j==(m+1)*K+m){
                if(BDRY[i+j*MM]<500){
                  if (wind[(i-p)+(j-m)*M]==4){
                    if (UFM[i+j*MM]<UFM[i+1+j*MM]){
                      UFM[i+j*MM]=UFM[i+1+j*MM];
                    }
                  }
                  else if (wind[(i-p)+(j-m)*M]==2){
                    if (UFM[i+j*MM]<UFM[i+(j-1)*MM]){
                      UFM[i+j*MM]=UFM[i+(j-1)*MM];
                    }
                  }
                }
              }
              else if (i==(p+1)*K+p && j==m*K+m){
                if(BDRY[i+j*MM]<500){
                  if (wind[(i-p)+(j-m)*M]==2){
                    if (UFM[i+j*MM]<UFM[i-1+j*MM]){
                      UFM[i+j*MM]=UFM[i-1+j*MM];
                    }
                  }
                  else if (wind[(i-p)+(j-m)*M]==4){
                    if (UFM[i+j*MM]<UFM[i+(j+1)*MM]){
                      UFM[i+j*MM]=UFM[i+(j+1)*MM];
                    }
                  }
                }
              }
              else if (i==(p+1)*K+p && j==(m+1)*K+m){
                if(BDRY[i+j*MM]<500){
                  if (wind[(i-p)+(j-m)*M]==1){
                    if (UFM[i+j*MM]<UFM[i+(j-1)*MM]){
                      UFM[i+j*MM]=UFM[i+(j-1)*MM];
                    }
                  }
                  else if (wind[(i-p)+(j-m)*M]==3){
                    if (UFM[i+j*MM]<UFM[i-1+j*MM]){
                      UFM[i+j*MM]=UFM[i-1+j*MM];
                    }
                  }
                }

              }
              else if (i==p*K+p){
                if(BDRY[i+j*MM]<500){
                  if (wind[(i-p)+(j-m)*M]==5 || wind[(i-p)+(j-m)*M]==2){
                    if (UFM[i+j*MM]<UFM[i+(j-1)*MM]){
                      UFM[i+j*MM]=UFM[i+(j-1)*MM];
                    }
                  }
                  else if (wind[(i-p)+(j-m)*M]==3 || wind[(i-p)+(j-m)*M]==7){
                    if (UFM[i+j*MM]<UFM[i+(j+1)*MM]){
                      UFM[i+j*MM]=UFM[i+(j+1)*MM];
                    }
                  }
                }

              }
              else if (i==(p+1)*K+p){
                if(BDRY[i+j*MM]<500){
                  if (wind[(i-p)+(j-m)*M]==5 || wind[(i-p)+(j-m)*M]==1){
                    if (UFM[i+j*MM]<UFM[i+(j-1)*MM]){
                      UFM[i+j*MM]=UFM[i+(j-1)*MM];
                    }
                  }
                  else if (wind[(i-p)+(j-m)*M]==4 || wind[(i-p)+(j-m)*M]==7){
                    if (UFM[i+j*MM]<UFM[i+(j+1)*MM]){
                      UFM[i+j*MM]=UFM[i+(j+1)*MM];
                    }
                  }
                }

              }
              else if (j==m*K+m){
                if(BDRY[i+j*MM]<500){
                  if (wind[(i-p)+(j-m)*M]==8 || wind[(i-p)+(j-m)*M]==1){
                    if (UFM[i+j*MM]<UFM[i+1+j*MM]){
                      UFM[i+j*MM]=UFM[i+1+j*MM];
                    }
                  }
                  else if (wind[(i-p)+(j-m)*M]==6 || wind[(i-p)+(j-m)*M]==2){
                    if (UFM[i+j*MM]<UFM[i-1+j*MM]){
                      UFM[i+j*MM]=UFM[i-1+j*MM];
                    }
                  }
                }
              }
              else if (j==(m+1)*K+m){
                if(BDRY[i+j*MM]<500){
                  if (wind[(i-p)+(j-m)*M]==8 || wind[(i-p)+(j-m)*M]==4){
                    if (UFM[i+j*MM]<UFM[i+1+j*MM]){
                      UFM[i+j*MM]=UFM[i+1+j*MM];
                    }
                  }
                  else if (wind[(i-p)+(j-m)*M]==6 || wind[(i-p)+(j-m)*M]==3){
                    if (UFM[i+j*MM]<UFM[i-1+j*MM]){
                      UFM[i+j*MM]=UFM[i-1+j*MM];
                    }
                  }
                }

              }
          }
        }
      }
    }
  }
}*/

    //fast sweep
    //s1=-1;
    //s2=1;
    for (m=0; m<=N-1; m++){ //should be parallel for loop
        for (p=0; p<=N-1; p++){ //should be parallel for loop
            for (int iter=0; iter<10;iter++){
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
                                    if (j==0){
                                      m1=UFM[i+(j+1)*MM];
                                      tempwindy=7;
                                    }
                                    else if (j==MM-1){
                                        m1=UFM[i+(j-1)*MM];
                                        tempwindy=5;
                                    }
                                    else if (j==m*K+m){
                                        m1=UFM[i+(j+1)*MM];//(s2<0 ? UFM[i+(j+1)*MM] : 1000);//UFM[i+(j+1)*MM];//
                                        //tempwindy = (s2<0 ? 7 : 5);;
                                        tempwindy=7;
                                    }
                                    else if (j==(m+1)*K+m){
                                        m1=UFM[i+(j-1)*MM];//(s2<0 ? 1000 : UFM[i+(j-1)*MM]);//
                                        //tempwindy = (s2<0 ? 7 : 5);
                                        tempwindy=5;
                                    }
                                    else{
                                        m1=min(UFM[i+(j-1)*MM],UFM[i+(j+1)*MM]);//(s2>0 ? UFM[i+(j-1)*MM] : UFM[i+(j+1)*MM]);//
                                        if (UFM[i+(j-1)*MM]<UFM[i+(j+1)*MM]){//(s2>0){//(UFM[i+(j-1)*MM]<UFM[i+(j+1)*MM])//
                                            tempwindy=5;
                                        }
                                        else{
                                            tempwindy=7;
                                        }
                                    }
                                    if (i==0){
                                      m2=UFM[(i+1)+j*MM];
                                      tempwindx=8;
                                    }
                                    else if (i==MM-1){
                                      m2=UFM[(i-1)+j*MM];
                                      tempwindx=6;
                                    }
                                    else if (i==p*K+p){
                                        m2=UFM[(i+1)+j*MM];//(s1<0 ? UFM[(i+1)+j*MM] : 1000);//
                                        //tempwindx = (s1<0 ? 8 : 6);
                                        tempwindx=8;
                                    }
                                    else if (i==(p+1)*K+p){
                                        m2=UFM[(i-1)+j*MM];// ;(s1<0 ? 1000 : UFM[(i-1)+j*MM]);//
                                        //tempwindx = (s1<0 ? 8 : 6);
                                        tempwindx=6;
                                    }
                                    else{
                                        m2=min(UFM[(i-1)+j*MM], UFM[(i+1)+j*MM]);//(s1>0 ? UFM[(i-1)+j*MM] : UFM[(i+1)+j*MM]);//
                                        if (UFM[(i-1)+j*MM]<UFM[(i+1)+j*MM]){//(s1>0){//
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


    for (m=0; m<=N-1; m++){
        for (p=0; p<=N-1; p++){
            for (i=p*K+p; i<=(p+1)*K+p; i++){
                for (j=m*K+m; j<=(m+1)*K+m; j++)
                {
                    if (BC[(i-p)+(j-m)*M]>=0)
                    {
                        FINE[(i-p)+(j-m)*M]=BC[(i-p)+(j-m)*M];
                    }
                    else{
                        //if (BDRY[i+j*MM]<0){
                        old=fine[(i-p)+(j-m)*M];
                        cont=1;
                        if (i==p*K+p && j==m*K+m){
                            if (WIND[i+j*MM]==7 || WIND[i+j*MM]==8 || WIND[i+j*MM]==4)
                                if (UFM[i+j*MM]<old){
                                    fine[(i-p)+(j-m)*M]=UFM[i+j*MM];
                                    FINE[(i-p)+(j-m)*M]=fine[(i-p)+(j-m)*M];
                                    fwind[(i-p)+(j-m)*M]=WIND[i+j*MM];
                                }
                            cont=0;

                        }
                        else if (i==p*K+p && j==(m+1)*K+m){
                            if (WIND[i+j*MM]==8 || WIND[i+j*MM]==1 || WIND[i+j*MM]==5)
                                if (UFM[i+j*MM]<old){
                                    fine[(i-p)+(j-m)*M]=UFM[i+j*MM];
                                    FINE[(i-p)+(j-m)*M]=fine[(i-p)+(j-m)*M];
                                    fwind[(i-p)+(j-m)*M]=WIND[i+j*MM];
                                }
                            cont=0;

                        }

                        else if (i==(p+1)*K+p && j==m*K+m){
                            if (WIND[i+j*MM]==6 || WIND[i+j*MM]==3 || WIND[i+j*MM]==7)
                                if (UFM[i+j*MM]<old){
                                     fine[(i-p)+(j-m)*M]=UFM[i+j*MM];
                                    FINE[(i-p)+(j-m)*M]=fine[(i-p)+(j-m)*M];
                                    fwind[(i-p)+(j-m)*M]=WIND[i+j*MM];
                                }
                            cont=0;

                        }
                        else if (i==(p+1)*K+p && j==(m+1)*K+m){
                            if (WIND[i+j*MM]==5 || WIND[i+j*MM]==2 || WIND[i+j*MM]==6)
                                if (UFM[i+j*MM]<old){
                                    fine[(i-p)+(j-m)*M]=UFM[i+j*MM];
                                    FINE[(i-p)+(j-m)*M]=fine[(i-p)+(j-m)*M];
                                    fwind[(i-p)+(j-m)*M]=WIND[i+j*MM];
                                }
                            cont=0;

                        }
                        else if (i==p*K+p){
                            if (WIND[i+j*MM]==7 || WIND[i+j*MM]==4 || WIND[i+j*MM]==8 || WIND[i+j*MM]==1 || WIND[i+j*MM]==5)
                                if (UFM[i+j*MM]<old){
                                    fine[(i-p)+(j-m)*M]=UFM[i+j*MM];
                                    FINE[(i-p)+(j-m)*M]=fine[(i-p)+(j-m)*M];
                                    fwind[(i-p)+(j-m)*M]=WIND[i+j*MM];
                                }
                            cont=0;

                        }
                        else if (i==(p+1)*K+p){
                            if (WIND[i+j*MM]==7 || WIND[i+j*MM]==3 || WIND[i+j*MM]==6 || WIND[i+j*MM]==2 || WIND[i+j*MM]==5)
                                if (UFM[i+j*MM]<old){
                                    fine[(i-p)+(j-m)*M]=UFM[i+j*MM];
                                    FINE[(i-p)+(j-m)*M]=fine[(i-p)+(j-m)*M];
                                    fwind[(i-p)+(j-m)*M]=WIND[i+j*MM];
                                }
                            cont=0;

                        }
                        else if (j==m*K+m){
                            if (WIND[i+j*MM]==8 || WIND[i+j*MM]==4 || WIND[i+j*MM]==7 || WIND[i+j*MM]==3 || WIND[i+j*MM]==6)
                                if (UFM[i+j*MM]<old){
                                    fine[(i-p)+(j-m)*M]=UFM[i+j*MM];
                                    FINE[(i-p)+(j-m)*M]=fine[(i-p)+(j-m)*M];
                                    fwind[(i-p)+(j-m)*M]=WIND[i+j*MM];
                                }
                            cont=0;

                        }
                        else if (j==(m+1)*K+m){
                            if (WIND[i+j*MM]==8 || WIND[i+j*MM]==1 || WIND[i+j*MM]==5 || WIND[i+j*MM]==2 || WIND[i+j*MM]==6)
                                if (UFM[i+j*MM]<old){
                                    fine[(i-p)+(j-m)*M]=UFM[i+j*MM];
                                    FINE[(i-p)+(j-m)*M]=fine[(i-p)+(j-m)*M];
                                    fwind[(i-p)+(j-m)*M]=WIND[i+j*MM];
                                }
                            cont=0;

                        }

                        if (cont==1){

                            if (UFM[i+j*MM]<old){
                                fine[(i-p)+(j-m)*M]=UFM[i+j*MM];
                                FINE[(i-p)+(j-m)*M]=fine[(i-p)+(j-m)*M];
                                fwind[(i-p)+(j-m)*M]=WIND[i+j*MM];
                            }
                        }
                    }
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
