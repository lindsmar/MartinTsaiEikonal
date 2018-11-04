//
//  fineupdate.h
//  Eikonal Project 5
//
//  Created by Lindsay Martin on 6/22/16.
//  Copyright © 2016 Lindsay Martin. All rights reserved.
//

/*#ifndef fineupdate_h
#define fineupdate_h
#include <algorithm>
#include <cmath>
using namespace std;
#define DIM1 2
#define DIM2 800
void fineupdate(double* FINE, double *Ucoarse[DIM1][DIM2], double *f, double *BC,double *obs, int N,int K, double *wind[DIM1][DIM2], double *fwind);*/
#include "sweepdefs.h"///scratch/04362/lmartin/paralleltests/coarsegridrevise/sweepdefs.h"

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


void fineupdate(double * FINE, double *Ucoarse[DIM1][DIM2], double *f,double *BC,double *obs, int N,int K, double *wind[DIM1][DIM2], double *fwind, int *totalsweeps)
{
    int M;
    int MM;
    double h;

    h=(double)1/(double) (N*K);
    M=N*K+1;
    MM=M+N-1;
   // int s1;
    //int s2;
    //int cont, tempwindx,tempwindy, windtemp;


    //double *FINE, *FINE2, *exwind, *cellbdry;
    double *fine, *WIND, *BDRY, *UFM;

    //fine2=(double *) mxCalloc(M*M, sizeof(double));
    fine=new double [M*M];
    //FINE=mxCreateDoubleMatrix(M,M,mxREAL);
    //fine=mxGetPr(FINE);
    //FINE2=mxCreateDoubleMatrix(M,M,mxREAL);
    //fine2=mxGetPr(FINE2);
    UFM=new double [MM*MM];
    WIND=new double [MM*MM];
    BDRY=new double [MM*MM];
    //int c;
    #pragma omp parallel for
    for (int c=0; c<(M+N-1)*(M+N-1); c++){
        UFM[c]=1000;
    }
    
#pragma omp parallel for
    for (int z=0; z<M*M; z++){
        fine[z]=1000;
        //FINE[z]=Ucoarse[z];
        
    }
    #pragma omp parallel for
    for (int z=0; z<(M+N-1)*(M+N-1); z++){
        WIND[z]=1000;
         BDRY[z]=-1;
    }
#pragma omp parallel for collapse(2)
    for (int dir=0; dir<2; dir++){
        for (int i=0; i<K; i++){
            for (int k=0; k<(N+1); k++){
                for (int m=0; m<(N+1); m++){
                    if (dir==0){
                        if (i==0){
                            fwind[k*K+m*K*M]=wind[dir][i][k+m*(N+1)];
                           //  FINE[k*K+m*K*M]=Ucoarse[dir][i][k+m*(N+1)];
                        }
                        else if (k!=14 && m!=14){
                            
                            fwind[k*K+(m*K+i)*M]=wind[dir][i][k+m*(N+1)];
                           // FINE[k*K+(m*K+i)*M]=Ucoarse[dir][i][k+m*(N+1)];
                        }
                    }
                    else{
                        if (i!=0 && k!=14 && m!=14){
                            fwind[k*K+i+m*K*M]=wind[dir][i][k+m*(N+1)];
                            //FINE[k*K+i+m*K*M]=Ucoarse[dir][i][k+m*(N+1)];
                        }
                    }
                }
            }
        }
    }
    
    
    
  
    //Set Bdry condition for each subdomain
  #pragma omp parallel for collapse(2)
    for (int m=0; m<=N-1; m++){
        for (int p=0; p<=N-1; p++){
            for (int i=p*K+p; i<=(p+1)*K+p; i++){
                for (int j=m*K+m; j<=(m+1)*K+m; j++)
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
                   
                    else{
                        if ((i-p)%K==0 || (j-m)%K==0){
                        int DIRkLIN[3];
                        int dir,k,lin;
                        cartto2CG(i-p,j-m,N+1,K,DIRkLIN);
                        dir=DIRkLIN[0];
                        k=DIRkLIN[1];
                        lin=DIRkLIN[2];
                        if (i==p*K+p && j==m*K+m){
                         
                                
                            if (wind[dir][k][lin]==5  ||  wind[dir][k][lin]==6 || wind[dir][k][lin]==3 || wind[dir][k][lin]==1 || wind[dir][k][lin]==2){
                                
                                
                                UFM[i+j*MM]=Ucoarse[dir][k][lin];
                                WIND[i+j*MM]=wind[dir][k][lin];
                                BDRY[i+j*MM]=Ucoarse[dir][k][lin];
                            }
                            else{
                                if (wind[dir][k][lin]!=0){
                                    UFM[i+j*MM]=1000;
                                }
                                
                            }
                        }
                        
                        else if (i==p*K+p && j==(m+1)*K+m){
                           
                            if (wind[dir][k][lin]==6  ||  wind[dir][k][lin]==7 || wind[dir][k][lin]==4 || wind[dir][k][lin]==2 || wind[dir][k][lin]==3){
                                UFM[i+j*MM]=Ucoarse[dir][k][lin];
                                WIND[i+j*MM]=wind[dir][k][lin];
                                BDRY[i+j*MM]=Ucoarse[dir][k][lin];
                            }
                            else{
                                if (wind[dir][k][lin]!=0){
                                    UFM[i+j*MM]=1000;
                                }
                            }
                        }
                        else if (i==(p+1)*K+p && j==m*K+m){
                           

                            if (wind[dir][k][lin]==8  ||  wind[dir][k][lin]==5 || wind[dir][k][lin]==4 || wind[dir][k][lin]==2|| wind[dir][k][lin]==1){
                                UFM[i+j*MM]=Ucoarse[dir][k][lin];
                                WIND[i+j*MM]=wind[dir][k][lin];
                                BDRY[i+j*MM]=Ucoarse[dir][k][lin];
                            }
                            else{
                                if (wind[dir][k][lin]!=0){
                                    UFM[i+j*MM]=1000;
                                }
                            }
                        }
                        else if (i==(p+1)*K+p && j==(m+1)*K+m){
                            
                            
                            if (wind[dir][k][lin]==8  ||  wind[dir][k][lin]==7 || wind[dir][k][lin]==3 || wind[dir][k][lin]==1 || wind[dir][k][lin]==4){
                                UFM[i+j*MM]=Ucoarse[dir][k][lin];
                                WIND[i+j*MM]=wind[dir][k][lin];
                                BDRY[i+j*MM]=Ucoarse[dir][k][lin];
                            }
                            else{
                                if (wind[dir][k][lin]!=0){
                                    UFM[i+j*MM]=1000;
                                }
                            }

                        }
                        else if (i==p*K+p){
                            
                            if (wind[dir][k][lin]==2 || wind[dir][k][lin]==3 || wind[dir][k][lin]==6){// || wind[0][(j-m)%K][(i-p)/K+(j-m)/K*(N+1)]==5 || wind[0][(j-m)%K][(i-p)/K+(j-m)/K*(N+1)]==7 ){
                                UFM[i+j*MM]=Ucoarse[dir][k][lin];
                                WIND[i+j*MM]=wind[dir][k][lin];
                                BDRY[i+j*MM]=Ucoarse[dir][k][lin];
                            }
                            else{
                                if (wind[dir][k][lin]==0){
                                UFM[i+j*MM]=1000;
                                }
                                
                            }
                        }
                        else if (i==(p+1)*K+p){
                            
                           
                            if (wind[dir][k][lin]==1 || wind[dir][k][lin]==4 || wind[dir][k][lin]==8){// || wind[0][(j-m)%K][(i-p)/K+(j-m)/K*(N+1)]==7 || wind[0][(j-m)%K][(i-p)/K+(j-m)/K*(N+1)]==5){
                                UFM[i+j*MM]=Ucoarse[dir][k][lin];
                                WIND[i+j*MM]=wind[dir][k][lin];
                                BDRY[i+j*MM]=Ucoarse[dir][k][lin];
                            }
                            else{
                                if (wind[dir][k][lin]==0){
                                    UFM[i+j*MM]=1000;
                                }
                            }
                        }
                        else if (j==m*K+m){
                           
                            if (wind[dir][k][lin]==1 || wind[dir][k][lin]==2 || wind[dir][k][lin]==5){// || wind[1][(i-p)%K][(i-p)/K+(j-m)/K*(N+1)]==8 || wind[1][(i-p)%K][(i-p)/K+(j-m)/K*(N+1)]==6){
                                UFM[i+j*MM]=Ucoarse[dir][k][lin];
                                WIND[i+j*MM]=wind[dir][k][lin];
                                BDRY[i+j*MM]=Ucoarse[dir][k][lin];
                            }
                            else{
                                if (wind[1][(i-p)%K][(i-p)/K+(j-m)/K*(N+1)]==0){
                                    UFM[p*K+p+j*MM]=1000;
                                }
                            }
                        }
                        else if (j==(m+1)*K+m){
                           
                            if (wind[dir][k][lin]==3 || wind[dir][k][lin]==4 || wind[dir][k][lin]==7){// || wind[1][(i-p)%K][(i-p)/K+(j-m)/K*(N+1)]==8 || wind[1][(i-p)%K][(i-p)/K+(j-m)/K*(N+1)]==6){
                                UFM[i+j*MM]=Ucoarse[dir][k][lin];
                                WIND[i+j*MM]=wind[dir][k][lin];
                                BDRY[i+j*MM]=Ucoarse[dir][k][lin];
                            }
                            else{
                                if (wind[dir][k][lin]==0){
                                    UFM[p*K+p+j*MM]=1000;
                                }
                                UFM[i+((m+1)*K+m)*MM]=1000;
                            }

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
#pragma omp parallel for collapse(2)
    for (int m=0; m<=N-1; m++){ //should be parallel for loop
        for (int p=0; p<=N-1; p++){ //should be parallel for loop
           
            //int full=1;
            for (int iter=0; iter<10;iter++){
                double maxError=0;
                 double currError;
               for (int s1=1; s1>=-1; s1-=2){
                 for (int s2=1; s2>=-1; s2-=2){
                 //s1=-1;
                 //s2=1;
                        for (int i =(s1<0 ? (p+1)*K+p : p*K+p); (s1<0 ? i>=p*K+p : i<= (p+1)*K+p); i+=s1){
                            for (int j=(s2<0 ? (m+1)*K+m : m*K+m); (s2<0 ? j>=m*K+m : j<= (m+1)*K+m); j+=s2)
                            {
                                if (obs[(i-p)+(j-m)*M]>1)
                                {
                                    //don't do anything; point is in obstacle
                                }
                                else if (BC[(i-p)+(j-m)*M]>=0)
                                {
                                    UFM[i+j*MM]=BC[(i-p)+(j-m)*M];
                                }
                               
                                
                                else if (BC[(i-p)+(j-m)*M]<0)
                                {
                                    double m1,m2, unew, old, discr;
                                    int tempwindx,tempwindy, windtemp;
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
                                            currError=0;
                                            UFM[i+j*MM]=old;
                                            maxError=max(currError,maxError);
                                        }
                                        else{

                                            WIND[i+j*MM]=windtemp;
                                            currError=abs(unew-old);
                                            UFM[i+j*MM]=unew;
                                            maxError=max(currError,maxError);
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



 
                                        old=UFM[i+j*MM];
                                        if (old<unew){
                                            currError=0;
                                            UFM[i+j*MM]=old;
                                            maxError=max(currError,maxError);
                                        }
                                        else{
                                            WIND[i+j*MM]=windtemp;
                                            currError=abs(unew-old);
                                            UFM[i+j*MM]=unew;
                                            maxError=max(currError,maxError);
                                        }
                                    }
                                }
                            }
                          }
                       }
                      }
                if (maxError<1e-12){
                    totalsweeps[p+m*(N)]=iter;
                    //full=0;
                    break;
                }
                else if (iter==9){
                    totalsweeps[p+m*(N)]=iter;
                    //full=0;
                }
                    }
          
                  }
                }



//    cout << "This is FINE[400,563] right after exiting function "<< FINE[399+562*1121]<< endl;
    //now take minimums along shared boundaries
 //cout << "This is FINE[0] right before entering loop "<< FINE[0]<< endl;
#pragma omp parallel for collapse(2)
    for (int m=0; m<=N-1; m++){
        for (int p=0; p<=N-1; p++){
            for (int i=p*K+p; i<=(p+1)*K+p; i++){
                for (int j=m*K+m; j<=(m+1)*K+m; j++)
                {
                    if (BC[(i-p)+(j-m)*M]>=0)
                    {
                        
                        FINE[(i-p)+(j-m)*M]=BC[(i-p)+(j-m)*M];
                    }
                    else{
                        //if (BDRY[i+j*MM]<0){
                        double old=fine[(i-p)+(j-m)*M];
                        int cont=1;
                        
                       // double old;
                        //int cont;
                        
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

   
 

  //  cout << "This is FINE[400,563] right after exiting loop "<< FINE[399+562*1121]<< endl;

    delete[] WIND;
    delete[] BDRY;
    delete[] fine;
    delete[] UFM;
}

//#endif /* fineupdate_h */
