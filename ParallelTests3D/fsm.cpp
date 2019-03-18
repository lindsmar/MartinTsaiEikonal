//
//  fineupdate.h
//  Eikonal Project 5
//
//  Created by Lindsay Martin on 6/22/16.
//  Copyright © 2016 Lindsay Martin. All rights reserved.
//

/*#ifndef fsm_h
#define fsm_h
#include <algorithm>
#include <cmath>

using namespace std;

int fsm(double *u, double *f, double *BC, double *obs, int iter, int N, double h);*/
#include "sweepdefs.h"


int fsm(double *u, double *f, double *BC, double *obs, int iter, int N, double h){
    double sumtime=0;
    
    double *prev;
    prev = new double [N*N*N];
    for (int v=0; v<N*N*N; v++){
        prev[v]=u[v];
    }
double maxError=0;
    double tic,toc;
int s1, s2,s3, i,j,k;
double unew;
double m1,m2,m3,M1,M2,M3;
double discr;
double old;
    int totaliters=0;
for (int s=0;s<iter;s++){
     tic = omp_get_wtime();
    totaliters++;
    
    //maxError=0;
    double currError;
    
    
    for (s1=1; s1>=-1; s1-=2){
        for (s2=1; s2>=-1; s2-=2){
            for (s3=1; s3>=-1; s3-=2){
                for( i=(s1<0 ? N -1 : 0 ); (s1<0 ? i >= 0 : i <= N-1); i+=s1 ){
                    for( j=(s2<0 ? N -1 : 0 ); (s2<0 ? j >= 0 : j <= N-1); j+=s2 ){
                        for( k=(s3<0 ? N -1 : 0 ); (s3<0 ? k >= 0 : k <= N-1); k+=s3 ){
                        if (obs[i+j*N+k*N*N]>1)
                        {
                            //don't do anything; point is in obstacle
                        }
                        else if (BC[i+j*N+k*N*N]>=0)
                        {
                            u[i+j*N+k*N*N]=BC[i+j*N+k*N*N];
                            
                        }
                        else if (BC[i+j*N+k*N*N]<0)
                        {
                            if (j==0){
                                m1=u[i+(j+1)*N+k*N*N];
                            }
                            else if (j==N-1){
                                m1=u[i+(j-1)*N+k*N*N];
                            }
                            else{
                                m1=min(u[i+(j-1)*N+k*N*N],u[i+(j+1)*N+k*N*N]);
                                
                            }
                            
                            if (i==0){
                                m2=u[(i+1)+j*N+k*N*N];
                            }
                            else if (i==N-1){
                                m2=u[(i-1)+j*N+k*N*N];
                            }
                            else{
                                m2=min(u[(i-1)+j*N+k*N*N],u[(i+1)+j*N+k*N*N]);
                                
                            }
                            
                            if (k==0){
                                m3=u[i+j*N+(k+1)*N*N];
                            }
                            else if (k==N-1){
                                m3=u[i+j*N+(k-1)*N*N];
                            }
                            else{
                                m3=min(u[i+j*N+(k-1)*N*N],u[i+j*N+(k+1)*N*N]);
                                
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
                            
                            unew=M1+f[i+j*N+k*N*N]*h;
                            
                        
                            if (unew <= M2){
                                old=u[i+j*N+k*N*N];
                                if (old<unew){
                                    u[i+j*N+k*N*N]=old;
                                }
                                else{
                                    
                                    u[i+j*N+k*N*N]=unew;
                                    
                                }
                            }
                            else{
                                discr=2*pow(f[i+j*N+k*N*N],2)*pow(h,2)-pow(M1-M2,2);
                                unew=0.5*(M1+M2+sqrt(discr));
                                if (unew <= M3){
                                    old=u[i+j*N+k*N*N];
                                    if (old<unew){
                                        u[i+j*N+k*N*N]=old;
                                    }
                                    else{
                                        
                                        u[i+j*N+k*N*N]=unew;
                                       
                                    }
                                }
                                else{
                                    discr=4*pow(M1+M2+M3,2)-12*(M1*M1+M2*M2+M3*M3-pow(f[i+j*N+k*N*N],2)*pow(h,2));
                                    unew=(2*(M1+M2+M3)+sqrt(discr))/6;
                                        if (unew<1000){
                                        old=u[i+j*N+k*N*N];
                                        if (old<unew){
                                            u[i+j*N+k*N*N]=old;
                                        }
                                        else{
                                            
                                            u[i+j*N+k*N*N]=unew;
                                            
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
    
    toc = (omp_get_wtime() - tic);
    sumtime+= toc;
    
    double diff=0;
    
    for (int v=0; v<N*N*N; v++){
        diff=diff + fabs(u[v]-prev[v]);
        prev[v]=u[v];
    }
    
    diff=diff/double (N*N*N);
     //cout << diff << endl;
  
    //cout << maxError << endl;
    if (diff<1e-14){
        //cout<< "fineupdate time: " <<sumtime << endl;
        break;
    }
    
}

    cout<< "Fastsweeping time: " <<sumtime << endl;
return totaliters;
}
//#endif /* fsm_h */
