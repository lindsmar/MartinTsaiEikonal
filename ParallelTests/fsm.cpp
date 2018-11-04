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
    prev = new double [N*N];
    for (int v=0; v<N*N; v++){
        prev[v]=u[v];
    }
double maxError=0;
    double tic,toc;
int s1, s2, i,j;
double unew;
double m1,m2;
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
            for( i=(s1<0 ? N -1 : 0 ); (s1<0 ? i >= 0 : i <= N-1); i+=s1 ){
                for( j=(s2<0 ? N -1 : 0 ); (s2<0 ? j >= 0 : j <= N-1); j+=s2 )
                {
                    if (obs[i+j*N]>1)
                    {
                        //don't do anything; point is in obstacle
                    }
                    else if (BC[i+j*N]>=0)
                    {
                        u[i+j*N]=BC[i+j*N];
                    }
                    else if (BC[i+j*N]<0)
                    {
                        if (j==0){
                            m1=u[i+(j+1)*N];
                        }
                        else if (j==N-1){
                            m1=u[i+(j-1)*N];
                        }
                        else{
                            m1=(s2>0 ? u[i+(j-1)*N] : u[i+(j+1)*N]);
                        }
                        if (i==0){
                            m2=u[(i+1)+j*N];
                        }
                        else if (i==N-1){
                            m2=u[(i-1)+j*N];
                        }
                        else{
                            m2=(s1>0 ? u[(i-1)+j*N] : u[(i+1)+j*N]);
                        }
                        
                        discr=2*pow(f[i+j*N],2)*pow(h,2)-pow(m1-m2,2);
                        
                        
                        if (abs(m1-m2)<f[i+j*N]*h)
                        {
                            unew=0.5*(m1+m2+sqrt(discr));
                            old=u[i+j*N];
                            if (old<unew){
                                u[i+j*N]=old;
                            }
                            else{
                              
                                u[i+j*N]=unew;
                            }
                            
                        }
                        else
                        {
                            unew=min(m1,m2)+f[i+j*N]*h;
                            old=u[i+j*N];
                            if (old<unew){
                                u[i+j*N]=old;
                            }
                            else{
                              
                                u[i+j*N]=unew;
                            }
                            
                            
                        }
                        
                        
                        
                       // currError=u[i+j*N]-unew;
                        //currError=fabs(currError);
                        //maxError=max(maxError,currError);
                        //u[i+j*N]=min(u[i+j*N],unew);
                    }
                }
            }
        }
       
    }
    
    toc = (omp_get_wtime() - tic);
    sumtime+= toc;
    
    double diff=0;
    for (int v=0; v<N*N; v++){
        diff=diff + fabs(u[v]-prev[v]);
        prev[v]=u[v];
    }
    diff=diff/double (N*N);
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
