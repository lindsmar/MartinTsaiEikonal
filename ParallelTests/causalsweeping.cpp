//#ifndef causalsweeping_H
//#define causalsweeping_H

//#define DIM1 2
//#define DIM2 800




//void causalsweeping(double *Ucoarse[DIM1][DIM2], double *BC,double *obs, int N,int K,double *wind[DIM1][DIM2],double* R);
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


void causalsweeping(double *Ucoarse[DIM1][DIM2], double *BC,double *obs, int N,int K,double *wind[DIM1][DIM2],double *R){
    int M;
    M=N*K+1;
    double h;
    h=(double)1/(double) (N*K);
    int s1,s2,i,j,i1;
    int DIRkLIN[3];
    int dir,k,lin;
    int DIRkLINnb[3];
    int dirnb,knb,linnb;
    double tic,toc,sumtime;
    sumtime=0;
    for (int sweep=0; sweep<1; sweep++){
        for (s2 = -1; s2 <= 1; s2 += 2){
            for (s1 = -1; s1 <= 1; s1 += 2){
                //s1=-1;
                //s2=1;
                for (i = (s1 < 0 ? M-1 : 0); (s1 < 0 ? i>=0 : i<= M-1);  i+=K*s1){
                    for (j = (s2 <0 ? M-1 : 0); (s2 < 0 ? j>=0 : j<= M-1); j+=s2){

                        if (BC[i+j*M]<0){
                           
                            cartto2CG(i,j,N+1,K,DIRkLIN);
                            dir=DIRkLIN[0];
                            k=DIRkLIN[1];
                            lin=DIRkLIN[2];
                           
                            if (wind[dir][k][lin]==1 || wind[dir][k][lin]==5 || wind[dir][k][lin]==2){
                                if ((j!=0)){
                                
                                     cartto2CG(i,j-1,N+1,K,DIRkLINnb);
                                    dirnb=DIRkLINnb[0];
                                    knb=DIRkLINnb[1];
                                    linnb=DIRkLINnb[2];
                                    
                                    
                                    if (Ucoarse[dir][k][lin]<Ucoarse[dirnb][knb][linnb]){
                                        Ucoarse[dir][k][lin]=Ucoarse[dirnb][knb][linnb];
                                    }
                                }
                            }
                            if (wind[dir][k][lin]==4 || wind[dir][k][lin]==7 || wind[dir][k][lin]==3){
                                if (j!=M-1){
                                    
                                    cartto2CG(i,j+1,N+1,K,DIRkLINnb);
                                    dirnb=DIRkLINnb[0];
                                    knb=DIRkLINnb[1];
                                    linnb=DIRkLINnb[2];
                                    
                                    if (Ucoarse[dir][k][lin]<Ucoarse[dirnb][knb][linnb]){
                                        Ucoarse[dir][k][lin]=Ucoarse[dirnb][knb][linnb];
                                    }
                                }
                            }
                            if (j%K==0){
                                if (i!=M-1){
                                    if (wind[dir][k][lin]==4 || wind[dir][k][lin]==8 || wind[dir][k][lin]==1){
                                    
                                        cartto2CG(i+1,j,N+1,K,DIRkLINnb);
                                        dirnb=DIRkLINnb[0];
                                        knb=DIRkLINnb[1];
                                        linnb=DIRkLINnb[2];
                                   
                                        if (Ucoarse[dir][k][lin]<Ucoarse[dirnb][knb][linnb]){
                                            Ucoarse[dir][k][lin]=Ucoarse[dirnb][knb][linnb];
                                        }
                                    }
                                }
                                if (i!=0){
                                    cartto2CG(i-1,j,N+1,K,DIRkLINnb);
                                  
                                    dirnb=DIRkLINnb[0];
                                    knb=DIRkLINnb[1];
                                    linnb=DIRkLINnb[2];
                                  
                                    if (wind[dir][k][lin]==3 || wind[dir][k][lin]==6 || wind[dir][k][lin]==2){
                                        if (Ucoarse[dir][k][lin]<Ucoarse[dirnb][knb][linnb]){
                                            Ucoarse[dir][k][lin]=Ucoarse[dirnb][knb][linnb];//+R[dir][k][lin]*h;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (s1<0 ? (i-K>=0) : (i+K <= M-1)){
                        for (i1 = (s1 < 0 ? i-1 : i+1); (s1 < 0 ? i1>i-K : i1< i+K);  i1+=s1){
                            for (j = (s2 <0 ? M-1 : 0); (s2 < 0 ? j>=0 : j<= M-1); j+=K*s2){

                                if (BC[i1+j*M]<0){
                                    
                                    cartto2CG(i1,j,N+1,K,DIRkLIN);
                                    dir=DIRkLIN[0];
                                    k=DIRkLIN[1];
                                    lin=DIRkLIN[2];
                                    if (wind[dir][k][lin]==4 || wind[dir][k][lin]==8 || wind[dir][k][lin]==1){// || wind[dir][k][lin]==5) //only added ==5 for maze problem
                                        if (i1!=M-1){
                                            cartto2CG(i1+1,j,N+1,K,DIRkLINnb);
                                            dirnb=DIRkLINnb[0];
                                            knb=DIRkLINnb[1];
                                            linnb=DIRkLINnb[2];
                                          
                                            if (Ucoarse[dir][k][lin]<Ucoarse[dirnb][knb][linnb]){
                                                Ucoarse[dir][k][lin]=Ucoarse[dirnb][knb][linnb];//+R[dir][k][lin]*h;
                                            }
                                        }
                                    }
                                    if (wind[dir][k][lin]==3 || wind[dir][k][lin]==6 || wind[dir][k][lin]==2 ){

                                        if (i1!=0){
                                          
                                            cartto2CG(i1-1,j,N+1,K,DIRkLINnb);
                                            dirnb=DIRkLINnb[0];
                                            knb=DIRkLINnb[1];
                                            linnb=DIRkLINnb[2];
                                          
                                            if (Ucoarse[dir][k][lin]<Ucoarse[dirnb][knb][linnb]){
                                                Ucoarse[dir][k][lin]=Ucoarse[dirnb][knb][linnb];//+R[dir][k][lin]*h;
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


   // cout << "Cartto2CG time in causal sweep: " << sumtime << endl;

}



//#endif
