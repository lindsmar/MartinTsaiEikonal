//
// Created by Lindsay Martin on 09/06/17.
//

//#ifndef EIKONALPROJECT_COARSEUPDATENEWTHETA_H
//#define EIKONALPROJECT_COARSEUPDATENEWTHETA_H
//
// Created by Lindsay Martin on 10/23/16.
//

//#include "consistent.h"
//#include "smooth.h"
//#define DIM1 2
//#define DIM2 80
/*void coarseupdatenewtheta(double *Ucoarse[DIM1][DIM2], double *Ucoarse1[DIM1][DIM2], double *Ufine1, double *Ufine2,  double *UK1[DIM1][DIM2], double *UK2[DIM1][DIM2], double *UK3[DIM1][DIM2],  double * R, double *BC,double *obs, int N,int K, int p, int z, double *wind[DIM1][DIM2], double *fwind, int sweepiter, double epsilon, double delta, double x0, double *theta[DIM1][DIM2], int wt1, int wt2, int wt3);*/
#include "sweepdefs.h"//"/scratch/04362/lmartin/paralleltests/coarsegridrevise/sweepdefs.h"

inline double smooth(double theta, double amplifier, double delta, double x0){
    
    double sigma;
    sigma=1/(1+exp((theta-x0)/amplifier));
    
    double thetaused;
    thetaused=sigma*theta+delta*(1-sigma)*theta;
    
    return thetaused;
    
}







inline void cartto2CG(int i, int j, int k, int N, int K,int *DHBlin){
    
    DHBlin[0]=i%K;
    DHBlin[1]=j%K;
    DHBlin[2]=k%K;
    DHBlin[3]=i/K+j/K*N+k/K*N*N;
    
    return;
}






void coarseupdatenewtheta(double *Ucoarse[DIM1][DIM2][DIM3], double *Ucoarse1[DIM1][DIM2][DIM3], double *Ufine1, double *Ufine2,  double *UK1[DIM1][DIM2][DIM3], double *UK2[DIM1][DIM2][DIM3], double *UK3[DIM1][DIM2][DIM3],  double * R, double *BC,double *obs, int N,int K, int p, int z, int r, double *wind[DIM1][DIM2][DIM3], double *fwind, int sweepiter, double epsilon, double delta, double x0,  double *theta[DIM1][DIM2][DIM3], int wt1, int wt2, int wt3,int *totalsweeps)
{
    int M;
    double H,h,old,m1,m2,m3,M1,M2,M3,discr,hx,hy,unew,F, utemp;
    int tempwindx, tempwindy, istart1,istart2,jstart1,jstart2,kstart1,kstart2;
    double prevW, newwind,prevw;
    double m11,m12,m21,m22,m31,m32;

    
    bool consis;
    double currError;

    H=(double)1/(double) N;
    h=(double)1/(double)(N*K);
    M=N*K+1;
    double T, wtdsum, lower;
    double *marked, *Uktemp1, *tempwind;
    Uktemp1= new double [(N+1)*(N+1)*(N+1)];
    marked= new double [(N+1)*(N+1)*(N+1)];
    
    tempwind= new double [(N+1)*(N+1)*(N+1)];
    /*for (int c=0; c<(N+1)*(N+1); c++){
        marked[c]=0;
        tempwind[c]=0;
    }*/
    
    int s1,s2,s3,i,j,k;


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
    if (r==0){
        kstart1=M-1;
        kstart2=0;
    }
    else{
        kstart1=M-1-(K-r);
        kstart2=r;
    }
    int dir,shift,ishift,jshift;
    
    int DHBlin[4];
    int Hs,Bs,Ds,lin;
    int DHBlinnb[4];
    int Hsnb,Bsnb,Dsnb,linnb;

    //for (int iter=0; iter<sweepiter; iter++) {
         double maxError=0;
      // for (s2 = -1; s2 <= 1; s2 += 2){
        //for (s1 = -1; s1 <= 1; s1 += 2){
          //for (s2 = 1; s2 >= -1; s2 -= 2)
          // for (s1 = 1; s1 >= -1; s1 -= 2)
        //s2=1;
        //s1=1;
    s1=1;
    s2=1;
    s3=1;
                for (i = (s1 < 0 ? istart1 : istart2); (s1 < 0 ? i >= 0 : i <= M - 1); i += K * s1){
                    for (j = (s2 < 0 ? jstart1 : jstart2); (s2 < 0 ? j >= 0 : j <= M - 1); j += K * s2) {
                        for (k = (s3 < 0 ? kstart1 : kstart2); (s3 < 0 ? k >= 0 : k <= M - 1); k += K * s3) {
                            if (obs[i + j * M + k*M*M] > 1) {
                                //don't do anything; point is in obstacle
                            }  else if (BC[i + j * M + k*M*M] >= 0) {
                                //tic=omp_get_wtime();
                                cartto2CG(i,j,k,N+1,K,DHBlin);
                                //toc=omp_get_wtime()-tic;
                                //toc;
                                Ds=DHBlin[0];
                                Hs=DHBlin[1];
                                Bs=DHBlin[2];
                                lin=DHBlin[3];
                                Ucoarse[Ds][Hs][Bs][lin] = BC[i + j * M + k*M*M];
                                
                                
                            } else if (BC[i + j * M + k*M*M] < 0) {
                                // tic=omp_get_wtime();
                                cartto2CG(i,j,k,N+1,K,DHBlin);
                                //toc=omp_get_wtime()-tic;
                                //carttime+=toc;
                                Ds=DHBlin[0];
                                Hs=DHBlin[1];
                                Bs=DHBlin[2];
                                lin=DHBlin[3];
                                
                                
                                
                                tempwindx = 0;
                                tempwindy = 0;
                                // tempwindz = 0;
                                if (j == z) {
                                    //tic=omp_get_wtime();
                                    cartto2CG(i,j+K,k,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m1 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    tempwindy = 7;
                                    
                                } else if (j == M - 1 - (K - z) && z != 0) {
                                    // tic=omp_get_wtime();
                                    cartto2CG(i,j-K,k,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    // carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m1 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    
                                    tempwindy = 5;
                                } else if (j == M - 1) {
                                    // tic=omp_get_wtime();
                                    cartto2CG(i,j-K,k,N+1,K,DHBlinnb);
                                    // toc=omp_get_wtime()-tic;
                                    // carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m1 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    
                                    
                                    
                                    tempwindy = 5;
                                } else {
                                    // tic=omp_get_wtime();
                                    cartto2CG(i,j+K,k,N+1,K,DHBlinnb);
                                    // toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m11 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    // tic=omp_get_wtime();
                                    cartto2CG(i,j-K,k,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    // carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m12 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    m1 = min(m11,m12);
                                    
                                    
                                }
                                
                                if (i == p) {
                                    //tic=omp_get_wtime();
                                    cartto2CG(i+K,j,k,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m2 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    
                                    
                                    tempwindx = 8;
                                } else if (i == M - 1 - (K - p) && p != 0) {
                                    //tic=omp_get_wtime();
                                    cartto2CG(i-K,j,k,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m2 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    
                                    tempwindx = 6;
                                    
                                } else if (i == M - 1) {
                                    //tic=omp_get_wtime();
                                    cartto2CG(i-K,j,k,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m2 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    tempwindx = 6;
                                    
                                } else {
                                    //tic=omp_get_wtime();
                                    cartto2CG(i-K,j,k,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m21 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    //tic=omp_get_wtime();
                                    cartto2CG(i+K,j,k,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m22 = Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    m2 =min(m21,m22);
                                    
                                    
                                }
                                
                                if (k == r) {
                                    //tic=omp_get_wtime();
                                    cartto2CG(i,j,k+K,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m3= Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    
                                    tempwindx = 8;
                                } else if (k == M - 1 - (K - r) && r != 0) {
                                    //tic=omp_get_wtime();
                                    cartto2CG(i,j,k-K,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m3= Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    
                                    tempwindx = 6;
                                    
                                } else if (k == M - 1) {
                                    //tic=omp_get_wtime();
                                    cartto2CG(i,j,k-K,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m3= Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    tempwindx = 6;
                                    
                                } else {
                                    //tic=omp_get_wtime();
                                    cartto2CG(i,j,k-K,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m31= Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    //tic=omp_get_wtime();
                                    cartto2CG(i,j,k+K,N+1,K,DHBlinnb);
                                    //toc=omp_get_wtime()-tic;
                                    //carttime+=toc;
                                    Dsnb=DHBlinnb[0];
                                    Hsnb=DHBlinnb[1];
                                    Bsnb=DHBlinnb[2];
                                    linnb=DHBlinnb[3];
                                    m32= Ucoarse[Dsnb][Hsnb][Bsnb][linnb];
                                    m3 =min(m31,m32);
                                    
                                    
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
                                
                                unew=M1+R[i + j * M + k*M*M]*H;
                            
                            if (unew <= M2){
                                old=Ucoarse[Ds][Hs][Bs][lin];
                                if (fabs(unew-UK1[Ds][Hs][Bs][lin])<1e-8){// || abs(Ufine[i+j*M]-UK[i+j*M])<1e-13){
                                    utemp =Ufine1[i + j * M+k*M*M];
                                }
                                else{
                                    wtdsum=(wt1*(unew-UK1[Ds][Hs][Bs][lin])+wt2*(UK1[Ds][Hs][Bs][lin]-UK2[Ds][Hs][Bs][lin])+wt3*(UK2[Ds][Hs][Bs][lin]-UK3[Ds][Hs][Bs][lin]))/(wt1+wt2+wt3);
                                    T=(Ufine1[i+j*M+k*M*M]-Ufine2[i+j*M+k*M*M])/wtdsum;
                                    //T=(ORIG[i+j*M]-Ufine[i+j*M])/(unew-UK[i+j*M]);
                                    //T=T-.01;
                                    lower=(Ucoarse1[Ds][Hs][Bs][lin]-Ufine1[i+j*M+k*M*M])/(unew-UK1[Ds][Hs][Bs][lin]);
                                    if (T<0){
                                        T=0;
                                    }
                                    else{
                                        T=smooth(T,epsilon,delta,x0);
                                    }
                                    if (epsilon==0){
                                        T=0;
                                    }
                                    theta[Ds][Hs][Bs][lin]=T;
                                    utemp = T * unew + Ufine1[i+j*M+k*M*M] - T * UK1[Ds][Hs][Bs][lin];
                                    
                                    
                                   
                                }
                                
                                Uktemp1[lin] = unew;
                                currError=abs(unew-old);
                                Ucoarse[Ds][Hs][Bs][lin]= utemp;
                                maxError=max(currError,maxError);
                            }
                                
                                
                                else{
                                    discr=2*pow(R[i+j*N+k*N*N],2)*pow(H,2)-pow(M1-M2,2);
                                    unew=0.5*(M1+M2+sqrt(discr));
                                    if (unew <= M3){
                                        old=Ucoarse[Ds][Hs][Bs][lin];
                                        if (fabs(unew-UK1[Ds][Hs][Bs][lin])<1e-8){// || abs(Ufine[i+j*M]-UK[i+j*M])<1e-13){
                                            utemp =Ufine1[i + j * M+k*M*M];
                                        }
                                        else{
                                            wtdsum=(wt1*(unew-UK1[Ds][Hs][Bs][lin])+wt2*(UK1[Ds][Hs][Bs][lin]-UK2[Ds][Hs][Bs][lin])+wt3*(UK2[Ds][Hs][Bs][lin]-UK3[Ds][Hs][Bs][lin]))/(wt1+wt2+wt3);
                                            T=(Ufine1[i+j*M+k*M*M]-Ufine2[i+j*M+k*M*M])/wtdsum;
                                            //T=(ORIG[i+j*M]-Ufine[i+j*M])/(unew-UK[i+j*M]);
                                            //T=T-.01;
                                            lower=(Ucoarse1[Ds][Hs][Bs][lin]-Ufine1[i+j*M+k*M*M])/(unew-UK1[Ds][Hs][Bs][lin]);
                                            if (T<0){
                                                T=0;
                                            }
                                            else{
                                                T=smooth(T,epsilon,delta,x0);
                                            }
                                            if (epsilon==0){
                                                T=0;
                                            }
                                            theta[Ds][Hs][Bs][lin]=T;
                                            utemp = T * unew + Ufine1[i+j*M+k*M*M] - T * UK1[Ds][Hs][Bs][lin];
                                            
                                            
                                           
                                        }
                                        
                                        Uktemp1[lin] = unew;
                                        currError=abs(unew-old);
                                        Ucoarse[Ds][Hs][Bs][lin]= utemp;
                                        maxError=max(currError,maxError);
                                        
                                        
                                    }
                                    else{
                                        discr=4*pow(M1+M2+M3,2)-12*(M1*M1+M2*M2+M3*M3-pow(R[i+j*N+k*N*N],2)*pow(H,2));
                                        unew=(2*(M1+M2+M3)+sqrt(discr))/6;
                                        if (unew<1000){
                                            old=Ucoarse[Ds][Hs][Bs][lin];
                                            if (fabs(unew-UK1[Ds][Hs][Bs][lin])<1e-8){// || abs(Ufine[i+j*M]-UK[i+j*M])<1e-13){
                                                utemp =Ufine1[i + j * M+k*M*M];
                                            }
                                            else{
                                                wtdsum=(wt1*(unew-UK1[Ds][Hs][Bs][lin])+wt2*(UK1[Ds][Hs][Bs][lin]-UK2[Ds][Hs][Bs][lin])+wt3*(UK2[Ds][Hs][Bs][lin]-UK3[Ds][Hs][Bs][lin]))/(wt1+wt2+wt3);
                                                T=(Ufine1[i+j*M+k*M*M]-Ufine2[i+j*M+k*M*M])/wtdsum;
                                                //T=(ORIG[i+j*M]-Ufine[i+j*M])/(unew-UK[i+j*M]);
                                                //T=T-.01;
                                                lower=(Ucoarse1[Ds][Hs][Bs][lin]-Ufine1[i+j*M+k*M*M])/(unew-UK1[Ds][Hs][Bs][lin]);
                                                if (T<0){
                                                    T=0;
                                                }
                                                else{
                                                    T=smooth(T,epsilon,delta,x0);
                                                }
                                                if (epsilon==0){
                                                    T=0;
                                                }
                                                theta[Ds][Hs][Bs][lin]=T;
                                                utemp = T * unew + Ufine1[i+j*M+k*M*M] - T * UK1[Ds][Hs][Bs][lin];
                                                
                                                
                                               
                                            }
                                            
                                            Uktemp1[lin] = unew;
                                            currError=abs(unew-old);
                                            Ucoarse[Ds][Hs][Bs][lin]= utemp;
                                            maxError=max(currError,maxError);
                                            
                                        }
                                    }
                                }
                                
                                
                            }
                            
                            
                      }
                    }
                    }
    



  


    for (int k=0; k<(N+1)*(N+1)*(N+1); k++){
        UK3[p][z][r][k]= UK2[p][z][r][k];
        UK2[p][z][r][k]= UK1[p][z][r][k];
        UK1[p][z][r][k]= Uktemp1[k];
        wind[p][z][r][k]=tempwind[k];
    }

    totalsweeps[0]=sweepiter;
    
    delete[] marked;
    delete[] Uktemp1;
    delete[] tempwind;
    return;
}


//#endif //EIKONALPROJECT_COARSEUPDATENEWTHETA_H
