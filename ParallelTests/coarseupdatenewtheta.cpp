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
inline bool consistent(double prevW, double prevw, int p, int z){
    
    bool ans;
    
    if (p==0 && z==0){
        
        if (prevW==1){
            if (prevw==1||prevw==5||prevw==8){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==2){
            if (prevw==2||prevw==5||prevw==6){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==3){
            if (prevw==3||prevw==6||prevw==7){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==4){
            if (prevw==4||prevw==7||prevw==8){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==8){
            if (prevw==1||prevw==4||prevw==8){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==5){
            if (prevw==1||prevw==2||prevw==5){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==6){
            if (prevw==2||prevw==3||prevw==6){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==7){
            if (prevw==3||prevw==4||prevw==7){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        //ans=(prevW==prevw);
    }
    else if (p==0){
        if (prevW==5){
            if (prevw==5||prevw==1 || prevw==2){//(prevw==5){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==7){
            if (prevw==7 ||prevw==3 || prevw==4){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==1){
            if (prevw==1||prevw==4||prevw==8|| prevw==5){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==2){
            if (prevw==2||prevw==3||prevw==6 || prevw==5){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==3){
            if (prevw==2||prevw==3||prevw==6||prevw==7){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==4){
            if (prevw==1||prevw==4||prevw==8||prevw==7){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==6){
            if (prevw==2||prevw==3||prevw==6){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==8){
            if (prevw==1||prevw==4||prevw==8){
                ans=true;
            }
            else{
                ans=false;
            }
        }
    }
    else if (z==0){
        if (prevW==6){
            if (prevw==6 ||prevw==3 || prevw==2){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==8){
            if (prevw==8 ||prevw==1 || prevw==4){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==1){
            if (prevw==1||prevw==2||prevw==5||prevw==8){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==2){
            if (prevw==1||prevw==2||prevw==5||prevw==6){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==3){
            if (prevw==3||prevw==4||prevw==7||prevw==6){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==4){
            if (prevw==3||prevw==4||prevw==7||prevw==8){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==7){
            if (prevw==3||prevw==4||prevw==7){
                ans=true;
            }
            else{
                ans=false;
            }
        }
        else if (prevW==5){
            if (prevw==1||prevw==2||prevw==5){
                ans=true;
            }
            else{
                ans=false;
            }
        }
    }
    return ans;
}



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




void coarseupdatenewtheta(double *Ucoarse[DIM1][DIM2], double *Ucoarse1[DIM1][DIM2], double *Ufine1, double *Ufine2,  double *UK1[DIM1][DIM2], double *UK2[DIM1][DIM2], double *UK3[DIM1][DIM2],  double * R, double *BC,double *obs, int N,int K, int p, int z, double *wind[DIM1][DIM2], double *fwind, int sweepiter, double epsilon, double delta, double x0,  double *theta[DIM1][DIM2], int wt1, int wt2, int wt3,int *totalsweeps)
{
    int M;
    double H,h,old,m1,m2,discr,hx,hy,unew,F, utemp, uknew;
    int tempwindx, tempwindy, istart1,istart2,jstart1,jstart2;
    double prevW, newwind,prevw;
    
    bool consis;
    double currError;

    H=(double)1/(double) N;
    h=(double)1/(double)(N*K);
    M=N*K+1;
    double T, wtdsum, lower;
    double *marked, *Uktemp1, *tempwind;

    marked= new double [(N+1)*(N+1)];
    Uktemp1= new double [(N+1)*(N+1)];
    tempwind= new double [(N+1)*(N+1)];
    for (int c=0; c<(N+1)*(N+1); c++){
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
    int dir,shift,ishift,jshift;
    if (p==0){
        dir=0;
        shift=z;
        ishift=0;
        jshift=z;
        
    }
    else{
        dir=1;
        shift=p;
        ishift=p;
        jshift=0;
    }
    
    int DIRkLIN[3];
    int k,lin;
    int DIRkLINnbx1[3],DIRkLINnbx2[3];
    int dirnbx1,knbx1,linnbx1,dirnbx2,knbx2,linnbx2;
    int DIRkLINnby1[3],DIRkLINnby2[3];
    int dirnby1,knby1,linnby1,dirnby2,knby2,linnby2;

    for (int iter=0; iter<sweepiter; iter++) {
         double maxError=0;
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
                            cartto2CG(i,j,N+1,K,DIRkLIN);
                            dir=DIRkLIN[0];
                            k=DIRkLIN[1];
                            lin=DIRkLIN[2];
                            Ucoarse[dir][k][lin] = BC[i + j * M];
                        }else if (BC[i + j * M] < 0){
                            cartto2CG(i,j,N+1,K,DIRkLIN);
                            dir=DIRkLIN[0];
                            k=DIRkLIN[1];
                            lin=DIRkLIN[2];
                            
                            tempwindx = 0;
                            tempwindy = 0;
                            if (j == z) {
                                cartto2CG(i,j+K,N+1,K,DIRkLINnby2);
                                dirnby2=DIRkLINnby2[0];
                                knby2=DIRkLINnby2[1];
                                linnby2=DIRkLINnby2[2];
                                m1 = Ucoarse[dirnby2][knby2][linnby2];
                                tempwindy = 7;
                                hx = H;
                            } else if (j == M - 1 - (K - z) && z != 0) {
                                cartto2CG(i,j-K,N+1,K,DIRkLINnby1);
                                dirnby1=DIRkLINnby1[0];
                                knby1=DIRkLINnby1[1];
                                linnby1=DIRkLINnby1[2];
                                
                                m1 =Ucoarse[dirnby1][knby1][linnby1];
                                hx = H;
                                tempwindy = 5;
                            } else if (j == M - 1) {
                                cartto2CG(i,j-K,N+1,K,DIRkLINnby1);
                                dirnby1=DIRkLINnby1[0];
                                knby1=DIRkLINnby1[1];
                                linnby1=DIRkLINnby1[2];
                                
                                m1 =Ucoarse[dirnby1][knby1][linnby1];
                                
                                
                                hx = H;
                                tempwindy = 5;
                            } else {
                                cartto2CG(i,j+K,N+1,K,DIRkLINnby2);
                                dirnby2=DIRkLINnby2[0];
                                knby2=DIRkLINnby2[1];
                                linnby2=DIRkLINnby2[2];
                                
                                cartto2CG(i,j-K,N+1,K,DIRkLINnby1);
                                dirnby1=DIRkLINnby1[0];
                                knby1=DIRkLINnby1[1];
                                linnby1=DIRkLINnby1[2];
                                m1 = min(Ucoarse[dirnby1][knby1][linnby1], Ucoarse[dirnby2][knby2][linnby2]);
                                if (Ucoarse[dirnby1][knby1][linnby1] <Ucoarse[dirnby2][knby2][linnby2]){
                                    tempwindy = 5;
                                } else {
                                    tempwindy = 7;
                                }
                                hx = H;
                            }
                            
                            if (i == p) {
                                
                                cartto2CG(i+K,j,N+1,K,DIRkLINnbx2);
                                dirnbx2=DIRkLINnbx2[0];
                                knbx2=DIRkLINnbx2[1];
                                linnbx2=DIRkLINnbx2[2];
                                
                                m2 =Ucoarse[dirnbx2][knbx2][linnbx2];
                                hy = H;
                                tempwindx = 8;
                            } else if (i == M - 1 - (K - p) && p != 0) {
                                
                                cartto2CG(i-K,j,N+1,K,DIRkLINnbx1);
                                dirnbx1=DIRkLINnbx1[0];
                                knbx1=DIRkLINnbx1[1];
                                linnbx1=DIRkLINnbx1[2];
                                
                                m2 =Ucoarse[dirnbx1][knbx1][linnbx1];
                                
                                tempwindx = 6;
                                hy = H;
                            } else if (i == M - 1) {
                                cartto2CG(i-K,j,N+1,K,DIRkLINnbx1);
                                dirnbx1=DIRkLINnbx1[0];
                                knbx1=DIRkLINnbx1[1];
                                linnbx1=DIRkLINnbx1[2];
                                
                                m2 =Ucoarse[dirnbx1][knbx1][linnbx1];
                                tempwindx = 6;
                                hy = H;
                            } else {
                                cartto2CG(i-K,j,N+1,K,DIRkLINnbx1);
                                dirnbx1=DIRkLINnbx1[0];
                                knbx1=DIRkLINnbx1[1];
                                linnbx1=DIRkLINnbx1[2];
                                
                                cartto2CG(i+K,j,N+1,K,DIRkLINnbx2);
                                dirnbx2=DIRkLINnbx2[0];
                                knbx2=DIRkLINnbx2[1];
                                linnbx2=DIRkLINnbx2[2];
                                m2 =min(Ucoarse[dirnbx1][knbx1][linnbx1],Ucoarse[dirnbx2][knbx2][linnbx2]);
                                hy = H;
                                if (Ucoarse[dirnbx1][knbx1][linnbx1] < Ucoarse[dirnbx2][knbx2][linnbx2]){
                                    tempwindx = 6;
                                } else {
                                    tempwindx = 8;
                                }
                            }
                            
                            


                            if (m1 < 500 || m2 < 500) {
                                F = R[i + j * M];
                                prevW=wind[dir][k][lin];
                                prevw=fwind[i + j * M];
                                consis=consistent(prevW,prevw,p,z);
                                double old=Ucoarse[dir][k][lin];

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
                                    if (abs(unew-UK[i+j*M])<1e-(N+1)){
                                      utemp = Ufine[i + j * M];
                                      Ucoarse[i + j * M] = utemp;
                                       UKtemp[i + j * M] = unew;
                                       tempwind[i + j * M] = newwind;
                                      continue;
                                    }*/



                                      if (marked[lin]==0){
                                      if (consis){
                                        if (prevW==newwind){
                                          unew = 0.5 * (m1 + m2 + sqrt(discr));
                                          if (fabs(unew-UK1[dir][k][lin])<1e-15){
                                            utemp = Ufine1[i + j * M];
                                          }
                                          else{
                                              wtdsum=(wt1*(unew-UK1[dir][k][lin])+wt2*(UK1[dir][k][lin]-UK2[dir][k][lin])+wt3*(UK2[dir][k][lin]-UK3[dir][k][lin]))/(wt1+wt2+wt3);
                                            T=(Ufine1[i+j*M]-Ufine2[i+j*M])/wtdsum;
                                            //T=(ORIG[i+j*M]-Ufine[i+j*M])/(unew-UK[i+j*M]);
                                            //T=T-.01;
                                            lower=(Ucoarse1[dir][k][lin]-Ufine1[i+j*M])/(unew-UK1[dir][k][lin]);
                                            if (T<0){
                                              T=0;
                                            }
                                            else{
                                              T=smooth(T,epsilon,delta,x0);
                                            }
                                              if (epsilon==0){
                                                  T=0;
                                              }
                                            theta[dir][k][lin]=T;
                                          utemp = T * unew + Ufine1[i + j * M] - T * UK1[dir][k][lin];
                                          }
                                          uknew=unew;
                                          //marked[i+j*M]=1;
                                            currError=abs(utemp-old);
                                            maxError=max(currError,maxError);
                                          Ucoarse[dir][k][lin] = utemp;
                                          Uktemp1[lin] = unew;
                                          tempwind[lin]= fwind[i+j*M];
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
                                          currError=abs(utemp-old);
                                          maxError=max(currError,maxError);
                                          Ucoarse[dir][k][lin] = utemp;
                                          Uktemp1[lin] = unew;

                                          tempwind[lin]= newwind;
                                        //}


                                    }
                                    else{
                                      //do nothing
                                    }


                                } else {
                                     double old=Ucoarse[dir][k][lin];
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

                                    if (marked[lin]==0){
                                      if (consis){
                                        if (prevW==newwind){
                                          unew = min(m1, m2) + F * H;
                                           if (fabs(unew-UK1[dir][k][lin])<1e-15){
                                            utemp = Ufine1[i + j * M];
                                          }
                                          else{
                                            wtdsum=(wt1*(unew-UK1[dir][k][lin])+wt2*(UK1[dir][k][lin]-UK2[dir][k][lin])+wt3*(UK2[dir][k][lin]-UK3[dir][k][lin]))/(wt1+wt2+wt3);
                                            T=(Ufine1[i+j*M]-Ufine2[i+j*M])/wtdsum;
                                            //T=(ORIG[i+j*M]-Ufine[i+j*M])/(unew-UK[i+j*M]);
                                            //T=T-.01;
                                            lower=(Ucoarse1[dir][k][lin]-Ufine1[i+j*M])/(unew-UK1[dir][k][lin]);
                                            if (T<0){
                                              T=0;
                                            }
                                            else{
                                              T=smooth(T,epsilon,delta,x0);
                                            }
                                              if (epsilon==0){
                                                  T=0;
                                              }
                                            theta[dir][k][lin]=T;
                                          utemp = T * unew + Ufine1[i + j * M] - T * UK1[dir][k][lin];
                                          }
                                          uknew=unew;
                                          //marked[i+j*M]=1;
                                            currError=abs(utemp-old);
                                            maxError=max(currError,maxError);
                                            Ucoarse[dir][k][lin] = utemp;
                                            Uktemp1[lin] = unew;
                                            tempwind[lin]= fwind[i+j*M];
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
                                        currError=abs(utemp-old);
                                        maxError=max(currError,maxError);
                                        Ucoarse[dir][k][lin] = utemp;
                                        Uktemp1[lin] = unew;
                                        
                                        tempwind[lin]= newwind;
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
        if (maxError<1e-12){
            totalsweeps[0]=iter;
            for (int k=0; k<(N+1)*(N+1); k++){
                
                UK3[dir][shift][k]= UK2[dir][shift][k];
                UK2[dir][shift][k]= UK1[dir][shift][k];
                UK1[dir][shift][k]= Uktemp1[k];
                wind[dir][shift][k]=tempwind[k];
            }
            
            
            delete[] marked;
            delete[] Uktemp1;
            delete[] tempwind;
            return;
        }

    }

  


    for (int k=0; k<(N+1)*(N+1); k++){
        UK3[dir][shift][k]= UK2[dir][shift][k];
        UK2[dir][shift][k]= UK1[dir][shift][k];
        UK1[dir][shift][k]= Uktemp1[k];
        wind[dir][shift][k]=tempwind[k];
    }
    totalsweeps[0]=sweepiter;
    
    delete[] marked;
    delete[] Uktemp1;
    delete[] tempwind;
    return;
}


//#endif //EIKONALPROJECT_COARSEUPDATENEWTHETA_H
