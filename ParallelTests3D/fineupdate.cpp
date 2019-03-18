//
//  fineupdate.h
//  Eikonal Project 5
//
//  Created by Lindsay Martin on 6/22/16.
//  Copyright © 2016 Lindsay Martin. All rights reserved.
//


#include "sweepdefs.h"///scratch/04362/lmartin/paralleltests/coarsegridrevise/sweepdefs.h"



inline void cartto2CG(int i, int j, int k, int N, int K,int *DHBlin){
    
    DHBlin[0]=i%K;
    DHBlin[1]=j%K;
    DHBlin[2]=k%K;
    DHBlin[3]=i/K+j/K*N+k/K*N*N;
    
    return;
}



void fineupdate(double * FINE, double *Ucoarse[DIM1][DIM2][DIM3], double *f,double *BC,double *obs, int N,int K, double *wind[DIM1][DIM2][DIM3], double *fwind, int *totalsweeps)
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
    ofstream valuef;


    //double *FINE, *FINE2, *exwind, *cellbdry;
    double *fine, *WIND, *BDRY, *UFM;

    //fine2=(double *) mxCalloc(M*M, sizeof(double));
    fine=new double [M*M*M];
    //FINE=mxCreateDoubleMatrix(M,M,mxREAL);
    //fine=mxGetPr(FINE);
    //FINE2=mxCreateDoubleMatrix(M,M,mxREAL);
    //fine2=mxGetPr(FINE2);
    UFM=new double [MM*MM*MM];
    WIND=new double [MM*MM*MM];
    BDRY=new double [MM*MM*MM];
    //int c;
    #pragma omp parallel for
    for (int c=0; c<(M+N-1)*(M+N-1)*(M+N-1); c++){
        UFM[c]=1000;
    }
    
#pragma omp parallel for
    for (int z=0; z<M*M*M; z++){
        fine[z]=1000;
        //FINE[z]=Ucoarse[z];
        
    }
    #pragma omp parallel for
    for (int z=0; z<(M+N-1)*(M+N-1)*(M+N-1); z++){
        WIND[z]=1000;
         BDRY[z]=-1;
    }

    
    
    
  
   
    
#pragma omp parallel for collapse(3)
    for (int m=0; m<=N-1; m++){
        for (int p=0; p<=N-1; p++){
             for (int v=0; v<=N-1; v++){
                 for (int i=p*K+p; i<=(p+1)*K+p; i++){
                     for (int j=m*K+m; j<=(m+1)*K+m; j++){
                         for (int k=v*K+v; k<=(v+1)*K+v; k++)
                         {
                    
                    
                            if (obs[(i-p)+(j-m)*M+(k-v)*M*M]>1)
                            {
                                //don't do anything; point is in obstacle
                            }
                            else if (BC[(i-p)+(j-m)*M+(k-v)*M*M]>=0)
                            {
                                UFM[i+j*MM+k*MM*MM]=BC[(i-p)+(j-m)*M+(k-v)*M*M];
                                BDRY[i+j*MM+k*MM*MM]=BC[(i-p)+(j-m)*M+(k-v)*M*M];
                            }
                           
                            
                            else{
                                int DHBlin[4];
                                int Hs,Bs,Ds,lin;
                                cartto2CG(i-p,j-m,k-v,N+1,K,DHBlin);
                                Ds=DHBlin[0];
                                Hs=DHBlin[1];
                                Bs=DHBlin[2];
                                lin=DHBlin[3];
                               
                                
                                if (i==p*K+p && j==m*K+m && k==v*M+v){
                                    //continue;
                                }
                                if (i==0 || j==0 || k==0){
                                    if (i==0){
                                        if (j==m*K+m ){
                                            UFM[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                            BDRY[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                        }
                                        if (k==v*K+v){
                                            UFM[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                            BDRY[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                        }
                                    }
                                    if (j==0){
                                        if (i==p*K+p){
                                            UFM[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                            BDRY[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                        }
                                        if (k==v*K+v){
                                            UFM[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                            BDRY[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                        }
                                    }
                                    if (k==0){
                                        if (i==p*K+p){
                                            UFM[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                            BDRY[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                        }
                                        if (j==m*K+m){
                                            UFM[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                            BDRY[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                        }
                                    }
                                }
                                
                                
                                else{
                                    if ((j-m)%K==0){
                                    if (j==m*K+m){
                                        UFM[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                        BDRY[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                    }
                                }
                                 if ((i-p)%K==0){
                                    if (i==p*K+p){
                                        UFM[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                        BDRY[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                    }
                                 }
                                   if ((k-v)%K==0){
                                      if (k==v*K+v){
                                          UFM[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
                                          BDRY[i+j*MM+k*MM*MM]=Ucoarse[Ds][Hs][Bs][lin];
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
    
    /*valuef.open("bdryfine.dat");
    
    
    if (valuef.is_open()) { // If file has correctly opened...
        // Output debug message
        cout << "Valuef correctly opened" << endl;
        
        // Dynamically store data into array
        for (int r=0; r<MM; r++) {
            for (int c=0; c<MM; c++){// ... and while there are no errors,
                for (int l=0; l<MM; l++){
                    valuef << UFM[r+c*MM+l*MM*MM] << " " ; // fill the row with col elements
                }
            }
            valuef << endl;
        }
    }
    else cout << "Unable to open valuef file" << endl;
    valuef.close();*/


    //fast sweep
    //s1=-1;
    //s2=1;
#pragma omp parallel for collapse(3)
    for (int m=0; m<=N-1; m++){ //should be parallel for loop
        for (int p=0; p<=N-1; p++){ //should be parallel for loop
           for (int v=0; v<=N-1; v++){
            //int full=1;
               
            for (int iter=0; iter<10;iter++){
                double maxError=0;
                 double currError;
              for (int s1=1; s1>=-1; s1-=2){
                 for (int s2=1; s2>=-1; s2-=2){
                     for (int s3=1; s3>=-1; s3-=2){
               //int s1=1;
                //int s2=1;
                //int s3=1;
                //s1=-1;
                 //s2=1;
                        for (int i =(s1<0 ? (p+1)*K+p : p*K+p); (s1<0 ? i>=p*K+p : i<= (p+1)*K+p); i+=s1){
                            for (int j=(s2<0 ? (m+1)*K+m : m*K+m); (s2<0 ? j>=m*K+m : j<= (m+1)*K+m); j+=s2)
                                for (int k =(s3<0 ? (v+1)*K+v : v*K+v); (s3<0 ? k>=v*K+v : k<= (v+1)*K+v); k+=s3){
                            {
                                if (obs[(i-p)+(j-m)*M+(k-v)*M*M]>1)
                                {
                                    //don't do anything; point is in obstacle
                                }
                                else if (BC[(i-p)+(j-m)*M+(k-v)*M*M]>=0)
                                {
                                    UFM[i+j*MM+k*MM*MM]=BC[(i-p)+(j-m)*M+(k-v)*M*M];
                                }
                                //else if ( BDRY[i+j*MM+k*MM*MM]>=0){
                                 //   UFM[i+j*MM+k*MM*MM]= BDRY[i+j*MM+k*MM*MM];
                                //}
                                
                                else if (BC[(i-p)+(j-m)*M+(k-v)*M*M]<0)// && BDRY[i+j*MM+k*MM*MM]<0)
                                {
                                    double m1,m2,m3,M1,M2,M3, unew, old, discr;
                                    int tempwindx,tempwindy,tempwindz, windtemp;
                                    if (j==0){
                                      m1=UFM[i+(j+1)*MM+k*MM*MM];
                                      tempwindy=7;
                                    }
                                    else if (j==MM-1){
                                        m1=UFM[i+(j-1)*MM+k*MM*MM];
                                        tempwindy=5;
                                    }
                                    else if (j==m*K+m){
                                        m1=UFM[i+(j+1)*MM+k*MM*MM];//(s2<0 ? UFM[i+(j+1)*MM] : 1000);//UFM[i+(j+1)*MM];//
                                        //tempwindy = (s2<0 ? 7 : 5);;
                                        tempwindy=7;
                                    }
                                    else if (j==(m+1)*K+m){
                                        m1=UFM[i+(j-1)*MM+k*MM*MM];//(s2<0 ? 1000 : UFM[i+(j-1)*MM]);//
                                        //tempwindy = (s2<0 ? 7 : 5);
                                        tempwindy=5;
                                    }
                                    else{
                                        m1=min(UFM[i+(j-1)*MM+k*MM*MM],UFM[i+(j+1)*MM+k*MM*MM]);//(s2>0 ? UFM[i+(j-1)*MM] : UFM[i+(j+1)*MM]);//
                                        if (UFM[i+(j-1)*MM+k*MM*MM]<UFM[i+(j+1)*MM+k*MM*MM]){//(s2>0){//(UFM[i+(j-1)*MM]<UFM[i+(j+1)*MM])//
                                            tempwindy=5;
                                        }
                                        else{
                                            tempwindy=7;
                                        }
                                    }
                                    if (i==0){
                                      m2=UFM[(i+1)+j*MM+k*MM*MM];
                                      tempwindx=8;
                                    }
                                    else if (i==MM-1){
                                      m2=UFM[(i-1)+j*MM+k*MM*MM];
                                      tempwindx=6;
                                    }
                                    else if (i==p*K+p){
                                        m2=UFM[(i+1)+j*MM+k*MM*MM];//(s1<0 ? UFM[(i+1)+j*MM] : 1000);//
                                        //tempwindx = (s1<0 ? 8 : 6);
                                        tempwindx=8;
                                    }
                                    else if (i==(p+1)*K+p){
                                        m2=UFM[(i-1)+j*MM+k*MM*MM];// ;(s1<0 ? 1000 : UFM[(i-1)+j*MM]);//
                                        //tempwindx = (s1<0 ? 8 : 6);
                                        tempwindx=6;
                                    }
                                    else{
                                        m2=min(UFM[(i-1)+j*MM+k*MM*MM], UFM[(i+1)+j*MM+k*MM*MM]);//(s1>0 ? UFM[(i-1)+j*MM] : UFM[(i+1)+j*MM]);//
                                        if (UFM[(i-1)+j*MM+k*MM*MM]<UFM[(i+1)+j*MM+k*MM*MM]){//(s1>0){//
                                            tempwindx=6;
                                        }
                                        else{
                                            tempwindx=8;
                                        }
                                    }
                                    
                                    if (k==0){
                                        m3=UFM[i+j*MM+(k+1)*MM*MM];
                                        tempwindx=8;
                                    }
                                    else if (k==MM-1){
                                        m3=UFM[i+j*MM+(k-1)*MM*MM];
                                        tempwindx=6;
                                    }
                                    else if (k==v*K+v){
                                        m3=UFM[i+j*MM+(k+1)*MM*MM];//(s1<0 ? UFM[(i+1)+j*MM] : 1000);//
                                        //tempwindx = (s1<0 ? 8 : 6);
                                        tempwindx=8;
                                    }
                                    else if (k==(v+1)*K+v){
                                        m3=UFM[i+j*MM+(k-1)*MM*MM];// ;(s1<0 ? 1000 : UFM[(i-1)+j*MM]);//
                                        //tempwindx = (s1<0 ? 8 : 6);
                                        tempwindx=6;
                                    }
                                    else{
                                        m3=min(UFM[i+j*MM+(k-1)*MM*MM], UFM[i+j*MM+(k+1)*MM*MM]);//(s1>0 ? UFM[(i-1)+j*MM] : UFM[(i+1)+j*MM]);//
                                        if (UFM[i+j*MM+(k-1)*MM*MM]<UFM[i+j*MM+(k+1)*MM*MM]){//(s1>0){//
                                            tempwindx=6;
                                        }
                                        else{
                                            tempwindx=8;
                                        }
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
                                    
                                     unew=M1+f[(i-p)+(j-m)*M+(k-v)*M*M]*h;
                                    
                                    if (unew <= M2){
                                         old=UFM[i+j*MM+k*MM*MM];
                                        if (old<unew){
                                            currError=0;
                                            UFM[i+j*MM+k*MM*MM]=old;
                                            maxError=max(currError,maxError);
                                        }
                                        else{
                                            
                                            WIND[i+j*MM+k*MM*MM]=windtemp;
                                            currError=abs(unew-old);
                                            UFM[i+j*MM+k*MM*MM]=unew;
                                            maxError=max(currError,maxError);
                                        }
                                        
                                    }
                                    else{
                                        discr=2*pow(f[(i-p)+(j-m)*M+(k-v)*M*M],2)*pow(h,2)-pow(M1-M2,2);
                                        unew=0.5*(M1+M2+sqrt(discr));
                                        if (unew <= M3){
                                           
                                            old=UFM[i+j*MM+k*MM*MM];
                                            if (old<unew){
                                                currError=0;
                                                UFM[i+j*MM+k*MM*MM]=old;
                                                maxError=max(currError,maxError);
                                            }
                                            else{
                                                
                                                WIND[i+j*MM+k*MM*MM]=windtemp;
                                                currError=abs(unew-old);
                                                UFM[i+j*MM+k*MM*MM]=unew;
                                                maxError=max(currError,maxError);
                                            }
                                           
                                        }
                                        else{
                                            discr=4*pow(M1+M2+M3,2)-12*(M1*M1+M2*M2+M3*M3-pow(f[(i-p)+(j-m)*M+(k-v)*M*M],2)*pow(h,2));
                                            unew=(2*(M1+M2+M3)+sqrt(discr))/6;
                                            if (unew<1000){
                                    
                                                old=UFM[i+j*MM+k*MM*MM];
                                                if (old<unew){
                                                    currError=0;
                                                    UFM[i+j*MM+k*MM*MM]=old;
                                                    maxError=max(currError,maxError);
                                                }
                                                else{
                                                    
                                                    WIND[i+j*MM+k*MM*MM]=windtemp;
                                                    currError=abs(unew-old);
                                                    UFM[i+j*MM+k*MM*MM]=unew;
                                                    maxError=max(currError,maxError);
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
                if (maxError<1e-12){
                    totalsweeps[p+m*(N)+v*N*N]=iter;
                    //full=0;
                    break;
                }
                else if (iter==9){
                    totalsweeps[p+m*(N)+v*N*N]=iter;
                    //full=0;
                }
                    }
          
                  }
                }
    }
    
  /*valuef.open("extfine.dat");
    
    
    if (valuef.is_open()) { // If file has correctly opened...
        // Output debug message
        cout << "Valuef correctly opened" << endl;
        
        // Dynamically store data into array
        for (int r=0; r<MM; r++) {
            for (int c=0; c<MM; c++){// ... and while there are no errors,
                for (int l=0; l<MM; l++){
                    valuef << UFM[r+c*MM+l*MM*MM] << " " ; // fill the row with col elements
                }
            }
            valuef << endl;
        }
    }
    else cout << "Unable to open valuef file" << endl;
    valuef.close();*/



#pragma omp parallel for collapse(3)
    for (int m=0; m<=N-1; m++){
        for (int p=0; p<=N-1; p++){
             for (int v=0; v<=N-1; v++){
            for (int i=p*K+p; i<=(p+1)*K+p; i++){
                for (int j=m*K+m; j<=(m+1)*K+m; j++){
                    for (int k=v*K+v; k<=(v+1)*K+v; k++)
                {
                    double old;
                    if (BC[(i-p)+(j-m)*M+(k-v)*M*M]>=0)
                    {
                        
                        FINE[(i-p)+(j-m)*M+(k-v)*M*M]=BC[(i-p)+(j-m)*M+(k-v)*M*M];
                        if (BC[(i-p)+(j-m)*M+(k-v)*M*M]>0 && BC[(i-p)+(j-m)*M+(k-v)*M*M]<1000){
                            //cout << "i = " << i << " j= " << j << " k= " << k << endl;
                            //cout << "BC= "<< BC[(i-p)+(j-m)*M+(k-v)*M*M] << endl;
                        }
                    }
                    else{
                        
                        /*if (i==0 || j==0 || k==0){
                            
                          
                            
                            old= FINE[(i-p)+(j-m)*M+(k-v)*M*M];
                            
                            
                            if (UFM[i+j*MM+k*MM*MM] < old)
                            {
                                FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                            }
                        }*/
                        
                       // else{
                        if (BDRY[i+j*MM+k*MM*MM]<0){
                            FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                        }
                        
                        
                        if  (i==0)
                            
                         {
                             if ((j-m)%K!=0 && (k-v)%K!=0)
                             {
                                 
                                 FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                             }
                             if ((j-m)%K==0){
                                 if ((j==(m+1)*K+m && (k-v)%K!=0) || (j==0 && k!=v*K+v)){
                                 FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                                 }
                             }
                             if ((k-v)%K==0){
                                 if  ((k==(v+1)*K+v && (j-m)%K!=0) || (k==0 && j!=m*K+m)){
                                     FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                                 }
                             }
                             if (k==(v+1)*K+v && j==(m+1)*K+m){
                                 FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                             }
                         }
                        if  (j==0)
                            {
                                if ((i-p)%K!=0 && (k-v)%K!=0)
                                {
                                    FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                                }
                                if ((i-p)%K==0){
                                    if  ((i==(p+1)*K+p && (k-v)%K!=0)  || (i==0 && k!=v*K+v)){
                                        FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                                    }
                                }
                                if ((k-v)%K==0){
                                    if ((k==(v+1)*K+v && (i-p)%K!=0) ||(k==0 && i!=p*K+p)){
                                        FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                                    }
                                }
                                
                                if (i==(p+1)*K+p && k==(v+1)*K+v){
                                     FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                                }
                            }
                        
                        
                        if  (k==0)
                        {
                            if ((j-m)%K!=0 && (i-p)%K!=0)
                            {
                                FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                            }
                            if ((j-m)%K==0){
                                if  ((j==(m+1)*K+m && (i-p)%K!=0) || (j==0 && i!=p*K+p)){
                                    FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                                }
                            }
                            if ((i-p)%K==0){
                                if  ((i==(p+1)*K+p && (j-m)%K!=0)|| (i==0 && j!=m*K+m)){
                                    FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                                }
                            }
                             if (i==(p+1)*K+p && j==(m+1)*K+m){
                                 FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                             }
                        }
                        
                        
                        //}
                        
                        /*if (BDRY[i+j*MM+k*MM*MM]<0){
                        old= FINE[(i-p)+(j-m)*M+(k-v)*M*M];
                        
                        
                        if (UFM[i+j*MM+k*MM*MM] < old)
                        {
                            FINE[(i-p)+(j-m)*M+(k-v)*M*M]=UFM[i+j*MM+k*MM*MM];
                        }
                        }
                        else{
                           if (FINE[(i-p)+(j-m)*M+(k-v)*M*M]==1000){
                                FINE[(i-p)+(j-m)*M+(k-v)*M*M]=BDRY[i+j*MM+k*MM*MM];
                            }
                        //}
                            //fwind[(i-p)+(j-m)*M+(k-v)*M*M]=WIND[i+j*MM+k*MM];
                        }*/
                        
                        //if (BDRY[i+j*MM]<0){
                        /*double old=fine[(i-p)+(j-m)*M+(k-v)*M*M];
                        int cont=1;
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
                        }*/
                       // double old;
                        //int cont;
                        
                        /*if (i==p*K+p && j==m*K+m){
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
                        }*/
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
