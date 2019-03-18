//
//
//  Created by Lindsay Martin on 9/23/18.
//
//
#include <iostream>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <stdio.h>
using namespace std;

#include <fstream>


#include "sweepdefs.h"

#define COLS 51

#define ROWS 51
#define HEIGHT 51
#define DIM1 10
#define DIM2 10
#define DIM3 10

int main(int argc, char *argv[])
{
    if (argc<7){
        printf("usage: main [N] [K] [T0] [Gamma] [delta] [x0] [paraiter]\n");
        return 0;
    }        //declare variables
    
    int N = atoi(argv[1]);
    int K = atoi(argv[2]);
    //int  nthreads = atoi(argv[3]);
    //int numTrials = 5000/(N*K+1);
    double tic, toc;
    int M = N*K+1;
    double  *BC, *F,*ufine, *Ufine, *obs,  *w, *fw, *fwind, *UNEW, *UFM,  *UF1, *UF2, *overallfine;//, *Ucoarse1;
    double h,H,C, T0 =atof(argv[3]), eps=atof(argv[4]), delta=atof(argv[5]), x0=atof(argv[6]);
    int paraiter=atoi(argv[7]);
    int sweepiter;
    //Ucoarse1 =new double [M*M];
    // Set number of threads.
    int nthreads=1;
    omp_set_num_threads(nthreads);
    
    printf("\n N= %d, K= %d, Threads: %d \n T0= %4.3f, eps= %4.3f, delta= %4.3f, x0= %4.3f\n", N, K,nthreads, T0, eps, delta, x0);

    //N  number of nodes on coarse grid
    //K  number of nodes on fine grid
    
    sweepiter= 5;


    H=1.0/(double) N;
    h=1.0/(double)(N*K);
    C=0.5;
    
    //associate pointers for inputs

    F=new double [M*M*M];
    BC=new double [M*M*M];
    obs=new double [M*M*M];
    overallfine=new double[M*M*M];
    
    UF2=new double [M*M*M];
    UF1=new double [M*M*M];
    UFM=new double [M*M*M];
    
    
    fwind=new double [M*M*M];
    
    //read in Marmousi velocity
    
    
    
    //serial FSM
    #pragma omp parallel for
    for (int c=0; c<M; c++)
        for (int b=0; b<M; b++)
            for (int a=0; a<M; a++){
                //U[c]=1000;
                overallfine[c+b*M+a*M*M]=1000;
                obs[c+b*M+a*M*M]=0;
                BC[c+b*M+a*M*M]=-1;
		double ep=50;
                F[c+b*M+a*M*M]=1;//+0.01*sin(M_PI*c*h*ep)*sin(M_PI*b*h*ep)*sin(M_PI*a*h*ep);
            	 double x,y,z;
                x=c*h;
                y=b*h;
                z=a*h;
	 	double r=sqrt(pow(x-.3,2)+pow(y-.3,2)+pow(z-.3,2));
		if (r<.025){
                     F[c+b*M+a*M*M]=.001;
                }
                
                 r=sqrt(pow(x-.7,2)+pow(y-.7,2)+pow(z-.5,2));
                if (r<.025){
                    F[c+b*M+a*M*M]=10;
                }
                
                 r=sqrt(pow(x-.5,2)+pow(y-.4,2)+pow(z-.4,2));
                if (r<.025){
                    F[c+b*M+a*M*M]=.001;
                }
                
                r=sqrt(pow(x-.6,2)+pow(y-.3,2)+pow(z-.6,2));
                if (r<.025){
                    F[c+b*M+a*M*M]=10;
                }
		}
    
    int  k1= 0;
    int  m1=0;
    int  l1=0;
    BC[k1+m1*M+l1*M*M]=0;
    int totaliters;
    tic = omp_get_wtime();
    
    totaliters=fsm(overallfine, F, BC, obs, 4, M, h);
    
    toc = (omp_get_wtime() - tic);
    
  printf("Iters for convergence:  %d \n", totaliters);
   printf("FSM time:  %4.8f \n", toc);
    //output value function
    //uncomment to output value function to file valuef.dat
    ofstream valuef;
    int r=0;
    /*valuef.open("overallfine.dat");
    
    
    if (valuef.is_open()) { // If file has correctly opened...
        // Output debug message
        cout << "Valuef correctly opened" << endl;
        
        // Dynamically store data into array
        for (int r=0; r<M; r++) {
            for (int c=0; c<M; c++){// ... and while there are no errors,
                for (int l=0; l<M; l++){
                    valuef << overallfine[r+c*M+l*M*M] << " " ; // fill the row with col elements
                }
            }
            valuef << endl;
        }
    }
    else cout << "Unable to open valuef file" << endl;
    valuef.close();
    */
    
    //NEW METHOD
    for (int threads=1; threads<=64; threads*=2){
        omp_set_num_threads(threads);
    
    int disp=0;
    //int k,c;
    double* CG[DIM1][DIM2][DIM3];
    double *UK[DIM1][DIM2][DIM3];
    double *UK1[DIM1][DIM2][DIM3];
    double *UK2[DIM1][DIM2][DIM3];
    double *UK3[DIM1][DIM2][DIM3];
     double *wind[DIM1][DIM2][DIM3];
    double *Ucoarse1[DIM1][DIM2][DIM3];
    double *theta[DIM1][DIM2][DIM3];

    for (int D=0; D<DIM1; D++)
        for (int H=0; H<DIM2; H++)
            for (int B=0; B<DIM3; B++)
        {      CG[D][H][B]=new double[(N+1)*(N+1)*(N+1)];
            UK[D][H][B]=new double[(N+1)*(N+1)*(N+1)];
            UK1[D][H][B]=new double[(N+1)*(N+1)*(N+1)];
            UK2[D][H][B]=new double[(N+1)*(N+1)*(N+1)];
            UK3[D][H][B]=new double[(N+1)*(N+1)*(N+1)];
            wind[D][H][B]=new double[(N+1)*(N+1)*(N+1)];
            Ucoarse1[D][H][B]=new double[(N+1)*(N+1)*(N+1)];
             theta[D][H][B]=new double[(N+1)*(N+1)*(N+1)];
            
            for(int k=0; k<(N+1)*(N+1)*(N+1); k++)
                CG[D][H][B][k]= 1000;
        }
    
     double ticcoarse;
    double ticfine;
    double toccoarse;
    double tocfine;
    double sumcoarsetime=0;
    double sumfinetime=0;
   #pragma omp parallel for
    for (int k=0; k<=k1+K; k++){
        for (int c=0; c<=m1+K; c++){
            for (int q=0; q<=l1+K; q++){
                BC[k+c*M+q*M*M]=overallfine[k+c*M+q*M*M];
            }
        }
    }

    //int k,c;
    int sweepcoarse[1];
    int totalsweepcoarse[3*K*(K-1)+1];
    int coarsecounter=0;
    int finecounter=0;
    int sumfinesweeps=0;
    ticcoarse = omp_get_wtime();
    
    coarseinitial(CG,F,BC,obs,N,K,0,0,0,wind,sweepiter, sweepcoarse);
    totalsweepcoarse[0]=sweepcoarse[0];
    
    
    
    
   
    
    //k=0
 
 #pragma omp parallel for
   for (int p=1;p<3*K-2;p++)
    {
        if (p<K){
            int sweepcoarse[1];
            coarseinitial(CG,F,BC,obs,N,K,0,p,0,wind, sweepiter,sweepcoarse);
             totalsweepcoarse[p]=sweepcoarse[0];
        }
        else if(p<2*K-1){
            int sweepcoarse[1];
            coarseinitial(CG,F,BC,obs,N,K,p-(K-1),0,0,wind, sweepiter,sweepcoarse);
            
            totalsweepcoarse[p]=sweepcoarse[0];
        }
        else{
             int sweepcoarse[1];
            coarseinitial(CG,F,BC,obs,N,K,0,0,p-(2*K-2),wind, sweepiter,sweepcoarse);
           
             totalsweepcoarse[p]=sweepcoarse[0];
        }
        
    }
    #pragma omp parallel for collapse (2)
    for (int p=1; p<K; p++)
        for (int r=1; r<K; r++)
        {
            int sweepcoarse[1];
            coarseinitial(CG,F,BC,obs,N,K,0,p,r,wind, sweepiter,sweepcoarse);
            totalsweepcoarse[p]=sweepcoarse[0];
        }
    
#pragma omp parallel for collapse (2)
    for (int p=1; p<K; p++)
        for (int r=1; r<K; r++)
        {
            int sweepcoarse[1];
            coarseinitial(CG,F,BC,obs,N,K,r,p,0,wind, sweepiter,sweepcoarse);
            totalsweepcoarse[p]=sweepcoarse[0];
        }
    
#pragma omp parallel for collapse (2)
    for (int p=1; p<K; p++)
        for (int r=1; r<K; r++)
        {
            int sweepcoarse[1];
            coarseinitial(CG,F,BC,obs,N,K,r,0,p,wind, sweepiter,sweepcoarse);
            totalsweepcoarse[p]=sweepcoarse[0];
        }
  
  
  

 
    
    
    
  toccoarse = (omp_get_wtime() - ticcoarse);
    sumcoarsetime+=toccoarse;
    int sumcoarsesweeps=0;
    for (int i=0; i<3*K-2; i++){
        sumcoarsesweeps+=totalsweepcoarse[i];
    }
    coarsecounter+=3*K-2;
    ticcoarse=omp_get_wtime();
    
  
    
   // cout << "CG[1][1][11^3]= " << CG[2][1][11*11*11-1] << endl;
    
   

    #pragma omp parallel for
    for (int D=0; D<DIM1; D++)
        for (int H=0; H<DIM2; H++)
            for (int B=0; B<DIM3; B++)
                for (int k=0; k<(N+1)*(N+1)*(N+1); k++){
                    UK1[D][H][B][k]=CG[D][H][B][k];
                     UK3[D][H][B][k]=CG[D][H][B][k];
                }

    
    
    //causalsweeping(CG, BC, obs, N, K, wind,F);
   

   toccoarse = (omp_get_wtime() - ticcoarse);
    

    ticfine= omp_get_wtime();
    #pragma omp parallel for
    for (int c=0; c<M*M*M; c++){
        UFM[c]=1000;
    }
    int *sweepfine;
    sweepfine = new int [N*N*N];
    fineupdate(UFM,CG, F,BC,obs,N,K,wind,fwind,sweepfine);
   
   
   
    
     //cout << "UFM[51,51,52]= " << UFM[50+50*M+51*M*M] << endl;
    tocfine=omp_get_wtime()-ticfine;
    sumfinetime+=tocfine;
    
    
    for (int p=0; p<N*N*N; p++){
        sumfinesweeps+=sweepfine[p];
 	}


    finecounter+=N*N;
    ticfine=omp_get_wtime();
    
    
    #pragma omp parallel for
    for (int c=0; c<M*M*M; c++){
        UF1[c]=UFM[c];
    }

    tocfine = (omp_get_wtime() - ticfine);
    sumfinetime+=tocfine;
    
    //calculate relative L1 error
    double diff=0;
    
    for (int v=0; v<M*M*M; v++){
            diff=diff + fabs(UFM[v]-overallfine[v]);
        }
        
    
    
    diff=diff/double (M*M*M);
     //cout  << diff << endl;
   
    
    //k=1
    
    if (paraiter>=1){
        ticcoarse = omp_get_wtime();
      
        #pragma omp parallel for
        for (int D=0; D<DIM1; D++)
            for (int H=0; H<DIM2; H++)
                for (int B=0; B<DIM3; B++)
                     for(int k=0; k<(N+1)*(N+1)*(N+1); k++)
                     {
                         CG[D][H][B][k]= 1000;
                     }
        
        
        coarseupdate(CG, UFM, UK1, F,BC,obs,N,K,0,0,0,wind, fwind,T0, sweepiter,sweepcoarse);
        totalsweepcoarse[0]=sweepcoarse[0];
        
        
        
       
        
        
#pragma omp parallel for
        for (int p=1;p<3*K-2;p++)
        {
            if (p<K){
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,0,p,0,wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            else if(p<2*K-1){
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,p-(K-1),0,0,wind, fwind,T0, sweepiter,sweepcoarse);
                
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            else{
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,0,0,p-(2*K-2),wind, fwind,T0, sweepiter,sweepcoarse);
                
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            
        }
#pragma omp parallel for collapse (2)
        for (int p=1; p<K; p++)
            for (int r=1; r<K; r++)
            {
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,0,p,r,wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
        
#pragma omp parallel for collapse (2)
        for (int p=1; p<K; p++)
            for (int r=1; r<K; r++)
            {
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,r,p,0,wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
        
#pragma omp parallel for collapse (2)
        for (int p=1; p<K; p++)
            for (int r=1; r<K; r++)
            {
                int sweepcoarse[1];
               coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,r,0,p,wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
        
     
      
        
        
       
        toccoarse = (omp_get_wtime() - ticcoarse);
        sumcoarsetime+=toccoarse;
        
        for (int i=0; i<3*K-2; i++){
            sumcoarsesweeps+=totalsweepcoarse[i];
        }
        coarsecounter+=3*K-2;
        ticcoarse=omp_get_wtime();
    
       
        #pragma omp parallel for
        for (int D=0; D<DIM1; D++)
            for (int H=0; H<DIM2; H++)
                for (int B=0; B<DIM3; B++)
                    for (int k=0; k<(N+1)*(N+1)*(N+1); k++)
                {
                    UK2[D][H][B][k]=UK1[D][H][B][k];
                }
    
    
        //causalsweeping(CG, BC, obs, N, K, wind,F);
        toccoarse = (omp_get_wtime() - ticcoarse);
        sumcoarsetime+=toccoarse;
       
        ticfine=omp_get_wtime();
        #pragma omp parallel for
        for (int c=0; c<M*M*M; c++){
            UFM[c]=1000;
        }
        
        fineupdate(UFM,CG, F,BC,obs,N,K,wind,fwind,sweepfine);
        tocfine=omp_get_wtime()-ticfine;
        sumfinetime+=tocfine;
      
        
        
        for (int p=0; p<N*N*N; p++){
            sumfinesweeps+=sweepfine[p];
        }
        finecounter+=N*N;
        ticfine=omp_get_wtime();
       
        #pragma omp parallel for
        for (int c=0; c<M*M*M; c++){
            UF2[c]=UF1[c];
            UF1[c]=UFM[c];
        }
        tocfine=(omp_get_wtime()-ticfine);
        sumfinetime+=tocfine;
         diff=0;
       
        for (int v=0; v<M*M*M; v++){
                diff=diff + fabs(UFM[v]-overallfine[v]);
            }
            
       
        
        diff=diff/double (M*M*M);
       //  cout << diff << endl;
      
    }
    //k=2
    if (paraiter>=2){
        
        
        ticcoarse=omp_get_wtime();
        #pragma omp parallel for
        for (int D=0; D<DIM1; D++)
            for (int H=0; H<DIM2; H++)
                for (int B=0; B<DIM3; B++)
                    for(int k=0; k<(N+1)*(N+1)*(N+1); k++)
                    {
                        CG[D][H][B][k]= 1000;
                    }
        
        
        coarseupdate(CG, UFM, UK1, F,BC,obs,N,K,0,0,0,wind, fwind,T0, sweepiter,sweepcoarse);
        totalsweepcoarse[0]=sweepcoarse[0];
        
        
       
       // cout << "after k=2 0,0 coarseupdate"<< endl;
        
        
#pragma omp parallel for
        for (int p=1;p<3*K-2;p++)
        {
            if (p<K){
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,0,p,0,wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            else if(p<2*K-1){
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,p-(K-1),0,0,wind, fwind,T0, sweepiter,sweepcoarse);
                
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            else{
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,0,0,p-(2*K-2),wind, fwind,T0, sweepiter,sweepcoarse);
                
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            
        }
#pragma omp parallel for collapse (2)
        for (int p=1; p<K; p++)
            for (int r=1; r<K; r++)
            {
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,0,p,r,wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
        
#pragma omp parallel for collapse (2)
        for (int p=1; p<K; p++)
            for (int r=1; r<K; r++)
            {
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,r,p,0,wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
        
#pragma omp parallel for collapse (2)
        for (int p=1; p<K; p++)
            for (int r=1; r<K; r++)
            {
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,r,0,p,wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
        
       
         
        toccoarse = (omp_get_wtime() - ticcoarse);
        sumcoarsetime+=toccoarse;
        
        for (int i=0; i<3*K-2; i++){
            sumcoarsesweeps+=totalsweepcoarse[i];
        }
        coarsecounter+=3*K-2;
        ticcoarse=omp_get_wtime();

       #pragma omp parallel for
        for (int D=0; D<DIM1; D++)
            for (int H=0; H<DIM2; H++)
                for (int B=0; B<DIM3; B++)
                    for (int k=0; k<(N+1)*(N+1)*(N+1); k++)
                    {
                    Ucoarse1[D][H][B][k]=CG[D][H][B][k];
                    }
        
        
        //causalsweeping(CG, BC, obs, N, K, wind,F);
 
        
        toccoarse=(omp_get_wtime()-ticcoarse);
        sumcoarsetime+=toccoarse;
        
        
        ticfine= omp_get_wtime();

        #pragma omp parallel for
        for (int c=0; c<M*M*M; c++){
            UFM[c]=1000;
        }
        
        fineupdate(UFM,CG, F,BC,obs,N,K,wind,fwind,sweepfine);
        tocfine=omp_get_wtime()-ticfine;
        sumfinetime+=tocfine;
     
        
        for (int p=0; p<N*N; p++){
            sumfinesweeps+=sweepfine[p];
        }
        finecounter+=N*N*N;
        ticfine=omp_get_wtime();
        #pragma omp parallel for
        for (int c=0; c<M*M*M; c++){
            UF2[c]=UF1[c];
            UF1[c]=UFM[c];
        }
        tocfine=(omp_get_wtime()-ticfine);
        sumfinetime+=tocfine;
         diff=0;
        for (int v=0; v<M*M*M; v++){
            
                diff=diff + fabs(UFM[v]-overallfine[v]);
            }
            
       
        
        diff=diff/double (M*M*M);
        // cout << diff << endl;
    }
    
    
    
    ////////////////ALGORITHM////////////////
    totaliters=2;
   
    for (int a=3; a<=paraiter; a++){
        totaliters++;
        ticcoarse=omp_get_wtime();
        #pragma omp parallel for
        for (int D=0; D<DIM1; D++)
            for (int H=0; H<DIM2; H++)
                for (int B=0; B<DIM3; B++)
                    for(int k=0; k<(N+1)*(N+1)*(N+1); k++)
                    {
                        CG[D][H][B][k]= 1000;
                    }
        
        int w1=3,w2=1,w3=1;
       
        coarseupdatenewtheta(CG, Ucoarse1,UF1,UF2,UK1,UK2,UK3, F,BC,obs,N,K,0,0,0,wind, fwind, sweepiter, eps, delta, x0,theta,w1,w2,w3,sweepcoarse);
        totalsweepcoarse[0]=sweepcoarse[0];
        
       
       
        

        
#pragma omp parallel for
        for (int p=1;p<3*K-2;p++)
        {
            if (p<K){
                int sweepcoarse[1];
                 coarseupdatenewtheta(CG, Ucoarse1,UF1,UF2, UK1,UK2,UK3, F,BC,obs,N,K,0,p,0,wind, fwind, sweepiter, eps, delta, x0,theta,w1,w2,w3,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            else if(p<2*K-1){
                int sweepcoarse[1];
                coarseupdatenewtheta(CG, Ucoarse1,UF1,UF2, UK1,UK2,UK3, F,BC,obs,N,K,p-(K-1),0,0,wind, fwind, sweepiter, eps, delta, x0,theta,w1,w2,w3,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            else{
                int sweepcoarse[1];
               coarseupdatenewtheta(CG, Ucoarse1,UF1,UF2, UK1,UK2,UK3, F,BC,obs,N,K,0,0,p-(2*K-2),wind, fwind, sweepiter, eps, delta, x0,theta,w1,w2,w3,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            
        }
#pragma omp parallel for collapse (2)
        for (int p=1; p<K; p++)
            for (int r=1; r<K; r++)
            {
                int sweepcoarse[1];
                 coarseupdatenewtheta(CG, Ucoarse1,UF1,UF2, UK1,UK2,UK3, F,BC,obs,N,K,0,p,r,wind, fwind, sweepiter, eps, delta, x0,theta,w1,w2,w3,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
        
#pragma omp parallel for collapse (2)
        for (int p=1; p<K; p++)
            for (int r=1; r<K; r++)
            {
                int sweepcoarse[1];
                 coarseupdatenewtheta(CG, Ucoarse1,UF1,UF2, UK1,UK2,UK3, F,BC,obs,N,K,r,p,0,wind, fwind, sweepiter, eps, delta, x0,theta,w1,w2,w3,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
        
#pragma omp parallel for collapse (2)
        for (int p=1; p<K; p++)
            for (int r=1; r<K; r++)
            {
                int sweepcoarse[1];
                coarseupdatenewtheta(CG, Ucoarse1,UF1,UF2, UK1,UK2,UK3, F,BC,obs,N,K,r,0,p,wind, fwind, sweepiter, eps, delta, x0,theta,w1,w2,w3,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
        
        
        
        toccoarse = (omp_get_wtime() - ticcoarse);
        sumcoarsetime+=toccoarse;
        
        for (int i=0; i<3*K-2; i++){
            sumcoarsesweeps+=totalsweepcoarse[i];
        }
        coarsecounter+=3*K-2;
        ticcoarse=omp_get_wtime();
        
       // causalsweeping(CG, BC, obs, N, K, wind,F);
        toccoarse=omp_get_wtime()-ticcoarse;
        sumcoarsetime+=toccoarse;
        
        ticfine=omp_get_wtime();
        
        /*valuef.open("coarsenewtheta.dat");
        if (valuef.is_open()) { // If file has correctly opened...
            // Output debug message
            cout << "Valuef correctly opened" << endl;
            
            // Dynamically store data into array
            for (int r=0; r<N+1; r++) {
                for (int c=0; c<N+1; c++){// ... and while there are no errors,
                    for (int l=0; l<N+1; l++){
                        valuef << CG[0][1][1][r+c*(N+1)+l*(N+1)*(N+1)] << " " ; // fill the row with col elements
                    }
                }
                valuef << endl;
            }
        }
        else cout << "Unable to open valuef file" << endl;
        valuef.close();*/
        
    #pragma omp parallel for
        for (int c=0; c<M*M*M; c++){
            UFM[c]=1000;
        }
        
        fineupdate(UFM,CG, F,BC,obs,N,K,wind,fwind,sweepfine);
        
        tocfine=omp_get_wtime()-ticfine;
        sumfinetime+=tocfine;
        
        
        for (int p=0; p<N*N*N; p++){
            sumfinesweeps+=sweepfine[p];
        }
        finecounter+=N*N;
        ticfine=omp_get_wtime();
        
#pragma omp parallel for
        for (int c=0; c<M*M*M; c++){
            UF2[c]=UF1[c];
            UF1[c]=UFM[c];
        }
        tocfine=(omp_get_wtime()-ticfine);
        sumfinetime+=tocfine;
        
         diff=0;
        for (int v=0; v<M*M*M; v++){
                diff=diff + fabs(UFM[v]-overallfine[v]);
            }
            
        
        
        diff=diff/double (M*M*M);
        
        
        //cout << diff << endl;
        
        /*if (diff<1e-5 && disp==0){
            //cout<< "fineupdate time: " <<sumtime << endl;
            // break;
            
            double totaltime= sumcoarsetime+sumfinetime;
            printf("Parallel For total time:        %4.8fs \n", totaltime);
            printf("Parallel For coarse time:        %4.8fs \n", sumcoarsetime);
            printf("Parallel For fine time:        %4.8fs \n", sumfinetime);
            printf("Parallel For Total iters:        %d \n", a);
            cout << "L1 relative error: " << diff << endl;
            disp=1;
            cout << "Average coarse sweeps: " <<  double (sumcoarsesweeps)/double (coarsecounter) << endl;
            cout << "Average fine sweeps: " <<  double (sumfinesweeps)/double (finecounter) << endl;
        }*/
        /*if (diff<1e-13){
            //cout<< "fineupdate time: " <<sumtime << endl;
            // break;
            
            double totaltime= sumcoarsetime+sumfinetime;
            printf("Parallel For total time:        %4.8fs \n", totaltime);
            printf("Parallel For coarse time:        %4.8fs \n", sumcoarsetime);
            printf("Parallel For fine time:        %4.8fs \n", sumfinetime);
            printf("Parallel For Total iters:        %d \n", a);
            cout << "L1 relative error: " << diff << endl;
            cout << "Average coarse sweeps: " <<  double (sumcoarsesweeps)/double (coarsecounter) << endl;
            cout << "Average fine sweeps: " <<  double (sumfinesweeps)/double (finecounter) << endl;
            break;
        }*/
        
        
        
    }
        double totaltime= sumcoarsetime+sumfinetime;
        cout << totaltime << endl;
   
    
   /* valuef.open("valuef.dat");
     
     
     if (valuef.is_open()) { // If file has correctly opened...
     // Output debug message
     cout << "Valuef correctly opened" << endl;
     
     // Dynamically store data into array
     for (int r=0; r<M; r++) {
     for (int c=0; c<M; c++){// ... and while there are no errors,
     for (int l=0; l<M; l++){
     valuef << UFM[r+c*M+l*M*M] << " " ; // fill the row with col elements
     }
     }
     valuef << endl;
     }
     }
     else cout << "Unable to open valuef file" << endl;
     valuef.close();
    */
    
    //output value function
    //uncomment to output value funcntion to file valueftheta.dat
   /* ofstream valueftheta;
   // r=0;
    valueftheta.open("valueftheta.dat");
    
    
    if (valueftheta.is_open()) { // If file has correctly opened...
        // Output debug message
        cout << "Valueftheta correctly opened" << endl;
        
        // Dynamically store data into array
        for (int c=0; c<(N+1); c++) {
            for (int r=0; r<(N+1); r++){// ... and while there are no errors,
                valueftheta << CG[0][0][c+r*(N+1)] << " " ; // fill the row with col elements
            }
            valueftheta << endl;
        }
    }
    else cout << "Unable to valueftheta file" << endl;
    valueftheta.close();*/
    
   
    
    
    for (int D=0; D<DIM1; D++)
        for (int H=0; H<DIM2; H++)
            for (int B=0; B<DIM3; B++)
        {
            delete[] CG[D][H][B];
            delete[] UK[D][H][B];
            delete[] UK1[D][H][B];
            delete[] UK2[D][H][B];
            delete[] UK3[D][H][B];
            delete[] wind[D][H][B];
            delete[] Ucoarse1[D][H][B];
            delete[] theta[D][H][B];
            
        }
    

     delete[] sweepfine;
   
   
    //delete[] prev;
   
    }
    
    delete[] overallfine;
    //delete[] U;
    delete[] F;
    delete[] BC;
    delete[] obs;
   
    
    
    delete[] UF2;
    delete[] UF1;
    
    delete[] UFM;
    
    delete[] fwind;
    
     return 0;
}

