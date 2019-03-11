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

#include "marmtest.h"
#include "sweepdefs.h"


#define filename "marminterp.dat"
#define COLS 11201
#define ROWS 11201
#define DIM1 2
#define DIM2 800

int main(int argc, char *argv[])
{
    if (argc<7){
        printf("usage: main [N] [K] [num threads] [T0] [Gamma] [delta] [x0] \n");
        return 0;
    }        //declare variables
    
    int N = atoi(argv[1]);
    int K = atoi(argv[2]);
    int  nthreads = atoi(argv[3]);
    //int numTrials = 5000/(N*K+1);
    double tic, toc;
    int M = N*K+1;
    double  *BC, *F,*ufine, *Ufine, *obs,  *w, *fw, *fwind, *UNEW, *UFM,  *UF1, *UF2, *overallfine;//, *Ucoarse1;
    double h,H,C, T0 =atof(argv[4]), eps=atof(argv[5]), delta=atof(argv[6]), x0=atof(argv[7]);
    int paraiter,sweepiter;
    //Ucoarse1 =new double [M*M];
    // Set number of threads.
    omp_set_num_threads(nthreads);
    
    printf("\n N= %d, K= %d, Threads: %d \n T0= %4.3f, eps= %4.3f, delta= %4.3f, x0= %4.3f\n", N, K,nthreads, T0, eps, delta, x0);

    //N  number of nodes on coarse grid
    //K  number of nodes on fine grid
    paraiter=50;
    sweepiter= 10;


    H=1.0/(double) N;
    h=1.0/(double)(N*K);
    C=0.5;
    
    //associate pointers for inputs

    F=new double [M*M];
    BC=new double [M*M];
    obs=new double [M*M];
    overallfine=new double[M*M];
    
    UF2=new double [M*M];
    UF1=new double [M*M];
    UFM=new double [M*M];
    
    
    fwind=new double [M*M];
    
    //read in Marmousi velocity
    
    double *marmarray;
    marmarray = new double [ROWS*COLS];
    
    for (int k=0; k<ROWS*COLS; k++){
        marmarray[k]=1;
    }
    
    fstream marmdat;
    
    marmdat.open(filename, ios::in);
    
    if (marmdat.is_open()){
        marmread(marmarray, ROWS, COLS, marmdat);
        
        marmdat.close();
    }
    else{
        cout << "ERROR reading marmdat" << endl;
        return 0;
    }
    
    //serial FSM
    #pragma omp parallel for
    for (int c=0; c<M*M; c++){
        //U[c]=1000;
        overallfine[c]=1000;
        obs[c]=0;
        BC[c]=-1;
        F[c]=1.0/(double) marmarray[c];
        
    }
    
    int  k1= 2400;
    int  m1=1600;
    BC[k1+m1*M]=0;
    int totaliters;
    tic = omp_get_wtime();
    
    totaliters=fsm(overallfine, F, BC, obs, 100, M, h);
    
    toc = (omp_get_wtime() - tic);
    
    printf("Iters for convergence:  %d \n", totaliters);
    printf("FSM time:  %4.8f \n", toc);
    //output value function
    //uncomment to output value function to file valuef.dat
   /*ofstream valuef;
    int r=0;
    valuef.open("valuef.dat");
    
    
    if (valuef.is_open()) { // If file has correctly opened...
        // Output debug message
        cout << "Valuef correctly opened" << endl;
        
        // Dynamically store data into array
        for (int r=0; r<ROWS; r++) {
            for (int c=0; c<COLS; c++){// ... and while there are no errors,
                valuef << overallfine[r+c*COLS] << " " ; // fill the row with col elements
            }
            valuef << endl;
        }
    }
    else cout << "Unable to valuef file" << endl;
    valuef.close();*/
    
    
    //NEW METHOD
  
    int disp=0;
    //int k,c;
    double* CG[DIM1][DIM2];
    double *UK[DIM1][DIM2];
    double *UK1[DIM1][DIM2];
    double *UK2[DIM1][DIM2];
    double *UK3[DIM1][DIM2];
     double *wind[DIM1][DIM2];
    double *Ucoarse1[DIM1][DIM2];
    double *theta[DIM1][DIM2];

    for (int dir=0; dir<2; dir++)
        for (int i=0; i<K; i++)
        {      CG[dir][i]=new double[(N+1)*(N+1)];
            UK[dir][i]=new double[(N+1)*(N+1)];
            UK1[dir][i]=new double[(N+1)*(N+1)];
            UK2[dir][i]=new double[(N+1)*(N+1)];
            UK3[dir][i]=new double[(N+1)*(N+1)];
            wind[dir][i]=new double[(N+1)*(N+1)];
            Ucoarse1[dir][i]=new double[(N+1)*(N+1)];
             theta[dir][i]=new double[(N+1)*(N+1)];
            
            for(int k=0; k<(N+1)*(N+1); k++)
                CG[dir][i][k]= 1000;
        }
    
    double ticcoarse;
    double ticfine;
    double toccoarse;
    double tocfine;
    double sumcoarsetime=0;
    double sumfinetime=0;
   #pragma omp parallel for
    for (int k=k1-K; k<=k1+K; k++){
        for (int c=m1-K; c<=m1+K; c++){
            BC[k+c*M]=overallfine[k+c*M];
        }
    }

    //int k,c;
    int sweepcoarse[1];
    int totalsweepcoarse[2*K-1];
    int coarsecounter=0;
    int finecounter=0;
    int sumfinesweeps=0;
    ticcoarse = omp_get_wtime();
    
    coarseinitial(CG,F,BC,obs,N,K,0,0,wind,sweepiter, sweepcoarse);
    totalsweepcoarse[0]=sweepcoarse[0];
   
    
    //k=0
 
    #pragma omp parallel for
   for (int p=1;p<2*K-1;p++)
    {
        if (p<K){
            int sweepcoarse[1];
            coarseinitial(CG,F,BC,obs,N,K,0,p,wind, sweepiter,sweepcoarse);
             totalsweepcoarse[p]=sweepcoarse[0];
        }
        else{
             int sweepcoarse[1];
            coarseinitial(CG,F,BC,obs,N,K,p-(K-1),0,wind, sweepiter,sweepcoarse);
           
             totalsweepcoarse[p]=sweepcoarse[0];
        }
        
    }
   toccoarse = (omp_get_wtime() - ticcoarse);
    sumcoarsetime+=toccoarse;
    int sumcoarsesweeps=0;
    for (int i=0; i<2*K-1; i++){
        sumcoarsesweeps+=totalsweepcoarse[i];
    }
    coarsecounter+=2*K-1;
    ticcoarse=omp_get_wtime();


    #pragma omp parallel for
    for (int dir=0; dir<2; dir++){
        for (int i=0; i<K; i++){
            for (int k=0; k<(N+1)*(N+1); k++){
                UK1[dir][i][k]=CG[dir][i][k];
                 UK3[dir][i][k]=CG[dir][i][k];
            }
        }
    }
    
    
    causalsweeping(CG, BC, obs, N, K, wind,F);
   

   toccoarse = (omp_get_wtime() - ticcoarse);
    

    ticfine= omp_get_wtime();
    #pragma omp parallel for
    for (int c=0; c<M*M; c++){
        UFM[c]=1000;
    }
    int *sweepfine;
    sweepfine = new int [N*N];
    fineupdate(UFM,CG, F,BC,obs,N,K,wind,fwind,sweepfine);
   
    
    tocfine=omp_get_wtime()-ticfine;
    sumfinetime+=tocfine;
    
    
    for (int p=0; p<N*N; p++){
        sumfinesweeps+=sweepfine[p];
 	}


    finecounter+=N*N;
    ticfine=omp_get_wtime();
    
    
    #pragma omp parallel for
    for (int c=0; c<M*M; c++){
        UF1[c]=UFM[c];
    }

    tocfine = (omp_get_wtime() - ticfine);
    sumfinetime+=tocfine;
    
    //calculate relative L1 error
    double diff=0;
    
    for (int v=0; v<M*M; v++){
        for (int r=0; r< M*M; r++){
            diff=diff + fabs(UFM[v+r*M]-overallfine[v+r*M]);
        }
        
    }
    
    diff=diff/double (M*M);
   
    
    //k=1
    
    if (paraiter>=1){
        ticcoarse = omp_get_wtime();
      
        #pragma omp parallel for
        for (int dir=0; dir<2; dir++){
            for (int i=0; i<K; i++){
                for (int k=0; k<(N+1)*(N+1); k++){
                    CG[dir][i][k]=1000;
                }
            }
        }
        
        
        coarseupdate(CG, UFM, UK1, F,BC,obs,N,K,0,0,wind, fwind,T0, sweepiter,sweepcoarse);
        totalsweepcoarse[0]=sweepcoarse[0];
      
        
     
        #pragma omp parallel for
        for (int p=1;p<2*K-1;p++)
        {
            
            if (p<K){
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,p,0,wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            else{
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,0,p-(K-1),wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            
        }
        toccoarse = (omp_get_wtime() - ticcoarse);
        sumcoarsetime+=toccoarse;
        
        for (int i=0; i<2*K-1; i++){
            sumcoarsesweeps+=totalsweepcoarse[i];
        }
        coarsecounter+=2*K-1;
        ticcoarse=omp_get_wtime();
    
       
        #pragma omp parallel for
        for (int dir=0; dir<2; dir++){
            for (int i=0; i<K; i++){
                for (int k=0; k<(N+1)*(N+1); k++){
                  
                    UK2[dir][i][k]=UK1[dir][i][k];
                }
            }
        }
        causalsweeping(CG, BC, obs, N, K, wind,F);
        toccoarse = (omp_get_wtime() - ticcoarse);
        sumcoarsetime+=toccoarse;
       
        ticfine=omp_get_wtime();
        #pragma omp parallel for
        for (int c=0; c<M*M; c++){
            UFM[c]=1000;
        }
        
        fineupdate(UFM,CG, F,BC,obs,N,K,wind,fwind,sweepfine);
        tocfine=omp_get_wtime()-ticfine;
        sumfinetime+=tocfine;
        
        
        for (int p=0; p<N*N; p++){
            sumfinesweeps+=sweepfine[p];
        }
        finecounter+=N*N;
        ticfine=omp_get_wtime();
       
        #pragma omp parallel for
        for (int c=0; c<M*M; c++){
            UF2[c]=UF1[c];
            UF1[c]=UFM[c];
        }
        tocfine=(omp_get_wtime()-ticfine);
        sumfinetime+=tocfine;
        double diff=0;
       
        for (int v=0; v<M*M; v++){
            for (int r=0; r< M*M; r++){
                diff=diff + fabs(UFM[v+r*M]-overallfine[v+r*M]);
            }
            
        }
        
        diff=diff/double (M*M);
      
    }
    //k=2
    if (paraiter>=2){
        
        
        ticcoarse=omp_get_wtime();
        #pragma omp parallel for
        for (int dir=0; dir<2; dir++){
            for (int i=0; i<K; i++){
                for (int k=0; k<(N+1)*(N+1); k++){
                    CG[dir][i][k]=1000;
                }
            }
        }
        
        
       
        coarseupdate(CG, UFM, UK1, F,BC,obs,N,K,0,0,wind, fwind,T0, sweepiter,sweepcoarse);
        totalsweepcoarse[0]=sweepcoarse[0];
       // cout << "after k=2 0,0 coarseupdate"<< endl;
        
        
    #pragma omp parallel for //shared(U,wind,UK1)
        for (int p=1;p<2*K-1;p++)
        {
            //cout << p<< endl;
            if (p<K){
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,p,0,wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            else{
                int sweepcoarse[1];
                coarseupdate(CG,UFM, UK1, F,BC,obs,N,K,0,p-(K-1),wind, fwind,T0, sweepiter,sweepcoarse);
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            
        }
        
       
         
        toccoarse = (omp_get_wtime() - ticcoarse);
        sumcoarsetime+=toccoarse;
        
        for (int i=0; i<2*K-1; i++){
            sumcoarsesweeps+=totalsweepcoarse[i];
        }
        coarsecounter+=2*K-1;
        ticcoarse=omp_get_wtime();

       #pragma omp parallel for
        for (int dir=0; dir<2; dir++){
            for (int i=0; i<K; i++){
                for (int k=0; k<(N+1)*(N+1); k++){
                    Ucoarse1[dir][i][k]=CG[dir][i][k];
                }
            }
        }
        
        causalsweeping(CG, BC, obs, N, K, wind,F);
 
        
        toccoarse=(omp_get_wtime()-ticcoarse);
        sumcoarsetime+=toccoarse;
        
        
        ticfine= omp_get_wtime();

        #pragma omp parallel for
        for (int c=0; c<M*M; c++){
            UFM[c]=1000;
        }
        
        fineupdate(UFM,CG, F,BC,obs,N,K,wind,fwind,sweepfine);
        tocfine=omp_get_wtime()-ticfine;
        sumfinetime+=tocfine;
        
        
        for (int p=0; p<N*N; p++){
            sumfinesweeps+=sweepfine[p];
        }
        finecounter+=N*N;
        ticfine=omp_get_wtime();
        #pragma omp parallel for
        for (int c=0; c<M*M; c++){
            UF2[c]=UF1[c];
            UF1[c]=UFM[c];
        }
        tocfine=(omp_get_wtime()-ticfine);
        sumfinetime+=tocfine;
        double diff=0;
        for (int v=0; v<M*M; v++){
            for (int r=0; r< M*M; r++){
                diff=diff + fabs(UFM[v+r*M]-overallfine[v+r*M]);
            }
            
        }
        
        diff=diff/double (M*M);
        
    }
    
    
    
    ////////////////ALGORITHM////////////////
    totaliters=2;
   
    for (int a=3; a<=paraiter; a++){
        totaliters++;
        ticcoarse=omp_get_wtime();
#pragma omp parallel for
        for (int dir=0; dir<2; dir++){
            for (int i=0; i<K; i++){
                for (int k=0; k<(N+1)*(N+1); k++){
                    CG[dir][i][k]=1000;
                }
            }
        }
        
        int w1=2,w2=1,w3=1;
       
        coarseupdatenewtheta(CG, Ucoarse1,UF1,UF2,UK1,UK2,UK3, F,BC,obs,N,K,0,0,wind, fwind, sweepiter, eps, delta, x0,theta,w1,w2,w3,sweepcoarse);
        totalsweepcoarse[0]=sweepcoarse[0];
        
#pragma omp parallel for
        for (int p=1;p<2*K-1;p++)
        {
            if (p<K){
                int sweepcoarse[1];
               coarseupdatenewtheta(CG, Ucoarse1,UF1,UF2, UK1,UK2,UK3, F,BC,obs,N,K,p,0,wind, fwind, sweepiter, eps, delta, x0,theta,w1,w2,w3,sweepcoarse);
                
                
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            else{
                int sweepcoarse[1];
                coarseupdatenewtheta(CG, Ucoarse1,UF1,UF2, UK1,UK2,UK3, F,BC,obs,N,K,0,p-(K-1),wind, fwind, sweepiter, eps, delta, x0,theta,w1,w2,w3,sweepcoarse);
                
               
                totalsweepcoarse[p]=sweepcoarse[0];
            }
            
        }
        toccoarse = (omp_get_wtime() - ticcoarse);
        sumcoarsetime+=toccoarse;
        
        for (int i=0; i<2*K-1; i++){
            sumcoarsesweeps+=totalsweepcoarse[i];
        }
        coarsecounter+=2*K-1;
        ticcoarse=omp_get_wtime();
        
        causalsweeping(CG, BC, obs, N, K, wind,F);
        toccoarse=omp_get_wtime()-ticcoarse;
        sumcoarsetime+=toccoarse;
        
        ticfine=omp_get_wtime();
#pragma omp parallel for
        for (int c=0; c<M*M; c++){
            UFM[c]=1000;
        }
        
        fineupdate(UFM,CG, F,BC,obs,N,K,wind,fwind,sweepfine);
        
        tocfine=omp_get_wtime()-ticfine;
        sumfinetime+=tocfine;
        
        
        for (int p=0; p<N*N; p++){
            sumfinesweeps+=sweepfine[p];
        }
        finecounter+=N*N;
        ticfine=omp_get_wtime();
        
#pragma omp parallel for
        for (int c=0; c<M*M; c++){
            UF2[c]=UF1[c];
            UF1[c]=UFM[c];
        }
        tocfine=(omp_get_wtime()-ticfine);
        sumfinetime+=tocfine;
        
        double diff=0;
        for (int v=0; v<M; v++){
            for (int r=0; r< M; r++){
                diff=diff + fabs(UFM[v+r*M]-overallfine[v+r*M]);
            }
            
        }
        
        diff=diff/double (M*M);
        cout << diff << endl;
        
        if (diff<1e-5 && disp==0){
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
        }
        if (diff<1e-14){
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
        }
        
        
        
    }
    
    
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
    
   
    
    
    for (int dir=0; dir<2; dir++)
        for (int i=0; i<K; i++)
        {
            delete[] CG[dir][i];
            delete[] UK[dir][i];
            delete[] UK1[dir][i];
            delete[] UK2[dir][i];
            delete[] UK3[dir][i];
            delete[] wind[dir][i];
            delete[] Ucoarse1[dir][i];
             delete[] theta[dir][i];
            
        }
    

    
    delete[] marmarray;
    delete[] overallfine;
    //delete[] U;
    delete[] F;
    delete[] BC;
    delete[] obs;
    delete[] sweepfine;
  
   
    delete[] UF2;
    delete[] UF1;
   
    delete[] UFM;

    delete[] fwind;
    //delete[] prev;
    return 0;
}


