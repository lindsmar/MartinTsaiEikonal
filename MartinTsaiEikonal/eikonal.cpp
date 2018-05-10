//
//  eikonal.cpp
//  Eikonal Project 5
//
//  Created by Lindsay Martin on 9/28/16.
//  Copyright © 2016 Lindsay Martin. All rights reserved.
//

#include <algorithm>
#include <cmath>
#include <matrix.h>
#include <mex.h>
#include "coarseinitial.h"
#include "fineupdate.h"
#include "coarseupdate.h"
#include "consistent.h"
#include "causalsweeping.h"
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //declare variables

    mxArray *U_C, *BDRY, *fVals, *OB, *U_C_new, *U_F_new, *UF, *WIND, *W, *finewind, *uNEW, *UFINEMIN, *Uk, *FW;
    double *U, *BC, *F, *ucoarse,*ufine, *Ufine, *obs,  *wind, *w, *fw, *fwind, *UNEW, *UFM, *UK;
    double h,H,C, T;
    int N, K, M,paraiter,sweepiter;


    //associate inputs

    U_C = mxDuplicateArray(prhs[0]);
    fVals = mxDuplicateArray(prhs[1]);
    BDRY = mxDuplicateArray(prhs[2]);
    OB= mxDuplicateArray(prhs[3]);
    N = (int) mxGetScalar(prhs[4]); // number of nodes on coarse grid
    K= (int) mxGetScalar(prhs[5]); // number of nodes on fine grid
    paraiter= (int) mxGetScalar(prhs[6]);
    sweepiter= (int) mxGetScalar(prhs[7]);
    T = mxGetScalar(prhs[8]);

    H=(double)1/(double) N;
    h=(double)1/(double)(N*K);
    C=0.5;

    //associate pointers for inputs

    U = mxGetPr(U_C);

    F=mxGetPr(fVals);
    BC=mxGetPr(BDRY);
    obs=mxGetPr(OB);

    //associate outputs

    M = N*K+1;

    U_C_new = plhs[0] = mxCreateDoubleMatrix(M,M,mxREAL);
    U_F_new = plhs[1] = mxCreateDoubleMatrix(M,M,mxREAL);
    //UF= plhs[2]=mxCreateDoubleMatrix(M+N-1,M+N-1,mxREAL);
    UFINEMIN=mxCreateDoubleMatrix(M,M,mxREAL);
    W = plhs[2] = mxCreateDoubleMatrix(M,M,mxREAL);
    FW = plhs[3] = mxCreateDoubleMatrix(M,M,mxREAL);


    //associate pointers for outputs

    ucoarse = mxGetPr(U_C_new);
    ufine = mxGetPr(U_F_new);
    w = mxGetPr(W);
    fw = mxGetPr(FW);

    Ufine=mxGetPr(UF);
    UFM=mxGetPr(UFINEMIN);






    //initialize wind direction matrix
    WIND=mxCreateDoubleMatrix(M,M,mxREAL);
    wind=mxGetPr(WIND);
    finewind=mxCreateDoubleMatrix(M,M,mxREAL);
    fwind=mxGetPr(finewind);

    Uk=mxCreateDoubleMatrix(M,M,mxREAL);
    UK=mxGetPr(Uk);



    //// first initialize the coarse grid ////
   int k,c;

        U=coarseinitial(U,F,BC,obs,N,K,0,0,wind,sweepiter);

        int p;
        for (p=1;p<K;p++)
        {
            U=coarseinitial(U,F,BC,obs,N,K,0,p,wind, sweepiter);
            U=coarseinitial(U,F,BC,obs,N,K,p,0,wind, sweepiter);

        }
        for (c=0; c<M*M; c++){
            UK[c]=U[c];
        }

        causalsweeping(U, BC, obs, N, K, wind,F);



        for (c=0; c<M*M; c++){
            UFM[c]=1000;
        }

        fineupdate(UFM,U, F,BC,obs,N,K,wind,fwind);


    ////////////////ALGORITHM////////////////

    for (int a=1; a<=paraiter; a++){


      if (a<=8){



            for (k=0; k<M*M; k++){
                U[k]=1000;

            }

            U=coarseupdate(U, UFM, UK, F,BC,obs,N,K,0,0,wind, fwind,T, sweepiter);

            int p;
            for (p=1;p<K;p++)
            {
                U=coarseupdate(U,UFM, UK, F,BC,obs,N,K,p,0,wind, fwind,T, sweepiter);
                U=coarseupdate(U,UFM, UK, F,BC,obs,N,K,0,p,wind, fwind,T, sweepiter);

            }
            /*for (p=25;p<K;p++)
            {
                U=coarseupdate(U,UFM, UK, F,BC,obs,N,K,p,0,wind, fwind,T, sweepiter);
                U=coarseupdate(U,UFM, UK, F,BC,obs,N,K,0,p,wind, fwind,T, sweepiter);

            }*/
            causalsweeping(U, BC, obs, N, K, wind, F);

            for (c=0; c<M*M; c++){
                UFM[c]=1000;
            }

            fineupdate(UFM,U, F,BC,obs,N,K,wind,fwind);
          }
          else{
            for (k=0; k<M*M; k++){
                U[k]=1000;

            }
            //T=0.001;
            U=coarseupdate(U, UFM, UK, F,BC,obs,N,K,0,0,wind, fwind,T, sweepiter);

            int p;
            for (p=1;p<K;p++)
            {
                U=coarseupdate(U,UFM, UK, F,BC,obs,N,K,p,0,wind, fwind,T, sweepiter);
                U=coarseupdate(U,UFM, UK, F,BC,obs,N,K,0,p,wind, fwind,T, sweepiter);

            }

            causalsweeping(U, BC, obs, N, K, wind, F);

            for (c=0; c<M*M; c++){
                UFM[c]=1000;
            }

            fineupdate(UFM,U, F,BC,obs,N,K,wind,fwind);

          }






    }



    int d;

    for (d=0; d<M*M; d++)
    {
        ufine[d]=UFM[d];//copy final fine and coarse functions to outpu
        ucoarse[d]=U[d];//copy final fine and coarse functions to output
        w[d]=wind[d];
        fw[d]=fwind[d];
    }



    /*for (d=0; d<(M+N-1)*(M+N-1); d++)
     {
     ufine[d]=Ufine[d];
     }*/








    return;
}
