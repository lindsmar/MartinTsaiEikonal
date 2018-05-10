//
//  eikonal_homog.cpp
//  Eikonal Project 5
//
//  Created by Lindsay Martin on 9/28/16.
//  Copyright © 2016 Lindsay Martin. All rights reserved.
//

#include <algorithm>
#include <cmath>
#include <matrix.h>
#include <mex.h>
#include "coarseinitial_ani.h"
#include "fineupdate.h"
#include "coarseupdate_ani.h"
#include "coarseupdate_ani2.h"
#include "consistent.h"
#include "causalsweeping.h"
using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //declare variables

    mxArray *U_C, *BDRY, *fVals, *OB, *U_C_new, *U_F_new, *UF, *WIND, *W, *finewind, *uNEW, *UFINEMIN, *Uk, *FW,*CBAR,*Uk1,  *Uk2, *uf1, *uf2, *ucoarse1, *Uk3;
    double *U, *BC, *F, *ucoarse,*ufine, *Ufine, *obs,  *wind, *w, *fw, *fwind, *UNEW, *UFM, *UK, *cbar, *UK1, *UK2, *UF1, *UF2, *Ucoarse1, *UK3;
    double h,H,C, T0, eps, delta, x0;
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
    T0=mxGetScalar(prhs[8]);
    eps=mxGetScalar(prhs[9]);
    delta=mxGetScalar(prhs[10]);
    x0=mxGetScalar(prhs[11]);

    CBAR=mxDuplicateArray(prhs[12]);
    H=(double)1/(double) N;
    h=(double)1/(double)(N*K);
    C=0.5;

    //associate pointers for inputs

    U = mxGetPr(U_C);

    F=mxGetPr(fVals);
    BC=mxGetPr(BDRY);
    obs=mxGetPr(OB);
    cbar=mxGetPr(CBAR);
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
    Uk=mxCreateDoubleMatrix(M,M,mxREAL);
    UK=mxGetPr(Uk);
    Uk1=mxCreateDoubleMatrix(M,M,mxREAL);
    UK1=mxGetPr(Uk1);
    Uk2=mxCreateDoubleMatrix(M,M,mxREAL);
    UK2=mxGetPr(Uk2);
    Uk3=mxCreateDoubleMatrix(M,M,mxREAL);
    UK3=mxGetPr(Uk3);
    uf2=mxCreateDoubleMatrix(M,M,mxREAL);
    UF2=mxGetPr(uf2);
    uf1=mxCreateDoubleMatrix(M,M,mxREAL);
    UF1=mxGetPr(uf1);
    ucoarse1=mxCreateDoubleMatrix(M,M,mxREAL);
    Ucoarse1=mxGetPr(ucoarse1);


    //// first initialize the coarse grid ////
   int k,c;

        U=coarseinitial_ani(U,cbar,BC,obs,N,K,0,0,wind,sweepiter);

        int p;
        for (p=1;p<K;p++)
        {
            U=coarseinitial_ani(U,cbar,BC,obs,N,K,0,p,wind, sweepiter);
            U=coarseinitial_ani(U,cbar,BC,obs,N,K,p,0,wind, sweepiter);

        }


        for (c=0; c<M*M; c++){
            UK1[c]=U[c];
            UK3[c]=U[c];
        }
    causalsweeping(U, BC, obs, N, K, wind,F);


    for (c=0; c<M*M; c++){
        UFM[c]=1000;
    }

    fineupdate(UFM,U, F,BC,obs,N,K,wind,fwind);

    for (c=0; c<M*M; c++){
        UF1[c]=UFM[c];
    }

    if (paraiter>=1){
        for (k=0; k<M*M; k++){
            U[k]=1000;

        }
        U=coarseupdate_ani(U, UFM, UK1, cbar,BC,obs,N,K,0,0,wind, fwind,T0, sweepiter);

        int p;
        for (p=1;p<K;p++)
        {
            U=coarseupdate_ani(U,UFM, UK1, cbar,BC,obs,N,K,p,0,wind, fwind,T0, sweepiter);
            U=coarseupdate_ani(U,UFM, UK1, cbar,BC,obs,N,K,0,p,wind, fwind,T0, sweepiter);

        }
        for (c=0; c<M*M; c++){
            Ucoarse1[c]=U[c];
            UK2[c]=UK1[c];
        }
        causalsweeping(U, BC, obs, N, K, wind,F);
        for (c=0; c<M*M; c++){
            UFM[c]=1000;
        }

        fineupdate(UFM,U, F,BC,obs,N,K,wind,fwind);

        for (c=0; c<M*M; c++){
            UF2[c]=UF1[c];
            UF1[c]=UFM[c];
        }
      }

      if (paraiter>=2){
          for (k=0; k<M*M; k++){
              U[k]=1000;

          }

          U=coarseupdate_ani(U, UFM, UK1, cbar,BC,obs,N,K,0,0,wind, fwind,T0, sweepiter);

          int p;
          for (p=1;p<K;p++)
          {
              U=coarseupdate_ani(U,UFM, UK1, cbar,BC,obs,N,K,p,0,wind, fwind,T0, sweepiter);
              U=coarseupdate_ani(U,UFM, UK1, cbar,BC,obs,N,K,0,p,wind, fwind,T0, sweepiter);

          }
          for (c=0; c<M*M; c++){
              Ucoarse1[c]=U[c];
          }

          causalsweeping(U, BC, obs, N, K, wind,F);

          for (c=0; c<M*M; c++){
              UFM[c]=1000;
          }

          fineupdate(UFM,U, F,BC,obs,N,K,wind,fwind);

          for (c=0; c<M*M; c++){
              UF2[c]=UF1[c];
              UF1[c]=UFM[c];
          }

      }


    ////////////////ALGORITHM////////////////

    for (int a=3; a<=paraiter; a++){

            for (k=0; k<M*M; k++){
                U[k]=1000;

            }



            U=coarseupdate_ani2(U, Ucoarse1,UF1,UF2,UK1,UK2,UK3, cbar,BC,obs,N,K,0,0,wind, fwind, sweepiter, eps, delta, x0,2,6,6);


            for (p=1;p<K;p++)
            {
                U=coarseupdate_ani2(U, Ucoarse1,UF1,UF2, UK1,UK2,UK3, cbar,BC,obs,N,K,p,0,wind, fwind, sweepiter, eps, delta, x0,2,6,6);
                U=coarseupdate_ani2(U, Ucoarse1,UF1,UF2, UK1,UK2,UK3, cbar,BC,obs,N,K,0,p,wind, fwind, sweepiter, eps, delta, x0,2,6,6);

            }
            causalsweeping(U, BC, obs, N, K, wind,F);

            for (c=0; c<M*M; c++){
                UFM[c]=1000;
            }

            fineupdate(UFM,U, F,BC,obs,N,K,wind,fwind);

            for (c=0; c<M*M; c++){
                UF2[c]=UF1[c];
                UF1[c]=UFM[c];
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
