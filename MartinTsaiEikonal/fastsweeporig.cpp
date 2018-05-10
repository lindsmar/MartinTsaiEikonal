
#include <algorithm>
#include <cmath>
#include <matrix.h>
#include <mex.h>


using namespace std;


/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    //declare variables
    mxArray *uVals,*BDRY, *fVals, *UNEW, *OB, *breakiter, *WIND, *W;
    const mwSize *dims;
    double *u, *f, *Unew, *BC, *obs, *BreakIter, *wind, *w;
    int N, numdims;
    double h, old;
    int s, m, n,iter;
    double TOL=0.00000001;





    //associate inputs
    uVals = mxDuplicateArray(prhs[0]);
    fVals = mxDuplicateArray(prhs[1]);
    BDRY = mxDuplicateArray(prhs[2]);
    OB = mxDuplicateArray(prhs[3]);
    iter=(int) mxGetScalar(prhs[4]);
    h = mxGetScalar(prhs[5]);

    //figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    N = (int)dims[0];

    //associate outputs
    WIND=mxCreateDoubleMatrix(N,N,mxREAL);
    wind=mxGetPr(WIND);

    UNEW = plhs[0] = mxCreateDoubleMatrix(N,N,mxREAL);
    breakiter = plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
     W = plhs[2] = mxCreateDoubleMatrix(N,N,mxREAL);


    //associate pointers
    u = mxGetPr(uVals);
    f = mxGetPr(fVals);
    BC = mxGetPr(BDRY);
    obs = mxGetPr(OB);
   w = mxGetPr(W);
    Unew = mxGetPr(UNEW);
    BreakIter= mxGetPr(breakiter);

    //algorithm
    //h=(double) 1/(double) (N-1);
    double exititer=0;
    double maxError=1;
    int tempwindx, tempwindy;
    for (s=0;s<iter;s++){
        exititer++;
        int s1, s2, i,j;
        double unew;
        double m1,m2;
        double discr;
        maxError=0;
        double currError;


        for (s1=1; s1>=-1; s1-=2){
           for (s2=1; s2>=-1; s2-=2){
             //for (s2 = -1; s2 <= 1; s2 += 2){
               //for (s1 = -1; s1 <= 1; s1 += 2){
        //s1=1;
        //s2=1;

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
                                tempwindy=7;
                            }
                            else if (j==N-1){
                                m1=u[i+(j-1)*N];
                                tempwindy=5;
                            }
                            else{
                                m1=min(u[i+(j-1)*N],u[i+(j+1)*N]);
                                if (u[i+(j-1)*N]<u[i+(j+1)*N]){
                                    tempwindy=5;
                                }
                                else{
                                    tempwindy=7;
                                }
                            }
                            if (i==0){
                                m2=u[(i+1)+j*N];
                                tempwindx=8;
                            }
                            else if (i==N-1){
                                m2=u[(i-1)+j*N];
                                tempwindx=6;
                            }
                            else{
                                m2=min(u[(i-1)+j*N], u[(i+1)+j*N]);
                                if (u[(i-1)+j*N]<u[(i+1)+j*N]){
                                    tempwindx=6;
                                }
                                else{
                                    tempwindx=8;
                                }
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
                              if (tempwindx==6 && tempwindy==5){
                                  wind[i+j*N]=2;
                              }
                              else if (tempwindx==6 && tempwindy==7){
                                  wind[i+j*N]=3;
                              }
                              else if (tempwindx==8 && tempwindy==5){
                                  wind[i+j*N]=1;
                              }
                              else if (tempwindx==8 && tempwindy==7){
                                  wind[i+j*N]=4;
                              }
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
                              if (m1<m2){
                                  wind[i+j*N]=tempwindy;
                              }
                              else{
                                  wind[i+j*N]=tempwindx;
                              }
                                u[i+j*N]=unew;
                            }


                        }



                            currError=u[j+i*N]-unew;
                            currError=abs(currError);
                            maxError=max(maxError,currError);
                            //u[i+j*N]=min(u[i+j*N],unew);
                        }
                    }
                }
           }
        }

        /*if (maxError <TOL){

            break;
        }*/

    }

    for (m=0; m<N*N; m++)
       {
            Unew[m]=u[m];
             w[m]=wind[m];
        }
    BreakIter[0]=exititer;

    return;
}
