#ifndef causalsweeping_H
#define causalsweeping_H
void causalsweeping(double *Ucoarse, double *BC,double *obs, int N,int K,double *wind,double* R);
void causalsweeping(double *Ucoarse, double *BC,double *obs, int N,int K,double *wind,double *R){
  int M;
  M=N*K+1;
  double h;
  h=(double)1/(double) (N*K);
  int s1,s2,i,j,i1;
  for (int sweep=0; sweep<4; sweep++){
   for (s2 = -1; s2 <= 1; s2 += 2){
     for (s1 = -1; s1 <= 1; s1 += 2){
    //s1=-1;
    //s2=1;
      for (i = (s1 < 0 ? M-1 : 0); (s1 < 0 ? i>=0 : i<= M-1);  i+=K*s1){
        for (j = (s2 <0 ? M-1 : 0); (s2 < 0 ? j>=0 : j<= M-1); j+=s2){
          if (BC[i+j*M]<0){
            if (wind[i+j*M]==1 || wind[i+j*M]==5 || wind[i+j*M]==2){// || wind[i+j*M]==6){
              if ((j!=0)){
                if (Ucoarse[i+j*M]<Ucoarse[i+(j-1)*M]){
                  Ucoarse[i+j*M]=Ucoarse[i+(j-1)*M];//+R[i+j*M]*h;
                }
              }
            }
            if (wind[i+j*M]==4 || wind[i+j*M]==7 || wind[i+j*M]==3){
              if (j!=M-1){
                if (Ucoarse[i+j*M]<Ucoarse[i+(j+1)*M]){
                  Ucoarse[i+j*M]=Ucoarse[i+(j+1)*M];//+R[i+j*M]*h;
                }
              }
            }
            if (j%K==0){
              if (i!=M-1){
                if (wind[i+j*M]==4 || wind[i+j*M]==8 || wind[i+j*M]==1){
                  if (Ucoarse[i+j*M]<Ucoarse[i+1+j*M]){
                    Ucoarse[i+j*M]=Ucoarse[i+1+j*M];//+R[i+j*M]*h;
                  }
                }
              }
              if (i!=0){
                if (wind[i+j*M]==3 || wind[i+j*M]==6 || wind[i+j*M]==2){
                  if (Ucoarse[i+j*M]<Ucoarse[i-1+j*M]){
                    Ucoarse[i+j*M]=Ucoarse[i-1+j*M];//+R[i+j*M]*h;
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
                if (wind[i1+j*M]==4 || wind[i1+j*M]==8 || wind[i1+j*M]==1){// || wind[i1+j*M]==5){ //only added ==5 for maze problem
                  if (i1!=M-1){
                    if (Ucoarse[i1+j*M]<Ucoarse[i1+1+j*M]){
                      Ucoarse[i1+j*M]=Ucoarse[i1+1+j*M];//+R[i+j*M]*h;
                    }
                  }
                }
                if (wind[i1+j*M]==3 || wind[i1+j*M]==6 || wind[i1+j*M]==2 ){
                  if (i1!=0){
                    if (Ucoarse[i1+j*M]<Ucoarse[i1-1+j*M]){
                      Ucoarse[i1+j*M]=Ucoarse[i1-1+j*M];//+R[i+j*M]*h;
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


}

#endif
