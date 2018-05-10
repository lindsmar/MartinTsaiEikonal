#ifndef consistent_h
#define consistent_h
#include <algorithm>
#include <cmath>
using namespace std;

bool consistent(double prevW, double prevw, int p, int z);

bool consistent(double prevW, double prevw, int p, int z){

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









#endif /* consistent_h */
