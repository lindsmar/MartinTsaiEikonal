#ifndef smooth_H
#define smooth_H

double smooth(double theta, double amplifier, double delta, double x0);
double smooth(double theta, double amplifier, double delta, double x0){

  double sigma;
  sigma=1/(1+exp((theta-x0)/amplifier));

  double thetaused;
  thetaused=sigma*theta+delta*(1-sigma)*theta;

  return thetaused;

}

#endif //smooth_H
