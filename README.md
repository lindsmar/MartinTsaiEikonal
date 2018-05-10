# MartinTsaiEikonal

  *  Code for reproduction of numerical results produced in the paper “A Multiscale Domain Decomposition
     Algorithm For Boundary Value Problems For Eikonal Equations”
  *  The numerical algorithm developed by Lindsay Martin and Richard Tsai.
 
  Instructions:
  1. Open “examples.m”
  2. Select parameters:
    a. “N” N=1/H where H is the coarse grid step size.
    b. “K” K=1/(Nh) where h is the overall fine grid step size.
    c. "test" accepts an integer value between 1 to 8 and is used to load examples
      d. gamma, delta, and x0 are parameters from stability analysis section of the paper
    e. T0 is an initial theta value
  3. Run the script.

  Note: if necessary, precompile the *.cpp files to generate the matlab mex files, by uncommenting lines at the top of “examples.m”
