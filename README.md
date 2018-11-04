# MartinTsaiEikonal

  *  Code for reproduction of numerical results produced in the paper “A Multiscale Domain Decomposition
     Algorithm For Boundary Value Problems For Eikonal Equations”
  *  The numerical algorithm developed by Lindsay Martin and Richard Tsai.
 
  Instructions:
		
  1. Open “examples.m"
  2. Select parameters:
  
    a. “N” N=1/H where H is the coarse grid step size.
				
    b. “K” K=1/(Nh) where h is the overall fine grid step size.
				
    c. "test" accepts an integer value between 1 to 8 and is used to load examples
				
    d. "paraiter" is the number of parareal iterations 
				
    e. "sweepiter" is the number of sweeping iterations for the coarse grids
				
    f. "gamma", "delta", and "x0" are parameters from stability analysis section of the paper
				
    g. "T0" is an initial theta value
    
    (a good example starting point for the above values is T0=.01, gamma=.01, delta=.001, x0=.01)
  3. Run the script.

  Note: if necessary, precompile the *.cpp files to generate the matlab mex files, by uncommenting lines at the top of “examples.m”
  
  To run the parallel tests:
  
  Prerequistes: OpenMP
  
  1. Change Makefile so that the correct directory to compiler and correct OpenMP flag is given.
  2. Compile with make.
  3. Run executable with inputs N,K,number of threads, T0,gamma,delta,x0 (as defined above).
  
  To perform scale test, update make file to compile scaletest.cpp instead of main.cpp.
  1. Compile with make.
  2. Run executable with inputs N,K, T0,gamma,delta,x0 (as defined above).
