% UNCOMMENT NEXT 4 LINES TO COMPILE MEX FILES %
% mex eikonal.cpp %-new method with fixed theta
% mex eikonalnewtheta.cpp %-new method with theta_{i,j}^{k,used}
% mex eikonalnvariabletheta.cpp %-new method with M_{i,j}^k
% mex eikonal_homog2.cpp %-new method with theta_{i,j}^{k,used} where
        % coarse solver solves the homogenized equation
% mex fastsweeporig.cpp %FSM for overall fine solution


%%%EXAMPLES$%%%
%N=1/H
%K=1/(Nh)
%test cases are the following:
%1: R=1+.99sin(2pix)sin(2piy)
%2: R=1+.5sin(20pix)sin(20piy)
%3: distance function to (0.5,0.5)
%4: maze with curved obstacles
%5: maze with fast obstacle
%6: generalization of test case 2
%7: homogenization example
%8: random slowness function 

%gamma,delta,x0 are parameters from stability analysis section
%T0 intial theta value or theta value for eikonal.cpp


N=10;
K=15;
test=7;
%compute overall fine solution%
[ORIG,bdrycond,R]=testcase(N,K,test);



%number of parareal iterations and sweeping iterations
paraiter=1; 
sweepiter=5;
gamma=0.01;
delta=0.01;
x0=.1;
T0=0.1;

%initialize solution
Uold=1000*ones(N*K+1,N*K+1);                                                                                                                                                                    
obs=zeros(N*K+1,N*K+1);

% homogenized speed function for test case 7
%   cbar=[1,sqrt(2)/2,1,sqrt(2)/2,1,sqrt(2)/2,1,sqrt(2)/2];

[Ucoarse,Ufine, wind,fwind]=eikonalnewtheta(Uold,R,bdrycond,obs,N,K,paraiter,sweepiter,T0,gamma,delta,x0);
%   [Ucoarse,Ufine, wind,fwind,thetavals]=eikonalvariabletheta(Uold,R,bdrycond,obs,N,K,paraiter,sweepiter,ORIG);
%   [Ucoarse,Ufine, wind,fwind]=eikonal(Uold,R,bdrycond,obs,N,K,paraiter,sweepiter,T0);
%   [Ucoarse_homog,Ufine_homog, wind_homog,fwind_homog]=eikonal_homog2(Uold,R,bdrycond,obs,N,K,paraiter,sweepiter,T,0.01,0.001,.02,cbar);

diff=abs(ORIG-Ufine);

Linferror=max(diff(:))

L1error=sum(diff(:))/(N*K+1)^2 

%contour of fine solutions patched together
figure
X=linspace(0,1,N*K+1);
Y=flip(X);
contourf(X,Y, Ufine,linspace(0,2.5,150))