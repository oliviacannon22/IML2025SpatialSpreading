% System parameters
% add system parameters here
par.d = 1; %death rate
par.sd_a=0.1; %width of altruism convolution kernel
par.c= 0.5; %cost of altruism to individual
par.g0=5; %general birth param (probably could be scaled out)
par.sd_rc=2; %width of competition kernel (should be wider than altruism kernel)
par.kD=1e-3; %diffusion coeff. of motility
par.kDy=.01*par.kD; %right now, modeling mutation of altruism by diffusion, much slower than motility diffusion
par.K=20; %carrying capacity, roughly speaking
par.b0=0.5; %param for saturating nonlin
par.b_max=2; %param for saturating nonlin
par.mu=1e-3; %mutation probability
par.m=1e-2; %
par.sd_m=sqrt(2*par.kD/par.d); %"scale of motility" (different way of expressing diffusion coeff)
%par.phi=0.08; 

numPar.tf = 500;  % Final time

% Numerical parameters (length of domain and number of grid points)
%currently, nx and ny have to be the same, but want to change that
%eventually
numPar.Ly = 1;
numPar.Lx = 30;
numPar.nx = 400;
numPar.ny = 150;
numPar.dx = numPar.Lx/(numPar.nx-1);
numPar.dy = numPar.Ly/(numPar.ny-1);

numPar.xgrid = 'FD_Periodic'; % FD_Periodic = finite differences periodic, F_Periodic = Fourier, Periodic BC (assumes 2pi periodic)
numPar.ygrid = 'FD'; %FD = finite differences Neumann 
numPar.order = '2'; % Order of numerical scheme
par.Ly = numPar.Ly;
par.Lx = numPar.Lx;

altruism_varying_phi(par,numPar,1)