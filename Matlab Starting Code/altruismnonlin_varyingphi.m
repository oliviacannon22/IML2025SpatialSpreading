function f = altruismnonlin_varyingphi(p,par,numPar) 
%takes a vector p representing a 2D function, grouped by X
%N is number of grid points in each of x and y (so assumes it is a square)
%sd_a is altruism standard dev and sd_rc is competition standard dev
%everything else is parameters from paper
dx=numPar.dx;
N=numPar.nx;
g0=par.g0;
c=par.c;
%phi=par.phi;
b_max=par.b_max;
b0=par.b0;
d=par.d;
K=par.K;
sd_a=par.sd_a;
sd_rc=par.sd_rc;
%cap=par.cap;



M=N*dx;

yvec=linspace(0,1,N);
yvec=repmat(yvec,1,N)';

%define convolution kernels
xx=linspace(-M/2,M/2,N);
gaus = @(x,sig,amp)amp*exp(-(((x).^2)/(2*sig.^2)));
%make gaussian for altruism
G_a=zeros(N,1);
for i = 1:N
    G_a(i)=gaus(xx(i),sd_a,1/(sqrt(2*pi)*sd_a)); 
end 
%make gaussian for competition
G_rc=zeros(N,1);
for i = 1:N
    G_rc(i)=gaus(xx(i),sd_rc,1/(sqrt(2*pi)*sd_rc)); 
end 


%ignore the commented out nonlinearities below -- from earlier simplified versions
f = g0*p.*(1-c*yvec + (b_max*twod_conv_yGa_vary_altruism(p,N,G_a,dx))./(b_max/b0+twod_conv_yGa_vary_altruism(p,N,G_a,dx))).*(1-1/K*twod_conv_Grc_vary_altruism(p,N,G_rc,dx))... %full model varying altruism
...%g0*p.*(1-c*phi + (b_max*phi*twodimconv(p,N,G_a,G_a,dx))./(b_max/b0 + phi*twodimconv(p,N,G_a,G_a,dx))).*(1-twodimconv(p,N,G_rc,G_rc,dx)/K)... %keeps altruism constant
...%g0*p.*(1-c*phi + (b0*phi*twodimconv(p,N,G_a,G_a,dx))-twodimconv(p,N,G_rc,G_rc,dx)/K-cap*p)... %quadratic only, constant altruism
-d*p; 


return 


