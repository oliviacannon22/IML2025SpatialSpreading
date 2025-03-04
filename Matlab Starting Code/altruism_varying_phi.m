% Evolve a 2D square domain (want to make rectangular)
% Boundary conditions set by the Laplacian operator (2nd derivative) in
% both space and altruism parameter 
%any time you see Y here, it is altruism parameter, and X is space

close all; clear;

video_on = 0; % change to 1 to turn on video
if video_on
    v = VideoWriter('Video_name','Uncompressed AVI');
end

% System parameters
% add system parameters here
par.d =1; %death rate
par.sd_a=0.5; %width of altruism convolution kernel
par.c=1; %cost of altruism to individual
par.g0=5; %general birth param (probably could be scaled out)
par.sd_rc=2; %width of competition kernel (should be wider than altruism kernel)
par.kD=1e-2; %diffusion coeff. of motility
par.kDy=.01*par.kD; %right now, modeling mutation of altruism by diffusion, much slower than motility diffusion
par.K=20; %carrying capacity, roughly speaking
par.b0=0.5; %param for saturating nonlin
par.b_max=2; %param for saturating nonlin
par.mu=1e-3; %mutation probability
par.m=5e-3; %not really used for us (but might later)
par.sd_m=sqrt(2*par.kD/par.d); %"scale of motility" (different way of expressing diffusion coeff)
%par.phi=0.08; 


% Set numerical and system parameters
dt = 0.005;
numPar.tf = 1000;  % Final time
t = 0:dt:numPar.tf;
iter = length(t);
n_plot=10; %how often we plot, or every ___ time steps

disp(['Iter: ' num2str(iter)]);

% Numerical parameters (length of domain and number of grid points)
%currently, nx and ny have to be the same, but want to change that
%eventually
numPar.Ly = 1;
numPar.Lx = 15;
numPar.nx = 200;
numPar.ny = 150;
numPar.dx = numPar.Lx/(numPar.nx-1);
numPar.dy = numPar.Ly/(numPar.ny-1);

numPar.xgrid = 'FD_Periodic'; % FD_Periodic = finite differences periodic, F_Periodic = Fourier, Periodic BC (assumes 2pi periodic)
numPar.ygrid = 'FD'; %FD = finite differences Neumann 
numPar.order = '2'; % Order of numerical scheme
par.Ly = numPar.Ly;
par.Lx = numPar.Lx;

x = 0:numPar.dx:numPar.Lx;  % Domain
y = 0:numPar.dy:numPar.Ly;


% Define initial condition 
%U = 20*ones(numPar.nx,numPar.ny);
%U(numPar.nx/2:end,numPar.ny/2:end) = 200;
%U(numPar.nx/2:end,1:numPar.ny/2) = 50;
%U = U(:);
 U = 1*ones(numPar.nx,numPar.ny);
 U(end-numPar.nx/10:end,round(5*numPar.ny/8):round(7*numPar.ny/8)) = 5;
 %U(numPar.nx/2:end,1:numPar.ny/2) = 0;
 %U(numPar.nx/2-numPar.nx/10:numPar.nx/2,round(numPar.ny/8):round(3*numPar.ny/8)) = 5;
% U(1:numPar.nx/2,numPar.ny/2:end) = 0;
U = U + normrnd(0,.2,numPar.nx,numPar.ny);

%now make it a vector 
U = U(:);


%Compute matrices that approximate 2nd derivative for x and y (for later, note U must be
%grouped by that variable)
[L2x,L2y]=Laplacians(numPar);
L2x = par.kD.*L2x; %multiply componentwise by diffusion coefficient
L2y = par.kDy.*L2y;

%Code for task 3
%Matrix for implicit one-time step
DUx = speye(numPar.nx*numPar.ny) - dt*L2x;

% Matrix for implicit
D2Ux = speye(numPar.nx*numPar.ny) - dt/2*L2x;  %  block matrix - this implementation assumes the same number of x and y gridpoints. 
D2Uy = speye(numPar.nx*numPar.ny) - dt/2*L2y;  %  block matrix - this implementation assumes the same number of x and y gridpoints.
% Prepare for first step
tmp_d2y = groupX(L2y*groupY(U,numPar),numPar);              % this is d^2(U)/dy^2 (explicit term) (Need to group by Y to take Y derivative, then regroup by X) . 

if video_on
    open(v);
end

   %plot initial condition 
   %this step makes a temporary U that is a proper square (vs vector) and plots it, 
   tmpU = reshape(U,numPar.nx,numPar.ny)';
                figure(1); pcolor(x,y,tmpU); shading interp; 
                colorbar;
                drawnow;


for k = 1:iter
        %this is the main step where all the work gets done

        % Grouped by x at beginning of each loop
       % Evaluate nonlinear term!
       fU = altruismnonlin_varyingphi(U,par,numPar); 

        
        
        %code for upgrade to one-time step

        %U = groupX(U,numPar);
       

        %After changing mutation term, we use one-time step

        % U = DUx \( U + dt * mu*( (fU + dU)*K - (fU + dU)) );
        
        U = DUx \( U + dt * fU); 
        %equation this is solving:
        %U_new = U_old + dt*(d^2/dx^2(U_new) + mu*( (fU + dU)*K - (fU + dU) ) )

        % Get ready for next step
        %tmp_d2y = groupX(L2y*U,numPar);              % this is d^2(U)/dy^2 again
        %U = groupX(U,numPar);   % Switch to X grouping 
       

        %plot current solution 
         if mod(k,n_plot) == 1
            
                %tmpU = reshape(U,numPar.nx,numPar.ny);
                tmpU = reshape(U,numPar.nx,numPar.ny)';
                figure(1); pcolor(x,y,tmpU); shading interp;
                title(['time=' num2str(k*dt), ' total pop = ' num2str(numPar.dx^2*sum(U)) ]);
                colorbar;
                drawnow;

                if video_on
                    fr = getframe(gcf);
                    writeVideo(v,fr);
                end

         end

end

 

if video_on

    close(v);
end


