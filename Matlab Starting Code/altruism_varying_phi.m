% Evolve a 2D square domain (want to make rectangular)
% Boundary conditions set by the Laplacian operator (2nd derivative) in
% both space and altruism parameter 
%any time you see Y here, it is altruism parameter, and X is space

close all; clear;

video_on = 0;
if video_on
    v = VideoWriter('varyingaltruism','MPEG-4');
end

% System parameters
% add system parameters here
par.d =1; %death rate
par.sd_a=1; %width of altruism convolution kernel
par.c=1; %cost of altruism to individual
par.g0=5; %general birth param (probably could be scaled out)
par.sd_rc=4; %width of competition kernel (should be wider than altruism kernel)
par.kD=3e-2; %diffusion coeff. of motility
par.kDy=.05*par.kD; %right now, modeling mutation of altruism by diffusion, much slower than motility diffusion
par.K=10; %carrying capacity, roughly speaking
par.b0=0.5; %param for saturating nonlin
par.b_max=2; %param for saturating nonlin
par.mu=1e-3; %mutation probability
par.m=5e-3; %not really used for us (but might later)
par.sd_m=sqrt(2*par.kD/par.d); %"scale of motility" (different way of expressing diffusion coeff)
%par.phi=0.08; 


% Set numerical and system parameters
dt = 0.01;
numPar.tf = 1000;  % Final time
t = 0:dt:numPar.tf;
iter = length(t);
n_plot=10; %how often we plot, or every ___ time steps

disp(['Iter: ' num2str(iter)]);

% Numerical parameters (length of domain and number of grid points)
%currently, nx and ny have to be the same, but want to change that
%eventually
numPar.Ly = 1;
numPar.Lx = 10;
numPar.nx = 150;
numPar.ny = 150;
numPar.dx = numPar.Lx/(numPar.nx-1);
numPar.dy = numPar.Ly/(numPar.ny-1);

numPar.xgrid = 'F_Periodic'; % FD = finite differences, F_Periodic = Fourier, Periodic BC
numPar.ygrid = 'FD'; %FD = finite differences Neumann 
numPar.order = '2'; % Order of numerical scheme
par.Ly = numPar.Ly;
par.Lx = numPar.Lx;

x = 0:numPar.dx:numPar.Lx;  % Domain
y = 0:numPar.dy:numPar.Ly;


% Define initial condition 
%U = 100*ones(numPar.nx,numPar.ny);
%U(numPar.nx/2:end,numPar.ny/2:end) = 200;
%U(numPar.nx/2:end,1:numPar.ny/2) = 50;
%U = U(:);
U = 1*ones(numPar.nx,numPar.ny);
U(end-numPar.nx/10:end,numPar.ny/2:end) = 75;
U(numPar.nx/2:end,1:numPar.ny/2) = 0;
U(numPar.nx/2-numPar.nx/10:numPar.nx/2,1:numPar.ny/2) = 50;
U(1:numPar.nx/2,numPar.ny/2:end) = 0;
U = U + normrnd(0,.2,numPar.nx,numPar.ny);
U = U(:);

%"Compute Linear Operator" outputs the matrix for the finite difference
%approximation of (here) the Laplacian (2nd derivative) 
%boundary conditions given by what we put in for numPar.xgrid (and ygrid)
%ComputeLinearOperator would output four outputs but we only want the
%3rd,so that's why the squiggles are there

% Compute Linear Operator for evolution
[~,~,L2x,~] = ComputeLinearOperator_rectangular(numPar); 
L2x = par.kD.*L2x; %multiply componentwise by diffusion coefficient 

%the following part is sort of hacky and is to get a Laplacian
%with different boundary conditions and coefficient for when we diffuse in
%altruism. but it's ok if it makes no sense
numPar.xgrid='FD';numPar.Lx=numPar.Ly;%change the xgrid bc the 2nd deriv always groups by x, but need one Neumann for y
%FD = finite differences Neumann 
[~,~,L2y,~] = ComputeLinearOperator_rectangular(numPar); 
L2y = par.kDy.*L2y;
numPar.Lx=10;

% Matrix for implicit
D2Ux = speye(numPar.nx*numPar.ny) - dt/2*L2x;  %  block matrix - this implementation assumes the same number of x and y gridpoints. 
D2Uy = speye(numPar.nx*numPar.ny) - dt/2*L2y;  %  block matrix - this implementation assumes the same number of x and y gridpoints.
% Prepare for first step
tmp_y1 = groupX(L2y*groupY(U,numPar),numPar);              % this is d^2(U)/dy^2 (explicit term) (Need to group by Y to take Y derivative, then regroup by X) . 

if video_on
    open(v);
end

   %plot initial condition 
   %this step makes a temporary U that is a proper square (vs vector) and plots it, 
   tmpU = reshape(U,numPar.nx,numPar.ny);
                figure(1); pcolor(x,y,tmpU); shading interp; 
                colorbar;
                drawnow;


for k = 1:iter
        %this is the main step where all the work gets done
        %we use an 'implicit method' to solve equations with 2nd
        %derivatives. Since we have two of these (d^2/dx^2 and d^2/dy^2), we split dt in two and do
        %two 'half steps': one for x, one for y

        % Grouped by x at beginning of each loop
        % First "half step": Implicit in X, explicit in Y   - grouped by x 
       % Evaluate nonlinear term!
       fU = altruismnonlin_varyingphi(U,par,numPar); 
        
       U = D2Ux \ ( U+ dt/2 * tmp_y1 + dt/2 * fU);   % Grouped by X 
       %equation this is solving:
       %U_new = U_old + dt/2*(d^2/dx^2(U_new) + d^2/dy^2(U_old) + f(U_old)

        % Second "half step": Implicit in Y, explicit in X - grouped by y
        tmp_x1 = groupY(L2x*U,numPar);      % this is d^2(U)/dx^2 (explicit term) - (calculate while grouped by X, then regroup in Y)
        
       % Evaluate nonlinear terms at the most recent values 
        fU = altruismnonlin_varyingphi(U,par,numPar); % put into function before switching grouping, bc function takes U grouped by X
        
        U = groupY(U,numPar);   % Switch to Y grouping 

        U = D2Uy \ (U + dt/2 * tmp_x1 + dt/2 * groupY(fU,numPar) );  % Grouped by Y
        %equation this is solving:
        %U_new = U_old + dt/2*(d^2/dx^2(U_old) + d^2/dy^2(U_new) + f(U_old)
       
        % Get ready for next step
        tmp_y1 = groupX(L2y*U,numPar);              % this is d^2(U)/dy^2 again
        U = groupX(U,numPar);   % Switch to X grouping 
       
        %plot current solution 
         if mod(k,n_plot) == 1
            
                tmpU = reshape(U,numPar.nx,numPar.ny);
                figure(1); pcolor(x,y,tmpU); shading interp;
                title(['time=' num2str(k*dt), ' total pop = ' num2str(numPar.dx^2*sum(U)) ]);
                colorbar;
                drawnow;

                if video_on
                    fr = getframe();
                    writeVideo(v,fr);
                end

         end

end

 

if video_on

    close(v);
end


