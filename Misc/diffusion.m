close all; clear;

video_on = 0;
if video_on
    v = VideoWriter('standard_diffusion','MPEG-4');
end

c = 5; %diffusion coefficient 

% Set numerical and system parameters
dt = 0.005;
numPar.tf = 100;  % Final time
t = 0:dt:numPar.tf;
iter = length(t);
n_plot=10; %how often we plot, or every ___ time steps

% Numerical parameters (length of domain and number of grid points)
%currently, nx and ny have to be the same, but want to change that
%eventually
numPar.Lx = 40;
numPar.nx = 10001;
numPar.dx = numPar.Lx/(numPar.nx-1);

x = 0:numPar.dx:numPar.Lx;

%Initial Conditions
U = zeros(numPar.nx,1);
U(5001) = 1000;

%Second Derrivative Operator with Neumann BC
D = sparse(1:numPar.nx,1:numPar.nx,[-1,-2*ones(1,numPar.nx-2),-1],numPar.nx,numPar.nx);
E = sparse(2:numPar.nx,1:numPar.nx-1,ones(1,numPar.nx-1),numPar.nx,numPar.nx);
S = E+D+E';
S = (c/numPar.dx^2).*S; %Second Derrivative Operator
A = speye(numPar.nx)-(dt).*S;

if video_on
    open(v);
end

   %plot initial condition 
   figure(1); plot(x,U); 
   drawnow;

for k = 1:iter 
       U = A \ U; %Implicit Euler

       if k < 20 %plot readjusts for first few timesteps
           M = max(U);
       end

        %plot current solution 
         if mod(k,n_plot) == 1
            
                figure(1); plot(x,U); ylim([-0.1 M]);
                title(['time = ' num2str(k*dt), ' total pop = ' num2str(numPar.dx^2*sum(U))]);
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


