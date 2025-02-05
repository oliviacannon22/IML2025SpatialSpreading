function [L2x,L2y] = Laplacians(numPar)
% Linear operators on a rectangular grid with specified boundary conditions
% Fourier periodic assumes on 2pi periodic frame.
%% rename parameters
	nx = numPar.nx;
	ny = numPar.ny;
    order = numPar.order;
	Lx = numPar.Lx;
	Ly = numPar.Ly;
    

%% X direction
switch numPar.xgrid
    case 'F_Periodic' % Fourier with periodic boundary conditions. Assumes 2pi periodic frame
        [~,D2x] = fourdif(nx,2);    % 2nd derivative matrix        

        
            
    case 'FD' % Finite Differences with Neumann Boundary conditions
        switch order 
            case '2'
            hx = Lx/(nx-1);  
            ex = ones(nx,1);

            D2x = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx) - sparse(1:nx,1:nx,ex,nx,nx); % 2nd derivative matrix
            D2x = (D2x + D2x'); 
            D2x(1,1:2) = [-2,2];
            D2x(nx,nx-1:nx)=[2, -2];     % Neumann boundary conditions     
            D2x = D2x/hx^2;

            case '4' % 4th order 
            hx = Lx/(nx-1);

            D2x = sparse(1:nx-1,[2:nx-1 nx],16*ones(nx-1,1),nx,nx) - sparse(1:nx-2,[3:nx-1 nx],ones(nx-2,1),nx,nx);
            D2x = (D2x + D2x' - 30*speye(nx))/12; 
            D2x(1:2,:) = 0; D2x(end-1:end,:) = 0; D2x(2,1:3) = [1, -2, 1];  % Neumann boundary conditions: use 2nd order at boundary
            D2x(end-1,end-2:end) = [1, -2, 1]; D2x(1,1:2) = [-2,2]; D2x(end,end-1:end) = [2,-2];
            D2x = D2x/hx^2; % 2nd derivative matrix  
        end


    case 'FD_Periodic' % Finite Differences with Periodic Boundary conditions
       
        switch order 
            case '2'
            hx = Lx/nx;  
            ex = ones(nx,1);

            D2x = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx) - sparse(1:nx,1:nx,ex,nx,nx); % 2nd derivative matrix
            D2x = (D2x + D2x'); 
            D2x(end,1) = 1; D2x(1,end) = 1; % Periodic boundary conditions    
            D2x = D2x/(hx^2);

            case '4' % 4th order 
            hx = Lx/nx;
            D2x = sparse(1:nx-1,[2:nx-1 nx],16*ones(nx-1,1),nx,nx) - sparse(1:nx-2,[3:nx-1 nx],ones(nx-2,1),nx,nx);
            D2x = (D2x + D2x' - 30*speye(nx))/12; 
            D2x(1,end-1:end) = [-1, 16]; D2x(2,end) = -1; D2x(end-1,1) = -1; D2x(end,1:2) = [16,-1]; % Periodic boundary conditions
            D2x = D2x/(12*hx^2); % 2nd derivative matrix  
        end
        
end
Ix  = speye(nx,nx);

% Y-direction
switch numPar.ygrid
    case 'F_Periodic'   % Fourier with periodic boundary conditions. Assumes 2pi periodic frame
        [~,D2y] = fourdif(nx,2);    % 2nd derivative matrix        


    case 'FD'   % Finite Differences with Neumann boundary conditions
        switch order 
            case '2'
            hy = Ly/(ny-1);  
            ey = ones(ny,1);

            D2y = sparse(1:ny-1,[2:ny-1 ny],ones(ny-1,1),ny,ny) - sparse(1:ny,1:ny,ey,ny,ny); % 2nd derivative matrix
            D2y = (D2y + D2y'); 
            D2y(1,1:2) = [-2,2];
            D2y(ny,ny-1:ny)=[2, -2];     % Neumann boundary conditions     
            D2y = D2y/hy^2;

            case '4' % 4th order 
 
            hy = Ly/(ny-1);
            D2y = sparse(1:ny-1,[2:ny-1 ny],16*ones(ny-1,1),ny,ny) - sparse(1:ny-2,[3:ny-1 ny],ones(ny-2,1),ny,ny);
            D2y = (D2y + D2y' - 30*speye(ny))/12; 
            D2y(1:2,:) = 0; D2y(end-1:end,:) = 0; D2y(2,1:3) = [1, -2, 1];  % Neumann boundary conditions: use 2nd order at boundary
            D2y(end-1,end-2:end) = [1, -2, 1]; D2y(1,1:2) = [-2,2]; D2y(end,end-1:end) = [2,-2];
            D2y = D2y/hy^2; % 2nd derivative matrix  
        end

end
Iy = speye(ny,ny);


% Second Derivative matrices 
L2x = kron(Iy,D2x);
L2y = kron(Ix,D2y);


L2x = sparse(L2x);
L2y = sparse(L2y);




