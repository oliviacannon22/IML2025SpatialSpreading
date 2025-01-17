function [L1x,L1y,L2x,L2y] = ComputeLinearOperator_rectangular(numPar)
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
        [~,D1x] = fourdif(nx,1);    % 1st derivative matrix
        
            
    case 'FD' % Finite Differences with Neumann Boundary conditions
        switch order 
            case '2'
            hx = Lx/(nx-1);  
            D1x = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx); 
            D1x = D1x - D1x';
            D1x(1,:) = 0; D1x(end,:) = 0;   % Neumann BC
            D1x = D1x/(2*hx);               % 2nd order FD centered difference first derivative matrix
            
            ex = ones(nx,1);

            D2x = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx) - sparse(1:nx,1:nx,ex,nx,nx); % 2nd derivative matrix
            D2x = (D2x + D2x'); 
            D2x(1,1:2) = [-2,2];
            D2x(nx,nx-1:nx)=[2, -2];     % Neumann boundary conditions     
            D2x = D2x/hx^2;

            case '4' % 4th order 
            hx = Lx/(nx-1);
            D1x = sparse(1:nx-1,[2:nx-1 nx],8*ones(nx-1,1),nx,nx) - sparse(1:nx-2,[3:nx-1 nx],ones(nx-2,1),nx,nx);
            D1x = (D1x - D1x')/12; D1x(1:2,:) = 0; D1x(end-1:end,:) = 0; D1x(2,1:4) = [-2, -3, 6, -1]/6; D1x(end-1,end-3:end) = [1, -6, 3, 2]/6; 
            D1x = D1x/hx; % First derivative matrix

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
            D1x = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx); 
            D1x = D1x - D1x';
            D1x(1,end) = -1; D1x(end,1) = 1; % Periodic boundary conditions
            D1x = D1x/(2*hx);               % 2nd order FD centered difference first derivative matrix
            
            ex = ones(nx,1);

            D2x = sparse(1:nx-1,[2:nx-1 nx],ones(nx-1,1),nx,nx) - sparse(1:nx,1:nx,ex,nx,nx); % 2nd derivative matrix
            D2x = (D2x + D2x'); 
            D2x(end,1) = 1; D2x(1,end) = 1; % Periodic boundary conditions    
            D2x = D2x/(hx^2);

            case '4' % 4th order 
            hx = Lx/nx;
            D1x = sparse(1:nx-1,[2:nx-1 nx],8*ones(nx-1,1),nx,nx) - sparse(1:nx-2,[3:nx-1 nx],ones(nx-2,1),nx,nx);
            D1x = (D1x - D1x'); 
            D1x(1,end-1:end) = [1, -8]; D1x(2,end) = 1; D1x(end-1,1) = -1; D1x(end,1:2) = [8,-1]; % Periodic boundary conditions
            D1x = D1x/(12*hx); % First derivative matrix

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
        [~,D1y] = fourdif(nx,1);    % 1st derivative matrix

    case 'FD'   % Finite Differences with Neumann boundary conditions
        switch order 
            case '2'
            hy = Ly/(ny-1);  
            D1y = sparse(1:ny-1,[2:ny-1 ny],ones(ny-1,1),ny,ny); 
            D1y = D1y - D1y';
            D1y(1,:) = 0; D1y(end,:) = 0;   % Neumann BC
            D1y = D1y/(2*hy);               % 2nd order FD centered difference first derivative matrix

            ey = ones(ny,1);

            D2y = sparse(1:ny-1,[2:ny-1 ny],ones(ny-1,1),ny,ny) - sparse(1:ny,1:ny,ey,ny,ny); % 2nd derivative matrix
            D2y = (D2y + D2y'); 
            D2y(1,1:2) = [-2,2];
            D2y(ny,ny-1:ny)=[2, -2];     % Neumann boundary conditions     
            D2y = D2y/hy^2;

            case '4' % 4th order 
            hy = Ly/(ny-1);
            D1y = sparse(1:ny-1,[2:ny-1 ny],8*ones(ny-1,1),ny,ny) - sparse(1:ny-2,[3:ny-1 ny],ones(ny-2,1),ny,ny);
            D1y = (D1y - D1y')/12; D1y(1:2,:) = 0; D1y(end-1:end,:) = 0; D1y(2,1:4) = [-2, -3, 6, -1]/6; D1y(end-1,end-3:end) = [1, -6, 3, 2]/6; 
            D1y = D1y/hy; % First derivative matrix

            D2y = sparse(1:ny-1,[2:ny-1 ny],16*ones(ny-1,1),ny,ny) - sparse(1:ny-2,[3:ny-1 ny],ones(ny-2,1),ny,ny);
            D2y = (D2y + D2y' - 30*speye(ny))/12; 
            D2y(1:2,:) = 0; D2y(end-1:end,:) = 0; D2y(2,1:3) = [1, -2, 1];  % Neumann boundary conditions: use 2nd order at boundary
            D2y(end-1,end-2:end) = [1, -2, 1]; D2y(1,1:2) = [-2,2]; D2y(end,end-1:end) = [2,-2];
            D2y = D2y/hy^2; % 2nd derivative matrix  
        end

end
Iy = speye(ny,ny);

% First derivative matrices
L1x = kron(Iy,D1x);  
L1y = kron(D1y,Ix);

% Second Derivative matrices 
L2x = kron(Iy,D2x);
L2y = kron(D2y,Ix);

L1x = sparse(L1x);
L1y = sparse(L1y);
L2x = sparse(L2x);
L2y = sparse(L2y);




