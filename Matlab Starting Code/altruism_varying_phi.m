% Evolve a 2D square domain (want to make rectangular)
% Boundary conditions set by the Laplacian operator (2nd derivative) in
% both space and altruism parameter 
%any time you see Y here, it is altruism parameter, and X is space
function f = altruism_varying_phi(par,numPar,video_on) 

close all;

Video_name='Some_stripes';
if video_on
    v = VideoWriter(Video_name,'Uncompressed AVI');
   v.FrameRate=20;
end

% Set numerical and system parameters
dt = 0.05;
t = 0:dt:numPar.tf;
iter = length(t);
n_plot=40; %how often we plot, or every ___ time steps

%disp(['Iter: ' num2str(iter)]);

x = 0:numPar.dx:numPar.Lx;  % Domain
y = 0:numPar.dy:numPar.Ly;

% Define initial condition 
%U = 20*ones(numPar.nx,numPar.ny);
%U(numPar.nx/2:end,numPar.ny/2:end) = 200;
%U(numPar.nx/2:end,1:numPar.ny/2) = 50;
%U = U(:);
% U = 1*ones(numPar.nx,numPar.ny);
U=zeros(numPar.nx,numPar.ny);

 %stripes
 %alt_stripes start with higher altruism, superalt_stripes even higher. 
 s_wid=round(numPar.nx/numPar.Lx); %width of one stripe
 n_stripes = 2; %total stripes
 n_superalt_stripes=2; %how many very altruistic stripes
 n_selfish_stripes=0; %how many selfish stripes
 n_alt_stripes=n_stripes-n_selfish_stripes; %how many middle-altruistic stripes
 s_starts=linspace(1,numPar.nx,n_stripes+1);
 s_starts=s_starts(1:end-1);
 
 %currently, layout left to right is superalt, then alt, then selfish, but could play
 %with having every other, for instance
 if n_superalt_stripes > 0
    for i = 1:n_superalt_stripes %altruism between 3/8 and 7/8
    U(s_starts(i):s_starts(i)+s_wid,round(1*numPar.ny/8):round(5*numPar.ny/8)) = 5;
    end
 end 

 for i = n_superalt_stripes+1:n_alt_stripes %altruism between 1/4 and 1/2
    U(s_starts(i):s_starts(i)+s_wid,round(2*numPar.ny/8):round(4*numPar.ny/8)) = 5;
 end
%selfish stripes

for i = n_alt_stripes+1:n_stripes %selfish stripes, altruism between 1/8 and 3/8
 U(s_starts(i):s_starts(i)+s_wid,round(1*numPar.ny/8):round(3*numPar.ny/8)) = 5;
end
 %U(numPar.nx/2:end,1:numPar.ny/2) = 0;
 %U(numPar.nx/2-numPar.nx/10:numPar.nx/2,round(numPar.ny/8):round(3*numPar.ny/8)) = 5;
% U(1:numPar.nx/2,numPar.ny/2:end) = 0;
%U = U + normrnd(0,.02,numPar.nx,numPar.ny);

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

for k = 1:iter
        %this is the main step where all the work gets done

        % Grouped by x at beginning of each loop
       % Evaluate nonlinear term!
       fU = altruismnonlin_varyingphi(U,par,numPar); 


        %code for task 3

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
                
             %this first section plots, in the top half of a figure, the total
             %population at each x value (to visualize size of different
             %colonies). It also calculates the average altruism level at
             %each x

              
                
                tmpUx=zeros(1,numPar.nx);
                avgAltx=zeros(1,numPar.nx);
                Ugy=groupY(U,numPar);
                yvec=linspace(0,1,numPar.ny)';
                for i = 1:numPar.nx
                start = numPar.ny*(i-1) + 1;
                tmpUx(i)=sum(Ugy(start:start+numPar.ny-1)*numPar.dy);
                avgAltx(i)=sum(yvec.*Ugy(start:start+numPar.ny-1)*numPar.dy);
                %v1(start:start+Ny-1)=sum(v1(start:start+Ny-1))*dy*ones(Ny,1);
                end 
                avgAltx=avgAltx./tmpUx;
                tiledlayout(3,1)
                nexttile
                plot(x,tmpUx)
                title(['time=' num2str(k*dt), 'Population at each x value'])
                ylim([0 300]) 
                nexttile([2 1])
       
            
                %plot the current solution at time t
                tmpU = reshape(U,numPar.nx,numPar.ny)';
                ymax=round(7*numPar.ny/8); %plot only part of picture to see more clearly
                tmpU=tmpU(1:ymax,:);
                pcolor(x,y(1:ymax),tmpU); shading interp;
                clim([0 10000]);
                hold on
                plot(x,avgAltx.*(tmpUx>0.7*max(tmpUx)),'o','Color','red') %plot the average altruism level on top of the figure, but only for where a lot of population is 
                title(['time=' num2str(k*dt), ' total pop = ' num2str(numPar.dx*numPar.dy*sum(U)) ]);
                
                colorbar;
                %drawnow;
                hold off

                if video_on
                    fr = getframe(figure(1));
                    writeVideo(v,fr);
                end
                
                %this plot is to see where new agents are being added.  
                % the log (ln) of the net new agents added (it is not plotting where agents are dying). 
                % the log is there to be able to distinguish between zero and very small
                %numbers better than a color plot usually could. Matlab has a machine precision of 10^-16, so note that a
                %log between -16 to -20ish can just mean 0. 
                % figure(3)
                % tmpfU = reshape(fU,numPar.nx,numPar.ny)';
                % tmpfU=tmpfU(1:ymax,:);
                % tmpfU=tmpfU.*(abs(tmpfU)>1e-14); 
                % tmpfUpos=log(tmpfU.*(tmpfU>0));
                % pcolor(x,y(1:ymax),tmpfUpos); shading interp;
                % title(['time=' num2str(k*dt), 'log of new agents added' ]);
                % colorbar;
                % drawnow;

                

         end

end

 

if video_on

    close(v);
end

return


