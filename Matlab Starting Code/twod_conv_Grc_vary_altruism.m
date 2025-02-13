function w = twod_conv_Grc_vary_altruism(v1,Nx,Gx,dx,dy)

%v1 is a vector representing a function of two variables, which is stacked
%at present by columns (grouped by X)
%Nx is number of x grid points

%calculate number of y grid points (in case it's different)
Ny=length(v1)/Nx; 

%initialize the eventual output vector 
wx=zeros(Nx,1); 
w=zeros(length(v1),1);

%define the 1D convolution for x (only different if not a square)
%uses fft or 'fast fourier transform', since convolution acts in fourier space by 
%multiplying fourier coefficients of the same frequency, so it's a fast way to compute it
convx =@(v1,v2) ifft(dx*fft(circshift(v2,Nx/2)).*fft(v1));

%for each fixed x, add up in the y direction
for i = 1:Nx
     start = Ny*(i-1) + 1;
     wx(i)=sum(v1(start:start+Ny-1))*dy;
     %v1(start:start+Ny-1)=sum(v1(start:start+Ny-1))*dy*ones(Ny,1);
 end 
%now for each y, convolve with the resulting function in the x
%direction
%for j = 1:Ny
%    w(j:Nx:end)=convx(v1(j:Nx:end),Gx);
%end 

wx=convx(wx,Gx);
w=repelem(wx,Ny);

return 
