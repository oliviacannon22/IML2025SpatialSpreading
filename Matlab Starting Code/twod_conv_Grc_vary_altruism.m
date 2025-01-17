function w = twod_conv_Grc_vary_altruism(v1,Nx,Gx,dx)
%convolve a function v1 with a function G which can be written as a product Gx*Gy
%v1 is a vector representing a function of two variables, which is stacked
%at present by columns (grouped by X)
%Nx is number of x grid points

%calculate number of y grid points (in case it's different)
Ny=length(v1)/Nx; 

%initialize the eventual output vector
w=zeros(length(v1),1);

%define the 1D convolution for x (only different if not a square)
convx =@(v1,v2) ifft(dx*fft(circshift(v2,Nx/2)).*fft(v1));

%for each fixed x, convolve in the y direction
% for i = 1:Nx
%     start = Ny*(i-1) + 1;
%     v1(start:start+Ny-1)=convy(v1(start:start+Ny-1),Gy);
% end 
%now for each y, convolve with the resulting function in the x
%direction
for j = 1:Ny
    w(j:Nx:end)=convx(v1(j:Nx:end),Gx);
end 

return 
