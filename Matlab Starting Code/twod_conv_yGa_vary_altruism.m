function w = twod_conv_yGa_vary_altruism(v1,Nx,Gx,dx)
%convolve a function v1 with a function G which can be written as a product Gx*Gy
%v1 is a vector representing a function of two variables, which is stacked
%at present by columns (grouped by X)
%Nx is number of x grid points

%calculate number of y grid points (in case it's different)
Ny=length(v1)/Nx; 

yvec=linspace(0,1,Ny);
yvec=repmat(yvec,1,Nx)';
dy=1/Ny;

%initialize the eventual output vector
w=zeros(length(v1),1);
wx=zeros(Nx,1);

%define the 1D convolution for x (only different if not a square)
convx =@(v1,v2) ifft(dx*fft(circshift(v2,Nx/2)).*fft(v1));

%multiply pointwise by y
v1=v1.*yvec;

%for each fixed x, add up all the y*v1(x,y) over y
for i = 1:Nx
    start = Ny*(i-1) + 1;
    wx(i)=sum(dy*v1(start:start+Ny-1)); %result is a vector of length Nx
end 


%now convolve with Gx in the x direction
wx=convx(wx,Gx);
w=repelem(wx,Ny);

return 
