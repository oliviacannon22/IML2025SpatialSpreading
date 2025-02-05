function v = conv_notper(v1,v2,dx)
N = length(v1);
v1Pad=[v1;zeros(size(v1))];
v2Pad=[v2;zeros(size(v2))];
convPad = ifft(dx*fft(circshift(v2Pad,N/2)).*fft(v1Pad));
v=convPad(N+1:end);
end