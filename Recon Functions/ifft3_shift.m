function [ y ] = ifft3_shift(x)
%IFFT3_SHIFT Summary of this function goes here
%   Detailed explanation goes here
y = ifftshift(ifftshift(ifftshift(x,1),2),3);
y = ifft(ifft(ifft(y,[],1),[],2),[],3);
y = sqrt(size(x,1)*size(x,2)*size(x,3))*fftshift(fftshift(fftshift(y,1),2),3);
end

