function [ y ] = ifft3(x)
%IFFT3_SHIFT Summary of this function goes here
%   Detailed explanation goes here
y = sqrt(size(x,1)*size(x,2)*size(x,3))*ifft(ifft(ifft(x,[],1),[],2),[],3);
end

