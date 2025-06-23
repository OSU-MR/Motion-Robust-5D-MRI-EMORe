function [ y ] = fft3(x)
%fft3_SHIFT Summary of this function goes here
%   Detailed explanation goes here
y = sqrt(size(x,1)*size(x,2)*size(x,3))*fft(fft(fft(x,[],1),[],2),[],3);
end

