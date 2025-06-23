function [ y ] = fft2_shift( x )
%FFT3_SHIFT take a unitary fft on the first three dimensions of an array with
%fftshifts

y = ifftshift(ifftshift(x,1),2);
y = fft(fft(y,[],1),[],2);
y = 1/sqrt(size(x,1)*size(x,2))*fftshift(fftshift(y,1),2);


end

