function [ y ] = fftshift3d( x )

    y= fftshift(fftshift(fftshift(x,1),2),3);
end

