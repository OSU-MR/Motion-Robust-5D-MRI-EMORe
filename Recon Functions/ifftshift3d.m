function [ y ] = ifftshift3d( x )

    y= ifftshift(ifftshift(ifftshift(x,1),2),3);
end

