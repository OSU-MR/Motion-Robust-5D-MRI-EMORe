function[val] = gauss2d(n,sig)

[x,y] = ndgrid(1:n(1), 1:n(2));
cntr=floor(n/2)+1;
val = 1/(2*pi*sig^2)*exp(-((x-cntr(1)).^2 + (y-cntr(2)).^2)./(2*sig^2));
val=val/sum(val(:));