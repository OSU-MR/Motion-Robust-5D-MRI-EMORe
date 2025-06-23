function y =  sth(x, lam)
% soft thresholding for wavelet bands
y = zeros(size(x),'like', x); 

for i = 1:size(x,6)
    xabs = abs(x(:,:,:,:,:,i));
    y(:,:,:,:,:,i) = max(xabs - lam(i), 0).* x(:,:,:,:,:,i);
    xabs(xabs == 0) = 1;
    y(:,:,:,:,:,i) = y(:,:,:,:,:,i)./xabs;
end
