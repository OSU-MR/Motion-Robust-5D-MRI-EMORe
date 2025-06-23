function k = scale_k(ref,x)
    k = (ref'*x)./norm(x).^2;
end