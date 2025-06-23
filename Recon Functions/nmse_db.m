function [y,k] = nmse_db(ref,x)
    k = (ref'*x)./norm(x).^2;
    y = 20.*log10(norm(ref-k.*x)./norm(ref));

end
