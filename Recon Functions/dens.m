function Y = dens(X, a, b, mu)
    Y = abs(X - mu) / a;
    Y = exp(-Y.^b);
    Y = Y * b / (2 * a);
    Y = Y / exp(gammaln(1 / b));
end