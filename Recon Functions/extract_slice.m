function y = extract_slice(x,slice)
    nDims = ndims(x);
    idx = repmat({':'}, 1, nDims);  % Create colons for each dimension
    idx{3} = slice;
    y = squeeze(x(idx{:}));
end