function [y] = TV(x, Dims)
    nDims = numel(Dims);
    y = zeros([Dims,nDims],'like', x);  
    for i = 1:nDims
        idx = repmat({':'}, 1, nDims);  % Create colons for each dimension
        idx{nDims+1} = i;  % Add the i-th index as the last dimension
        y(idx{:}) = x-circshift(x,-1,i);
        % Set the boundary differences to zero along the current dimension
        idx{i} = size(x, i);
        y(idx{:}) = 0;
    end
end