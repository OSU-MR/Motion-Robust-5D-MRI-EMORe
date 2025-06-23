function [y]=adjTV(x,Dims)
    nDims=numel(Dims);
    y=zeros(Dims,'like', x); 
    idx = repmat({':'}, 1, nDims);  % Create colons for all dimensions
    for i=1:nDims
        idx{nDims+1} = i;  % Add the i-th index as the last dimension
        y=y+x(idx{:})-circshift(x(idx{:}),1,i);
    end
end