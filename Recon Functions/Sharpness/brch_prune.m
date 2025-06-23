function [Io] = brch_prune(I, type)

% Removes braching from a binary image

% Input ================
% I: Input binary image
% type: Methods to remove braches
%       'A' value deletes all pixels (and its 8 neighbors) residing at 
%       the center of braching
%       'B' value deletes minimum pixels and trims shorter of the two
%       braches

% Output ===============
% Io: Image without any braches

% Check the number of arguments
if nargin==1
    type='A'; % Default value
end


Io=I; % Initialize output 
brch_ind= brch_find(Io); % Find pixels that reside at the junctions
if strcmp(type,'A')
    for i=1:numel(brch_ind)
        ind = nbr_ind(Io, brch_ind(i)); % Find 8 neighbors of the pixel
        Io(ind)=0;
    end
elseif strcmp(type, 'B')
    it_max=1000; % Maximum iterations allowed
    it=0;
    while numel(brch_ind) > 0 & it<it_max
        it=it+1;
        ind = nbr_ind(Io, brch_ind(1));
        brch_num=[];
        for j=1:numel(ind)
            I_tmp=Io;
            I_tmp(ind(j))=0;
            brch_num(j)=numel(brch_find(I_tmp));
            rgn_len=regionprops(I_tmp,'Area');

            for k=1:size(rgn_len,1)
                s(k)=rgn_len(k).Area;
            end
            cost(j) = sum(s.^2); % Cost function that needs to be maximized
        end
        tmp=find(brch_num < numel(brch_ind));
        [a,b]=max(cost(tmp));
        J=tmp(b);
        Io(ind(J))=0;
    %     figure; imagesc(Io);
        brch_ind=brch_find(Io);    
    end
else
    error('Only type A and B are valid inputs');
end






