function [I] = len_thresh(I,L)

% Throws away regions (edges) that are smaller that length L

I=bwlabel(I);
for i=1:max(I(:))
    len=sum(sum(I==i));
    if len < L
        I(I==i)=0;
    else 
%         lab(lab==i)=1;
    end
end
I=logical(I);