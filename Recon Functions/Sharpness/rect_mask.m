function [rct] = rect_mask(n,r)

rct=zeros(n);
r=floor(r);
rct(r+1:end-r, r+1:end-r)=1;
rct=logical(rct);