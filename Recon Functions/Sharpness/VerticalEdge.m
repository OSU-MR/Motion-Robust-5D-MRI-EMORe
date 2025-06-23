clear;
% Creates a 2D phantom used to compare fitting based shapreness measure
% against the gradient based results

% For fitting use Im(k).Data(i).Param(:,3) while for gradient based method
% use Im(k).Data(i).Param instead in sharp_main

% For fitting use use P=[P;p] and option 'three' in sigmoid_fit.m
% For gradient use tmpp=sort(abs(diff(dt))); P=[P;sqrt(mean(tmpp.^2))]; and
% use option 'one' in sigmoid_fit.m
m=19;
n=256+floor(m-1);
I=ones(n,n);
I(:,1:floor(n/2))=0;
I=conv2(I,gauss2d([m,m],0.75),'valid'); 
VerticalEdge=I;
figure; imagesc(I); colormap('gray'); axis('image');
% size(I),
save VerticalEdge2 VerticalEdge;

