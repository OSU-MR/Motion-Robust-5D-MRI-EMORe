clear;
clc;
close all;

x=imread('cardiac02.bmp'); 
x=rgb2gray(x); 
x=im2double(x);
% x=x(151:160, 151:160)
% x=phantom('Modified Shepp-Logan',128);
G=gauss2d(size(x),1);
n=64*1;
sig=0.05;

y=x+sig*randn(size(x));
y2=medfilt2(y,[3,3]);
y3=conv2(y,G,'same');

figure; subplot(221); imagesc(x);  axis('image');
        subplot(222); imagesc(y);  axis('image');
        subplot(223); imagesc(y2); axis('image');
        subplot(224); imagesc(y3); axis('image');
        
Mn = min([x(:); y(:); y2(:); y3(:)]);
Mx = max([x(:); y(:); y2(:); y3(:)]);
        
figure; subplot(221); hist(x(:),n); axis([Mn, Mx, 0,  5*numel(x)/n]);
        subplot(222); hist(y(:),n); axis([Mn, Mx, 0,  5*numel(x)/n]);
        subplot(223); hist(y2(:),n); axis([Mn, Mx, 0, 5*numel(x)/n]);
        subplot(224); hist(y3(:),n); axis([Mn, Mx, 0, 5*numel(x)/n]);
        
%%        
[a2,a1]=hist(x(:),n);
W=gausswin(numel(a1),0.525/sig)'; W=W/sum(W(:));
ac2=1.30*circshift(conv(a2,W,'full'),[0,-n/2+1]);
ac1=linspace(a1(1),2*a1(end), numel(ac2));
[b2,b1]=hist(y(:),n);
figure; plot(a1,a2); hold on; plot(ac1, ac2,'r'); hold on; plot(b1,b2,'g');

