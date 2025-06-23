clc;
% close all;
clear;


p.e_sig   = 1; % Width of smoothing filter for edge detection
p.e_th    = 0.03; % Thresholding for edge detection
p.e_tl    = p.e_th*0.5;
p.c_L     = 7; % Length used to compute curvature for normals
p.e_L     = p.c_L+1; % Discard all edges below length e_L
p.n_L     = 7; % Length (in pixels) of lines normal to the selected edges
p.n_s     = 6*p.n_L+1; % Number of samples along each normal line; keep it odd
p.hist    = 'no';
p.h_dc    = 0.02; % Histogram parameters (dc)
p.h_ln    = 1024; % Number of bins in the histogram
p.h_sig   = p.h_ln/3; % std of histogram
p.i_fn    = 1; % Frame number to be considered for analysis
p.i_ref   = [1];
% p.AngR    = [-90 -75; 75 90]*pi/180;
p.AngR    = [-90 90]*pi/180; % E

E=imread('Morph01.bmp');
E=rgb2gray(E);
E=im2double(E);
E=logical(E);
figure; imagesc(E); colormap(gray); axis('image');


% E=edge(I,'canny',[0.01 0.02],0.1);
% figure; imagesc(E); colormap(gray); axis('image');


E_b = brch_prune(E,'A');
figure; imagesc(E_b); colormap(gray); axis('image');

E_b = len_thresh(E_b,20);
figure; imagesc(E_b); colormap(gray); axis('image');

E_bs = len_thresh(E_b,50);
figure; imagesc(E_bs); colormap(gray); axis('image');

E_bsm   = morph_mask(E).*E_bs;
figure; imagesc(E_bsm); colormap(gray); axis('image');





