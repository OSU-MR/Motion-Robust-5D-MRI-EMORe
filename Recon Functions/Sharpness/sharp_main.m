%% Copy Rights and Credits
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
%
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.
%-------------------------------------------------------------------------- 
% The code is written by Rizwan Ahmad (Rizwan.Ahmad@osumc.edu).  


% clear;
% clc;
% close all;


%% User defined parameters
p.e_sig   = 1.7; % Width of smoothing filter for edge detection              % Imp
p.e_th    = 0.05; % Thresholding for edge detection
p.e_tl    = p.e_th*0.4;
p.c_L     = 3; % Length used to compute curvature for normals (in option B only)
p.e_L     = p.c_L+1; % Discard all edges below length e_L
p.n_L     = 7; % Length (in pixels) of lines normal to the selected edges % Imp
p.n_s     = 8*p.n_L+1; % Number of samples along each normal line; keep it odd
p.hist    = 'no'; % Histogram tweaking for the image;
p.h_dc    = 0.02; % Histogram parameters (dc)
p.h_ln    = 1024; % Number of bins in the histogram
p.h_sig   = p.h_ln/3; % std of histogram
p.i_fn    = 1; % Frame number to be considered for analysis
p.i_ref   = [1]; % The image index used to detect edges; use a number or 'all' for averaged
p.AngR    = [-90, 90]*pi/180; % 
% p.AngR    = [-90, -72; 72, 90]*pi/180; % 
p.sig     = 1e-3*0; % Noise added to the image
p.n_f     = 1; % 1/std of profile blurring filter; smaller number=more blurring
p.trim    = 'one'; % Method of profile truncation, 'one', 'two', 'three', 'four'
p.typ     = 'cf';    % Sharpness measurement method, '1090', 'grad', 'cf'
p.alph    = 0.01; % Alpha for t-testing


%% Add path
% addpath('.\common'); % Helper functions
% addpath('.\NFFT');


%% Load Image
%name=uipickfiles('filterspec', 'C:\Users\Rizwan\Research\SHARPNESS');%\Sharpness_Perfusion\Sharpness');
%name=uipickfiles('filterspec', 'C:\MATLAB7\work\Work_Space_Dual_CPU\2013\2013_06_13_Try_Boundary_Sharpness');%\Sharpness_Perfusion\Sharpness');
% name = 'C:\MATLAB7\work\Work_Space_Dual_CPU\2013\2013_KW_Filter\sharpness_temp_img.mat'
% for i=1:size(name,2)
%fid = fopen(char(name)); 
%Raw=load(char(name));
%Raw = struct2cell(Raw);
%if numel(Raw)>1
%    p.i_fn=Raw{2}(2);% ??? if empty
%end
%Raw=Raw{1};
%fclose(fid);
%if numel(size(Raw))==4
%    I=squeeze(Raw(:,:,p.i_fn,:)); % space, space, time, method
%elseif numel(size(Raw))==2 || numel(size(Raw))==3
%    I=Raw;
%else
%    error('Uploaded MAT file should be an array of size 2D to 4D');
%end
% end

% I=sqrt((p.sig*randn(size(I))).^2+ (p.sig*randn(size(I))).^2) + I;
clear Raw;
fclose('all');

% I=permute(I,[2,1,3]); % Swap phase-encoding and frequency-encoding
% directions
%% 
% cs = 0;
% if(cs == 1)
%     I = double(gather(extract_slice(abs(permute(flip(flip(x_cs(:,:,:,:,1),3),2),[2,3,1])),round(size(x,1)./2))));
% else
%     I = double(gather(extract_slice(abs(permute(flip(flip(x(:,:,:,:,1),3),2),[2,3,1])),round(size(x,1)./2))));
% end
% I=10*I/max(I(:));


cs = 0;
if(cs == 1)
    I = double(gather(abs(squeeze(flip(flip(x_cs(round(size(x,1)./2),:,:,:,1),3),2)))));
else
    I = double(gather(abs(squeeze(flip(flip(x(round(size(x,1)./2),:,:,:,1),3),2)))));
end
I=10*I/max(I(:));
%% Histogram Tweaking
hgram=hshape(p);
if strcmp(p.i_ref,'all')
    Ic=mean(I,3);
    if strcmp(p.hist,'yes'), Ic=histeq(Ic/max(Ic(:)),hgram); end
else
    Ic=I(:,:,p.i_ref);
    if strcmp(p.hist,'yes'), Ic=histeq(Ic/max(Ic(:)),hgram); end
end
% for i=1:size(I,3)
%     I(:,:,i)=histeq(I(:,:,i)/max(max(I(:,:,i))),hgram);
% end

%% Edge Detection
E = edge(Ic,'canny',[p.e_tl p.e_th],p.e_sig);


%% Remove branches
E_b = brch_prune(E,'A');
E_bl = len_thresh(E_b,p.e_L);


%% Plotting
% figure(1);
% for i=1:size(I,3);
%     subplot(round(sqrt(size(I,3))),ceil(sqrt(size(I,3))),i); 
%     imagesc(I(:,:,i)); colormap(gray); axis('image','off');
%     title(['Method ', num2str(i)]);
% end

figure(2);  subplot(221); imagesc(Ic); colormap(gray); axis('image','off'); title('composite image');
            subplot(222); imagesc(E_bl); colormap(gray); axis('image','off'); title('edges');
figure(3);  imagesc(E_bl); colormap(gray); axis('image','off'); title('using mouse, select appropriate edges');
E_bs=bwselect(E_bl,8);
%E_bs=bwselect(E_bl, x_boundary, y_boundary, 8);
figure(2);  subplot(223); imagesc(E_bs); colormap(gray); axis('image','off'); title('user selected edges');
%             subplot(224); imagesc(slct_e); colormap(gray); axis('square','off');



%% Morphological filtering
E_bsm   = morph_mask(E).*E_bs;
% E_bsm = len_thresh(E_bsm,p);


%% Remove pixels close to the edge
E_bsm = rect_mask(size(Ic),ceil(p.n_L/2+1)).*E_bsm;
E_bsm = len_thresh(E_bsm,p.e_L);


%% Sorting Edge Indices
St = edge_sort(E_bsm,p);


%% Find normals to local curvature
St = edge_norm(St,Ic,p,'A');
E_bsms=edge_final(St, E_bsm); 
figure(2);  subplot(224); imagesc(E_bsms); colormap(gray); axis('image','off'); title('edges for sharpness measurement');



%% Curve fitting
av=zeros(size(I,3),1);
mx=zeros(size(I,3),1);
mn=zeros(size(I,3),1);
H=zeros(size(I,3), size(I,3)); % hypothesis results
pval=zeros(size(I,3), size(I,3)); % p values
for k=1:size(I,3)
    %tic
    % Read Image along the lines
    Im(k).Data = read_int(St,I(:,:,k)); % ???? index


    % Curve fitting
    Im(k).Data = sigmoid_fit(Im(k).Data,p);
    Im(k).Sh=[];
    Im(k).IntS=[];
    Im(k).Fit=[];
    for i=1:size(St,1)
        if ~isempty(Im(k).Data(i).Param)
%             Im(k).Sh=[Im(k).Sh;Im(k).Data(i).Param];
            Im(k).Sh=[Im(k).Sh;Im(k).Data(i).Param(:,1)]; %??
            Im(k).IntS=[Im(k).IntS;Im(k).Data(i).IntS];
            Im(k).Fit=[Im(k).Fit;Im(k).Data(i).Fit];
        end
    end 
%     if k>1
        for g=1:1:k
%             [H(g,k),pval(g,k)]=ttest2(Im(g).Sh, Im(k).Sh,p.alph,'both','unequal'); % t-test
            [pval(g,k),H(g,k)]=ranksum(Im(g).Sh, Im(k).Sh, 'alpha',p.alph); % Wilcoxon ranksum

        end
%     end
    md(k)=median(Im(k).Sh)
    av(k)=mean(Im(k).Sh)
    sd(k)=std(Im(k).Sh)
    mn(k)=min(Im(k).Sh)
    mx(k)=max(Im(k).Sh)
    %{k, toc}
end
H=(H+H')./(eye(size(H,1))+ones(size(H)));
pval=(pval+pval')./(eye(size(H,1))+ones(size(H)));






% Ding 2013-07-24, start
% %% Plotting
% load cmap;
% c=round(linspace(1,size(cmap,1),size(I,3)));
% n1=4+round(100*rand);
% n2=104+round(100*rand);
% for k=1:size(I,3)
%     xscale=linspace(0,numel(Im(k).IntS)*p.n_L/p.n_s-1/2*p.n_L/p.n_s,numel(Im(k).IntS));
%     figure(n1);
%     subplot(ceil(sqrt(size(I,3))),round(sqrt(size(I,3))),k);
%     if ~isempty(Im(k).Fit)
%         plot(xscale, Im(k).IntS,'.',xscale, Im(k).Fit,'r');
%     elseif isempty(Im(k).Fit)
%         plot(xscale, Im(k).IntS,'.');
%     end
%     title(['Method ', num2str(k)]);
% 
%     figure(n2);
%     subplot(ceil(sqrt(size(I,3))),round(sqrt(size(I,3))),k);
%     hist(Im(k).Sh,linspace(min(mn),max(mx),32)); 
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor',cmap(c(k),:), 'EdgeColor','w');
%     title(['Method ', num2str(k)]);
%     legend(['Median: ' num2str(md(k)) '\pm' num2str(sd(k)) ', H: ' num2str(H(k,:))]);
%     [md(k), std(Im(k).Sh)], % Average and std of final sharpness
% end
% Ding 2013-07-24, end
% [av, md, std(Im(k).Sh)]
%% Save Relevant Data
% for i=1:k
%     Sav(i).Sh=Im(i).Sh;
% end
% Sav(1).p=p;
% Sav(1).t=t;
% 
% fname=name{1}; % Read file name and location
% tmp=find(fname=='\', 1, 'last' ); % Find last '\'
% dr=fname(1:tmp); % Location
% rname= [fname(tmp+1:end-4) '_Res_f' num2str(p.i_fn) '.mat']; % Append the name with 'Res'
% save([dr rname], 'Sav'); % Save results with the same directory
% % save(fullfile([dr rname]),'Sav')

% Temp
% figure; hold on;
% c=jet(8);
% ind=0;
% for i=1:8
%     ind=ind+1
%     plot(linspace(0,360, numel(Im(i).Sh)), circshift(Im(i).Sh,[-9,0]),'color',c(ind,:));
% end


% %% Temp
% xx=St.Angle; 
% xx(158:end)=pi+xx(158:end); 
% xx(450:end)=xx(450:end)+pi;
% c=jet(3);
% ind=0;
% figure; hold on;
% for i=[1,6,8]
%     ind=ind+1;
%     plot(xx,Im(i).Sh,'color',c(ind,:));
% end

% fiure; plot(av(1),'-o'); hold on; plot(av(2:13),'-rv'); hold on; plot(av(14:25),'-g^'); hold on; plot(av(26:end),'-s'); legend('original', 'R3', 'R5', 'R7');