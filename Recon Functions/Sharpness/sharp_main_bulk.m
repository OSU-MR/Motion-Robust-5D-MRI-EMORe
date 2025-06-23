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
% of this software or its fitness for any particular purpose."
%-------------------------------------------------------------------------- 
% The code is written by Rizwan Ahmad (Rizwan.Ahmad@osumc.edu).  


clear;
clc;
close all;


%% User defined parameters
p.e_sig   = 0.7; % Width of smoothing filter for edge detection
p.e_th    = 0.04; % Thresholding for edge detection
p.e_tl    = p.e_th*0.5;
p.c_L     = 3; % Length used to compute curvature for normals
p.e_L     = p.c_L+1; % Discard all edges below length e_L
p.n_L     = 8; % Length (in pixels) of lines normal to the selected edges
p.n_s     = 6*p.n_L+1; % Number of samples along each normal line; keep it odd
p.h_dc    = 0.02; % Histogram parameters (dc)
p.h_ln    = 1024; % Number of bins in the histogram
p.h_sig   = p.h_ln/3; % std of histogram
p.i_fn    = 1; % Frame number to be considered for analysis
p.i_ref   = 2;


%% Add path
% addpath('.\common'); % Helper functions
% addpath('.\NFFT');


%% Load Image
for d=1:51
    if d<10
        name=['C:\Users\Rizwan\Research\SHARPNESS\Sharpness_Cine\SharpnessR5\dataForSharpness_0' num2str(d) '.mat'];
    elseif d>=10
        name=['C:\Users\Rizwan\Research\SHARPNESS\Sharpness_Cine\SharpnessR5\dataForSharpness_' num2str(d) '.mat'];
    end
%     for i=1:size(name,2)
    fid = fopen(char(name)); 
    Raw=load(char(name));
    Raw = struct2cell(Raw);
    if numel(Raw)>1
        p.i_fn=Raw{2}(2);% ??? if empty
    end
    Raw=Raw{1};
    fclose(fid);
    if numel(size(Raw))==4
        I=squeeze(Raw(:,:,p.i_fn,:));
    elseif numel(size(Raw))==2 || numel(size(Raw))==3
        I=Raw;
    else
        error('Uploaded MAT file should be an array of size 2D to 4D');
    end
%     end
    I=double(I);
    clear Raw;
    fclose('all');

    %% Or create phantom
    % load I;



    %% Edge Detection
    hgram=hshape(p);
    if strcmp(p.i_ref,'all')
        Ic=mean(I,3);
        Ic=histeq(Ic/max(Ic(:)),hgram);
    else
        Ic=I(:,:,p.i_ref);
        Ic=histeq(Ic/max(Ic(:)),hgram);
    end
    for i=1:size(I,3)
    %     I(:,:,i)=histeq(I(:,:,i)/max(max(I(:,:,i))),hgram);
    end

    %% Edge Detection
    E = edge(Ic,'canny',[p.e_tl p.e_th],p.e_sig);


    %% Remove branches
    E_b = brch_prune(E,'A');
    E_b = len_thresh(E_b,p.e_L);


    %% Plotting
    figure(1);
    for i=1:size(I,3);
        subplot(floor(sqrt(size(I,3))),ceil(sqrt(size(I,3)))+1,i); 
        imagesc(I(:,:,i)); colormap(gray); axis('image','off');
        title(['Method ', num2str(i)]);
    end

    hfig2=figure(2);  
        subplot(221); imagesc(Ic); colormap(gray); axis('image','off'); title('composite image');
        subplot(222); imagesc(E_b); colormap(gray); axis('image','off'); title('edges');
        set(hfig2,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 
    hfig3=figure(3); 
        imagesc(E_b); colormap(gray); axis('image','off'); title('using mouse, select appropriate edges');
        set(hfig3,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); 

    E_bs=bwselect(E_b,8);
    figure(2);  subplot(223); imagesc(E_bs); colormap(gray); axis('image','off'); title('user selected edges');
    %             subplot(224); imagesc(slct_e); colormap(gray); axis('square','off');



    %% Morphological filtering
    E_bsm   = morph_mask(E_b).*E_bs;
    % E_bsm = len_thresh(E_bsm,p);


    %% Remove pixels close to the edge
    E_bsm = rect_mask(size(Ic),p.n_L/2).*E_bsm;
    E_bsm = len_thresh(E_bsm,p.e_L);
    figure(2);  subplot(224); imagesc(E_bsm); colormap(gray); axis('image','off'); title('edges for sharpness measurement');


    %% Sorting Edge Indices
    St = edge_sort(E_bsm);


    %% Find normals to local curvature
    St = edge_norm(St,size(Ic),p);


    %% Curve fitting
    av=zeros(size(I,3),1);
    mx=zeros(size(I,3),1);
    mn=zeros(size(I,3),1);
    t=zeros(size(I,3)-1, size(I,3)-1);
    for k=1:size(I,3)
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
                Im(k).Sh=[Im(k).Sh;Im(k).Data(i).Param(:,3)];
                Im(k).IntS=[Im(k).IntS;Im(k).Data(i).IntS];
                Im(k).Fit=[Im(k).Fit;Im(k).Data(i).Fit];
            end
        end 
        if k>1
            for g=1:1:k-1
                t(g,k-1)=ttest2(Im(g).Sh, Im(k).Sh,0.05,'both','unequal');
            end
        end
        av(k)=mean(Im(k).Sh);
        mn(k)=min(Im(k).Sh);
        mx(k)=max(Im(k).Sh);
    end
    t=[t; zeros(1,k-1)]+[zeros(1,k-1); t'];

    %% Plotting
    load cmap;
    c=round(linspace(1,size(cmap,1),size(I,3)));
    n1=4+round(100*rand);
    n2=104+round(100*rand);
    for k=1:size(I,3)
        xscale=linspace(0,numel(Im(k).IntS)*p.n_L/p.n_s-1/2*p.n_L/p.n_s,numel(Im(k).IntS));
        figure(n1);
        subplot(size(I,3),1,k);
        plot(xscale, Im(k).IntS,'.',xscale, Im(k).Fit,'r');
        title(['Method ', num2str(k)]);

        figure(n2);
        subplot(size(I,3),1,k);
        hist(Im(k).Sh,linspace(min(mn),max(mx),16)); 
        h = findobj(gca,'Type','patch');
        set(h,'FaceColor',cmap(c(k),:), 'EdgeColor','w');
        title(['Method ', num2str(k)]);
        legend(['Mean: ' num2str(av(k)) '\pm' num2str(std(Im(k).Sh)) ', H: ' num2str(t(k,:))]);
    end


    %% Save Relevant Data
    for i=1:k
        Sav(i).Sh=Im(i).Sh;
    end
    Sav(1).p=p;
    Sav(1).t=t;

    fname=name; % Read file name and location
    tmp=find(fname=='\', 1, 'last' ); % Find last '\'
    dr=fname(1:tmp); % Location
    rname= [fname(tmp+1:end-4) '_Res_f' num2str(p.i_fn) '.mat']; % Append the name with 'Res'
    save([dr rname], 'Sav'); % Save results with the same directory
    % save(fullfile([dr rname]),'Sav')
    close all;
end