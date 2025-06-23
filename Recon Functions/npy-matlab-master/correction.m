% clear;
close all;
addpath(genpath('D:\MRI Data\Murtaza\MATLAB Codes\npy-matlab-master'));
[name,path] = uigetfile('D:\MRI Data\Murtaza\Saved Data\comparison\rest\_4D Flow CORe', 'Select data file');
outputs = importdata(fullfile(path, name)); % Read the data
sag=outputs.xHat;

data = readNPY('D:\MRI Data\Murtaza\MATLAB Codes\SCC\correction_map.npy');
%['Sli', 'Lin', 'Col']
correction_sag=flip(flip(flip(permute(data,[2,1,3]),3),2),1);

crop1 = ceil((size(correction_sag,1) -size(sag,1))./2);
crop2 = floor((size(correction_sag,1) -size(sag,1))./2);

correction_sag=correction_sag(crop1+1:end-crop2,:,:,:);



corrected_sag=sag.*correction_sag;

ax=permute(sag, [2,3,1,4]);
correction_ax=permute(correction_sag, [2,3,1,4]);
corrected_ax=permute(corrected_sag, [2,3,1,4]);

cor=permute(sag, [1,3,2,4]);
correction_cor=permute(correction_sag, [1,3,2,4]);
corrected_cor=permute(corrected_sag, [1,3,2,4]);

% sliceViewer(correction_sag)
% clip=0.5;
% figure;
% for slice=1:72
%     subplot(331)
%     imagesc(sag(:,:,slice),[min(sag(:)) clip.*max(sag(:))]);colormap("gray");colorbar;
%     title("Uncorrected");
%     subplot(332)
%     imagesc(correction_sag(:,:,slice),[min(correction_sag(:)) clip.*max(correction_sag(:))]);colormap("gray");colorbar;
%     title("Map");
%     subplot(333)
%     imagesc(corrected_sag(:,:,slice),[min(corrected_sag(:)) clip.*max(corrected_sag(:))]);colormap("gray");colorbar;
%     title("Corrected"); 
% 
%     subplot(334)
%     imagesc(ax(:,:,slice),[min(sag(:)) clip.*max(sag(:))]);colormap("gray");colorbar;
%     subplot(335)
%     imagesc(correction_ax(:,:,slice),[min(correction_ax(:)) clip.*max(correction_ax(:))]);colormap("gray");colorbar;
%     subplot(336)
%     imagesc(corrected_ax(:,:,slice),[min(corrected_ax(:)) clip.*max(corrected_ax(:))]);colormap("gray");colorbar;
% 
%     subplot(337)
%     imagesc(cor(:,:,slice),[min(sag(:)) clip.*max(sag(:))]);colormap("gray");colorbar;
%     subplot(338)
%     imagesc(correction_cor(:,:,slice),[min(correction_cor(:)) clip.*max(correction_cor(:))]);colormap("gray");colorbar;
%     subplot(339)
%     imagesc(corrected_cor(:,:,slice),[min(corrected_cor(:)) clip.*max(corrected_cor(:))]);colormap("gray");colorbar;
%     pause();
% end
% 
