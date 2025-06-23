function [ResCol,ResLin,ResPar,Rate_s,Rate_w,NImageCols,NImageLins,NImagePars] = constructArray_v2_5d(FlowSG4D_Outputs, opt, name, saveLocation, param, raw_header,acq_param)


SG_Reveal4D_Data = FlowSG4D_Outputs;
clear FlowSG4D_Outputs


ar_size = SG_Reveal4D_Data.ArraySize;
ar_ind = SG_Reveal4D_Data.Indices;
ar_ind = [ar_ind,(1:size(ar_ind,1)).'];
data = SG_Reveal4D_Data.Data;
rWeights = SG_Reveal4D_Data.rWeights;
rWeights = squeeze(single(rWeights));
    
%     ar_size(end) = 1;

% crop = max([round(ar_size(1) / 4), round((ar_size(1)-96)/2)]); % 42 for patient 3d cine
% crop = max([round(ar_size(1) / 4), round((ar_size(1)-384)/2)]);
crop1 = round(opt.crop(1)*ar_size(1));
crop2 = round(opt.crop(2)*ar_size(1));

% crop2_1 = round(opt.crop2(1)*ar_size(2)); %side crop murtaza
% crop2_2 = round(opt.crop2(2)*ar_size(2)); %side crop murtaza

asym_size = max(find(data(1:end-1,1,1)==0));
asym_percent = asym_size/size(data,1);

% 
% asym_size2 = max(find(data(1,1:end-1,1)==0)); %side crop murtaza
% asym_percent2 = asym_size2/size(data,2); %side crop murtaza

% ar_size(1) = ar_size(1)-crop*2;
ar_size(1) = ar_size(1)-(crop1+crop2);
% ar_size(2) = ar_size(2)-(crop2_1+crop2_2);



data = ifftshift(ifft(fftshift(data),[],1));
% data = data(crop+1:end-crop,:,:);
data = data(crop1+1:end-crop2,:,:);
data = ifftshift(fft(fftshift(data),[],1));

    %%
    if ~opt.flow
        ar_size(5) = 1;
    end
    for r_ind = 1:ar_size(end)
        k_space = (zeros(ar_size([1,2,3,4,5,6]),'single'));
        weights = (zeros(1,1,ar_size(3),ar_size(4),ar_size(5),ar_size(6),'single'));
        r1_rows = find(ar_ind(:,5)==r_ind);
%         ar_ind_r1 = ar_ind(r1_rows,:);
        ar_ind_r1 = ar_ind;
        
        replace = zeros(1,length(ar_ind_r1));
        counter = 0;
        counter2 = 0;
        % k_space(:,:,ar_ind(ind,1),ar_ind(ind,2),ar_ind(ind,4)) = data(:,:,ind);
        
        rWeights_tmp = rWeights(:,r_ind);
        for ind = 1:length(ar_ind_r1)
            
            if rWeights(ar_ind_r1(ind,6),r_ind) > weights(:,:,ar_ind_r1(ind,1),ar_ind_r1(ind,2),ar_ind_r1(ind,3),ar_ind_r1(ind,4))
                k_space(:,:,ar_ind_r1(ind,1),ar_ind_r1(ind,2),ar_ind_r1(ind,3),ar_ind_r1(ind,4)) = data(:,:,ar_ind_r1(ind,6));
                weights(:,:,ar_ind_r1(ind,1),ar_ind_r1(ind,2),ar_ind_r1(ind,3),ar_ind_r1(ind,4)) = rWeights_tmp(ar_ind_r1(ind,6));
%                 replace(ar_ind_r1(ind,6),r_ind) = 1;
%                 counter = counter +1;
%                 replaceW(counter,r_ind) = rWeights(ar_ind_r1(ind,6),r_ind);
            else
%                 counter2 = counter2 +1;
%                 skip(counter2,r_ind) = rWeights(ar_ind_r1(ind,6),r_ind);
            end
            
            
        end
    
        k_space = permute(k_space,[1,3,4,2,6,5]);
        kb = k_space(:,:,:,:,:,1);
        if opt.flow
            kx = k_space(:,:,:,:,:,2);
            ky = k_space(:,:,:,:,:,3);
            kz = k_space(:,:,:,:,:,4);
        else
            kx = kb;
            ky = kb;
            kz = kb;
        end
        clear k_space
        
        weights = repmat(weights, ar_size(1), ar_size(2), 1, 1, 1, 1);
        weights = permute(weights, [1,3,4,2,6,5]);
        weightsB = weights(:,:,:,:,:,1);
        if opt.flow
            weightsX = weights(:,:,:,:,:,2);
            weightsY = weights(:,:,:,:,:,3);
            weightsZ = weights(:,:,:,:,:,4);
        else
            weightsX = weightsB;
            weightsY = weightsB;
            weightsZ = weightsB;
        end
        clear weights

        sampB = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');
        sampX = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');
        sampY = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');
        sampZ = zeros(size(squeeze(kx(:,:,:,1,:))),'logical');

        sampB(abs(squeeze(kb(:,:,:,1,:)))>0) =1;
        sampX(abs(squeeze(kx(:,:,:,1,:)))>0) =1;
        sampY(abs(squeeze(ky(:,:,:,1,:)))>0) =1;
        sampZ(abs(squeeze(kz(:,:,:,1,:)))>0) =1;
    
        % Enforce Asymmetric Echo
        sampB(1:round(size(sampB,1)*asym_percent),:,:,:,:) = 0;
        sampX(1:round(size(sampX,1)*asym_percent),:,:,:,:) = 0;
        sampY(1:round(size(sampY,1)*asym_percent),:,:,:,:) = 0;
        sampZ(1:round(size(sampZ,1)*asym_percent),:,:,:,:) = 0;
        
        kb = bsxfun(@times,kb,permute(sampB,[1,2,3,5,4]));
        kx = bsxfun(@times,kx,permute(sampX,[1,2,3,5,4]));
        ky = bsxfun(@times,ky,permute(sampY,[1,2,3,5,4]));
        kz = bsxfun(@times,kz,permute(sampZ,[1,2,3,5,4]));
        
        weightsB = bsxfun(@times,weightsB,permute(sampB,[1,2,3,5,4]));
        weightsX = bsxfun(@times,weightsX,permute(sampX,[1,2,3,5,4]));
        weightsY = bsxfun(@times,weightsY,permute(sampY,[1,2,3,5,4]));
        weightsZ = bsxfun(@times,weightsZ,permute(sampZ,[1,2,3,5,4]));
        
        weightsB = squeeze(weightsB(:,:,:,1,:));
        weightsX = squeeze(weightsX(:,:,:,1,:));
        weightsY = squeeze(weightsY(:,:,:,1,:));
        weightsZ = squeeze(weightsZ(:,:,:,1,:));
        
        % ==========================================================================================
        % MAKE RESOLUTION ISOTROPIC (OR AS CLOSE AS POSSIBLE)
        % (This needs more work)
        % ==========================================================================================
        % Maximum allowable matrix size as contraint
        % ******************************************
        maxCols = 96;
        maxLins = 96;
        maxPars = 72;
        
        NImageCols = maxCols;
        NImageLins = maxLins;
        NImagePars = maxPars;
        
%         % ******************************************
%         % Current resolution
        ResCol = SG_Reveal4D_Data.param.ResCol;
        ResLin = SG_Reveal4D_Data.param.ResLin;
        ResPar = SG_Reveal4D_Data.param.ResPar;
        
        
        
    
        scanParam = SG_Reveal4D_Data.param;
        scanParam.totalbeats=opt.numofbeats;
        scanParam.totaltime= opt.totalscantime;
        scanParam.flow= opt.flow;
        scanParam.flow= opt.flow;
        scanParam.flow= opt.flow;
        scanParam.NImageCols = size(kb,1);
        scanParam.NImageLins = size(kb,2);
        scanParam.NImagePars = size(kb,3);
        scanParam.ResCol = ResCol;
        scanParam.ResLin = ResLin;
        scanParam.ResPar = ResPar;

        scanParam.rawfile= acq_param.rawfile;
        scanParam.ind=acq_param.ind;
        if exist('acq_param.scanner','var')
        scanParam.scanner=acq_param.scanner;
        end
        scanParam.BW=acq_param.BW;
        scanParam.TE=acq_param.TE;
        scanParam.TR=acq_param.TR;
        scanParam.FA=acq_param.FA;
        if(flow)
        scanParam.venc=acq_param.venc;
        end
        scanParam.Nx= acq_param.Nx;
        scanParam.Ny=acq_param.Ny;
        scanParam.Nz1=acq_param.Nz1;
        scanParam.Nz2=acq_param.Nz2;
        scanParam.FOVx=acq_param.FOVx;
        scanParam.FOVy=acq_param.FOVy;
        scanParam.FOVz=acq_param.FOVz;
        scanParam.resx=acq_param.resx;
        scanParam.resy=acq_param.resy;
        scanParam.resz=acq_param.resz;
        
        
        dirName = [saveLocation,'\',name];
        if ~exist(dirName,'dir')
            mkdir(dirName);
        end
        
                
        % estimate acceleration rate
        % do not include asymmetric echo in rate calculation
        Rate_s = (size(sampB,2)*size(sampB,3)*size(sampB,4)) / sum(sum(sum(sampB(round(size(sampB,1)/2),:,:,:),2),3),4);
        Rate_w = (size(sampB,2)*size(sampB,3)*size(sampB,4)) / sum(sum(sum(weightsB(round(size(sampB,1)/2),:,:,:).^2,2),3),4);
        
        scanParam.AccelerationRate_S = Rate_s;
        scanParam.AccelerationRate_W = Rate_w;

        name_r=[name,'_res',num2str(r_ind)];
        convert_binning_output(dirName, name_r, kb, kx, ky, kz, sampB, sampX, sampY, sampZ, weightsB, weightsX, weightsY, weightsZ, scanParam);

%         Rate_s = prod(size(sampB)) / sum(sampB(:));
%         Rate_w = prod(size(sampB)) / sum(weightsB(:).^2);
        disp(['Acceleration rate...']);
        disp(['Samples: ',num2str(Rate_s)]);
        disp(['Weighted: ',num2str(Rate_w)]);        
        


    end
    
    clear kx ky kz
    
%     k_space = padarray(k_space,[padLength1,padLength+padLength2,padLength3,0,0,0],0,'both');
    k_space_n = kb;
    k_space_n(abs(k_space_n)==0) = NaN;
    k_space_n = mean(k_space_n,5,'omitnan');
%     k_space_n = mean(k_space_n,6,'omitnan');
    k_space_n(isnan(k_space_n))=0;
    k_space_n = squeeze(k_space_n);
    im = ifft3_shift(k_space_n);
    im_c = im;
    im = sos_combine(im);

    %%
    figure;
    mn = min(abs(im(:)));
    mx = max(abs(im(:)))/5;
    imagesc(abs(im(:,:,round(size(im,3)/2))),[mn,mx]); axis('off','image'); colormap('gray')
    figure;
    imagesc(squeeze(abs(im(round(size(im,1)/2),:,:))),[mn,mx]); axis('off','image'); colormap('gray');
    figure;
    imagesc(squeeze(abs(im(:,round(size(im,2)/2),:))),[mn,mx]); axis('off','image'); colormap('gray');
    
    figure;
    nC = size(im_c);
    mn = min(abs(im_c(:)))/2; mx = max(abs(im_c(:)))/3;
    dim1 = 6; dim2 = 6;
    for i = 1 : size(im_c, 4)
        subplot(dim1, dim2, i)
        imagesc(abs(im_c(:,:,round(size(im,3)/2),i)),[mn, mx]); axis('off','image'); colormap('gray');
    end
        
        
        
        
        
        
    
end



