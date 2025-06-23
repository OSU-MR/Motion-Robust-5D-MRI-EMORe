function [] = convert_binning_output(dirName, name, kb, kx, ky, kz, sampB, sampX, sampY, sampZ, weightsB, weightsX, weightsY, weightsZ, scanParam)
        % Construct output array
        % [FE, PE, SL, CH, CPh, RPh, ENC]
        K(:,:,:,:,:,1) = kb;
        K(:,:,:,:,:,2) = kx;
        K(:,:,:,:,:,3) = ky;
        K(:,:,:,:,:,4) = kz;
        % samples
        % [FE, PE, SL, CPh, RPh, ENC]
        S(:,:,:,:,1) = sampB;
        S(:,:,:,:,2) = sampX;
        S(:,:,:,:,3) = sampY;
        S(:,:,:,:,4) = sampZ;
        % weights
        % [FE, PE, SL, CPh, RPh, ENC]
        W(:,:,:,:,1) = weightsB;
        W(:,:,:,:,2) = weightsX;
        W(:,:,:,:,3) = weightsY;
        W(:,:,:,:,4) = weightsZ;
       
        tmp_samp = squeeze(S(round(end/2),:,:,:,:));
        AccelerationRate = numel(squeeze(S(round(end/2),:,:,:,:)))/sum(tmp_samp(:));
        scanParam.AccelerationRate = AccelerationRate;

        D.kb = squeeze(K(:,:,:,:,:,1));
        D.kx = squeeze(K(:,:,:,:,:,2));
        D.ky = squeeze(K(:,:,:,:,:,3));
        D.kz = squeeze(K(:,:,:,:,:,4));
        D.sampB = squeeze(S(:,:,:,:,1));
        D.sampX = squeeze(S(:,:,:,:,2));
        D.sampY = squeeze(S(:,:,:,:,3));
        D.sampZ = squeeze(S(:,:,:,:,4));
        D.weightsB = squeeze(W(:,:,:,:,1));
        D.weightsX = squeeze(W(:,:,:,:,2));
        D.weightsY = squeeze(W(:,:,:,:,3));
        D.weightsZ = squeeze(W(:,:,:,:,4));
        D.scanParam = scanParam;
            
        save([dirName,'/',name], 'D','-v7.3'); 
end

