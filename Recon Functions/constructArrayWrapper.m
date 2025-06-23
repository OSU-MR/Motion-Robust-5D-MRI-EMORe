function [ResCol,ResLin,ResPar,rate_samps,rate_weights,NImageCols,NImageLins,NImagePars] = constructArrayWrapper(FlowSG4D_Outputs, opt, name, saveLocation, param, raw_header)


% switch opt.rSoft
%     
%     % Basically, want to do svd on array after hardbinning, then switch to
%     % soft gating
%     case 0
%         [ResCol,ResLin,ResPar] = constructArray_2(FlowSG4D_Outputs, opt, name, saveLocation, param);
%         
%         
%     case 1
%         [ResCol,ResLin,ResPar,rate_samps,rate_weights,NImageCols,NImageLins,NImagePars] = constructArray_2(FlowSG4D_Outputs, opt, name, saveLocation, param);
        [ResCol,ResLin,ResPar,rate_samps,rate_weights,NImageCols,NImageLins,NImagePars] = constructArray_v2(FlowSG4D_Outputs, opt, name, saveLocation, param, raw_header);




end

