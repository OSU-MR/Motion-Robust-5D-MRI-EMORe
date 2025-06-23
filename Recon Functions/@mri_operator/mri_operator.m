classdef mri_operator
    %PMRI_OP_3D_T Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        frame_size;     % The Size of the image e.g. [256,256]
        N;              % Input Dimension Size
        M;              % Output Dimension Size
        maps;           % The coil Sensitivity maps
        compute;        % How to compute the wavelet transform
        uniform_var;
        mask_patterns;
        mask_weights;   % Weights for weighted least squares
    end
    
    properties(Access = protected, Hidden = true)
        
        nddwt;
        s;              % 2D wavelet signiture matrix 
    end
    
    methods
        % Object Constructor
        function obj = mri_operator(maps,samp,varargin)
                       
            % Set Default Options
            obj.compute = 'mat';
            obj.uniform_var = 0;
            obj.mask_weights = [];
            weights = [];
            
            % Copy any optional inputs 
            if mod(length(varargin),2)
                error('Optional inputs must come in pairs')
            end
            for ind = 1:2:length(varargin)
                switch lower(varargin{ind})
                    case 'uniform_var'
                        obj.uniform_var = varargin{ind+1};
                    case 'compute'
                        obj.compute = varargin{ind+1};
                    case 'weights'
                        weights = varargin{ind+1};
                    otherwise
                        warning(sprintf('Unknown optional input #%d ingoring!',ind))
                end
            end
            
            obj.maps = maps;
            % obj.frame_size(1) = size(maps,1);
            % obj.frame_size(2) = size(maps,2);
            % obj.frame_size(3) = size(maps,3);
            % obj.Q = size(samp,4);
            obj.frame_size = size(samp);

            
            % Find the size of the measured Data
            ndims_samp = ndims(samp);
            patterns = permute(samp,[1:3,ndims_samp+1,4:ndims_samp]);
            reps = ones(1, ndims_samp);
            reps(4) = size(maps,4);
            patterns = repmat(patterns,reps);
            obj.mask_patterns = find(patterns==1);
            obj.M = length(obj.mask_patterns);
            
            % Find the size of the measured Data
            if ~isempty(weights)
                weights = permute(weights,[1:3,ndims_samp+1,4:ndims_samp]);
                weights = repmat(weights,reps);
%                warning('look at this again')
                obj.mask_weights = weights(patterns==1);
            end
            
            % Find Image size
            obj.N = prod(obj.frame_size);
            
            % Get Frobenius Norm
            % obj.Fro2 = 1/(obj.N*size(obj.maps,4));
            
            % if ~obj.uniform_var
            %     obj.CSq = abs(obj.maps.^2)/(obj.N/obj.Q);
            %     obj.CSqTr = obj.CSq(:,:,:,:,1);
            %     obj.CSq = obj.CSq(:,:,:,:,1);
            %     obj.CSq = reshape(obj.CSq,[size(obj.CSq,1)*size(obj.CSq,2)*size(obj.CSq,3),size(obj.CSq,4)]);
            %     obj.CSq = obj.CSq.';
            % else
            %     obj.CSq = [];
            % end
            
%             % Check to make sure some inputs are the right size
%             if sum(size(samp(:,:,1)) ~= obj.frame_size) >=1
%                 error('Image size and sampling pattern must be the same size')
%             elseif size(maps,1) ~= obj.frame_size(1) || size(maps,2) ~= obj.frame_size(2)
%                 error('Image size and Coils must be the same size')
%             end
           

            
            if strcmpi(obj.compute,'gpu_off') || strcmpi(obj.compute,'gpu')
                obj.maps = gpuArray(obj.maps);
                % obj.CSq = gpuArray(obj.CSq);
                % obj.CSqTr = gpuArray(obj.CSqTr);
                obj.uniform_var = gpuArray(obj.uniform_var);
                obj.mask_patterns = gpuArray(obj.mask_patterns);
                obj.mask_weights = gpuArray(obj.mask_weights);
%                 obj.M = gpuArray(obj.M );
%                 obj.N = gpuArray(obj.N ); 
%                 obj.frame_size = gpuArray(obj.frame_size);
%                 obj.Q = gpuArray(obj.Q);
            end
        end

        % Return the Dimensions of the Matrix
        function [m,n] = size(obj)
            n = obj.N;	
            m = obj.M;
        end
        
        % Forward 2D t Operator
        y = mult(obj,x)

        % Hermitian 2D t Operator
        y = multTr(obj,x)

    end
    
end





