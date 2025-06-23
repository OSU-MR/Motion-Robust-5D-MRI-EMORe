function [selectedGPUIndex, gpuInfo] = selectGPUInteractive()
%SELECTGPUINTERACTIVE  List available GPUs and let user pick one.
%
%   [idx, info] = selectGPUInteractive()
%   • Queries gpuDeviceCount and gpuDevice()
%   • Prints a table of Name / Mem. Usage / (estimated) Util.
%   • Prompts for a valid index
%   • Calls gpuDevice(idx)
%
%   Outputs:
%     idx    – integer index of the chosen GPU
%     info   – 1×N cell array of structs with fields:
%                Name, TotalMemoryGB, AvailableMemoryGB, …
%
%   (c) Syed M. Arshad, 2025

    numGPUs = gpuDeviceCount;
    gpuInfo = cell(1, numGPUs);    % always initialize

    if numGPUs == 0
        error('No CUDA-enabled GPUs found.');
    end

    if numGPUs == 1
        % Single GPU: auto-select it
        fprintf('Only one GPU detected. Using GPU 1 by default.\n');
        g = gpuDevice(1);

        % Fill info struct
        totalGB     = g.TotalMemory/(1024^3);
        availGB     = g.AvailableMemory/(1024^3);
        usedGB      = totalGB - availGB;
        utilPercent = min(100,round((usedGB/totalGB)*100));

        gpuInfo{1} = struct(...
            'Name',            g.Name, ...
            'TotalMemoryGB',   totalGB, ...
            'AvailableMemoryGB',availGB, ...
            'UsedMemoryGB',    usedGB, ...
            'Utilization',     utilPercent, ...
            'ComputeCapability',g.ComputeCapability, ...
            'DriverVersion',   g.DriverVersion ...
        );
        selectedGPUIndex = 1;
        reset(g);  % clear state for subsequent recon
        
    else
        % Multiple GPUs: list them and ask user
        fprintf('Multiple GPUs detected (%d):\n', numGPUs);
        fprintf('-----------------------------------------------------------\n');
        fprintf('| GPU | %-35s |   Mem (GB)   | Util |\n', 'Name');
        fprintf('|-----|-------------------------------------|--------------|------|\n');

        % Gather info from each device
        for i = 1:numGPUs
            try
                g = gpuDevice(i);
                totalGB     = g.TotalMemory/(1024^3);
                availGB     = g.AvailableMemory/(1024^3);
                usedGB      = totalGB - availGB;
                utilPercent = min(100,round((usedGB/totalGB)*100));

                fprintf('| %3d | %-35s | %5.1f/%5.1f | %3d%% |\n', ...
                        i, g.Name, usedGB, totalGB, utilPercent);

                gpuInfo{i} = struct(...
                    'Name',            g.Name, ...
                    'TotalMemoryGB',   totalGB, ...
                    'AvailableMemoryGB',availGB, ...
                    'UsedMemoryGB',    usedGB, ...
                    'Utilization',     utilPercent, ...
                    'ComputeCapability',g.ComputeCapability, ...
                    'DriverVersion',   g.DriverVersion ...
                );
                reset(g);

            catch ME
                fprintf('| %3d | (error retrieving info)           |              |      |\n', i);
                warning('GPU %d info error: %s', i, ME.message);
                gpuInfo{i} = [];
            end
        end
        fprintf('-----------------------------------------------------------\n');

        % Prompt for selection
        while true
            idx = input(sprintf('Select GPU index (1–%d): ', numGPUs));
            if isnumeric(idx) && isscalar(idx) && idx>=1 && idx<=numGPUs ...
               && ~isempty(gpuInfo{idx})
                selectedGPUIndex = idx;
                gpuDevice(idx);
                fprintf('Selected GPU %d: %s\n', idx, gpuInfo{idx}.Name);
                break
            else
                fprintf('Invalid selection. Please try again.\n');
            end
        end
    end
end
