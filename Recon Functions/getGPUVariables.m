% Get a list of all variables in the current workspace.

% Get the current GPU device.
try
    currentGPU = gpuDevice();
catch ME
    if strcmp(ME.identifier, 'parallel:gpu:device:DeviceNotFound')
        error('No CUDA-enabled GPU found.  Check your drivers and installation.');
    else
        rethrow(ME); % Re-throw other errors
    end
end

% Display total and available memory.
totalMemoryGB = currentGPU.TotalMemory / (1024^3);
availableMemoryGB = currentGPU.AvailableMemory / (1024^3);
fprintf('Total GPU Memory: %.2f GB\n', totalMemoryGB);
fprintf('Available GPU Memory: %.2f GB\n', availableMemoryGB);
fprintf('-----------------------------------\n');


% Get a list of all variables in the current workspace.
allVars = whos;

% Initialize an empty cell array to store the names of GPU variables.
gpuVars = {};
totalUsedMemoryGB = 0; % Keep track of total memory used by GPU vars

% Iterate through all variables.
for i = 1:length(allVars)
    % Check if the 'class' of the variable is 'gpuArray'.
    if strcmp(allVars(i).class, 'gpuArray')
        % If it is, add the variable name to the gpuVars cell array.
        gpuVars{end+1} = allVars(i).name; %#ok<AGROW>

        % Calculate and display the memory used by this variable (in GB).
        variableMemoryBytes = allVars(i).bytes;
        variableMemoryGB = variableMemoryBytes / (1024^3);
        fprintf('Variable: %-20s  Size: %.2f GB\n', allVars(i).name, variableMemoryGB);

        totalUsedMemoryGB = totalUsedMemoryGB + variableMemoryGB;
    end
end

fprintf('-----------------------------------\n');
fprintf('Total Memory Used by GPU Variables: %.2f GB\n', totalUsedMemoryGB);

if isempty(gpuVars)
    fprintf('No variables found on the GPU.\n');
end


