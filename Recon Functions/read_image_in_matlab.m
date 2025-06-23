clear;
% Parameters
filename = './pat6/pat6_1x1x1mm_409x409x201x5x8_snr40_fa45_fb.cpx'; % Replace with your file path
width = 409;
height = 409;
depth = 201;
channels = 5;
timepoints = 8;
dataType = 'float32'; % Assuming 32-bit float

% Open the file
fileID = fopen(filename, 'r');

% Check if the file opened successfully
if fileID == -1
    error('Failed to open file: %s', filename);
end

% Read the entire content of the file
data = fread(fileID, inf, dataType);

% Close the file
fclose(fileID);

% Calculate the total number of elements per frame
numElementsPerFrame = width * height * depth * channels * timepoints;

% Verify if the data length matches the expected number of elements
if mod(length(data), 2) ~= 0
    error('The data length is not consistent with complex data format.');
end

% Separate the real and imaginary parts
numElements = length(data) / 2;
realPart = data(1:2:end);
imaginaryPart = data(2:2:end);

% Verify if the number of elements match
if length(realPart) ~= numElementsPerFrame || length(imaginaryPart) ~= numElementsPerFrame
    error('The number of elements does not match the expected size.');
end

% Combine real and imaginary parts to form the complex image
complexData = complex(realPart, imaginaryPart);

% Reshape the data to [width, height, depth, channels, timepoints]
complexImage = reshape(complexData, [width, height, depth, channels, timepoints]);

% Display the size of the complex image to verify
disp(size(complexImage)); % Expected output: [409, 409, 201, 5, 8]

% Optionally, visualize a slice (e.g., the first timepoint and channel)
sliceIndex = 1;
 orthosliceViewer(abs(complexImage(:,:,:,1,1)));