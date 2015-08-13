function [MIP_SAC, imgauss_SAC, finalred] = open_SAC(imR)
%% Nathaniel I. Krefman 2014-05-06

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   kernelsize_SAC & sigma_SAC
%   threshold & im_max (NOTE: These variables may be partly redundant with
%   adjRangeAndBitDepth.m)

% Dependencies:
%   openTiffStack.m
%   autogain_uint16.m

%%

% Open the spindle SAC image
image_SAC = openTiffStack(imR); % Need to have openTiffStack.m
% Convert the image to 8-bit
image_SAC = uint8(image_SAC);

% Pad the spindle SAC image to match the 
% dimensions of the other images
pad = 50;
[sizeY, sizeX, sizeZ] = size(image_SAC);
paddedY = sizeY+pad*2;
paddedX = sizeX+pad*2;
padded3D = zeros(paddedY,paddedX,sizeZ);
padded = zeros(paddedY,paddedX);

% Define the Gaussian filter for the SACs
kernelsize_SAC = 5; % 6 finds my peaks well, 5 is more conservative (ignores dimmest spots)
sigma_SAC = 2; % 4 finds my peaks well
kernelgauss = fspecial('gaussian', kernelsize_SAC, sigma_SAC);

% Make a Gaussian filtered image
imtempgauss = imfilter(double(image_SAC(:,:,:)), kernelgauss, 'symmetric', 'conv');
padded3D(pad+1:end-pad, pad+1:end-pad,:) = imtempgauss(:,:,:);
imgauss_SAC = padded3D;

%%

% Convert to a uint16 image
% with a fixed display range

% make sure the image is a double
finalred = double(imgauss_SAC(:,:,:));

% set a threshold
im_min = 50;
finalred = finalred - im_min;
threshold = 0;
finalred(find(finalred < threshold)) = threshold;

% select a max instensity (im_max)
im_max = 160;

% adjust the display range for uint16 
% (i.e. 65535 values), setting any 
% pixel over im_max to saturated
finalred = 65535*finalred/im_max;

% turn the image into a uint16 array
finalred = uint16(finalred);

%%

% Create a maximum intensity projection of the spots image
padded(pad+1:end-pad, pad+1:end-pad) = max(double(image_SAC),[],3);
MIP_SAC = padded;

% Convert to a uint16 image
% with a fixed display range

% make sure the image is a double
MIP_SAC = double(MIP_SAC);

% set a threshold
im_min = 60;
MIP_SAC = MIP_SAC - im_min;
threshold = 0;
MIP_SAC(find(MIP_SAC < threshold)) = threshold;

% select a max instensity (im_max)
im_max = 130;

% adjust the display range for uint16 
% (i.e. 65535 values), setting any 
% pixel over im_max to saturated
MIP_SAC = 65535*MIP_SAC/im_max;

% turn the image into a uint16 array
MIP_SAC = uint16(MIP_SAC);

MIP_SAC = autogain_uint16(MIP_SAC);