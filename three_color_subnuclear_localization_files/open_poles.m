function [MIP_pole, imgauss_pole, finalblue] = open_poles(imB, offsetBlueYdisplay, offsetBlueXdisplay)

%% Nathaniel I. Krefman 2014-05-06

% Opens the blue image, creates a max-intensity projection, 
% Gaussian-filtered image, and a contrast-adjusted image for 
% supervised analysis.

% Input variables:
%   imB:                    blue image
%   offsetBlueYdisplay:     Y offset to match blue to red
%   offsetBlueXdisplay:     X offset to match blue to red

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   kernelsize_pole & sigma_pole
%   subtract & gain (NOTE: These variables may be partly redundant with
%   adjRangeAndBitDepth.m)

% Dependencies:
%   openTiffStack.m
%   autogain_uint16.m

%%

% Open the spindle pole image
image_pole = openTiffStack(imB); % Need to have openTiffStack.m
% Convert the image to 8-bit
image_pole = uint8(image_pole);

% Pad the spindle pole image to match the 
% dimensions of the other images
pad = 50;
[sizeY, sizeX, sizeZ] = size(image_pole);
paddedY = sizeY+pad*2;
paddedX = sizeX+pad*2;
padded3D = zeros(paddedY,paddedX,sizeZ);
padded = zeros(paddedY,paddedX);

% Define the Gaussian filter for the poles
kernelsize_pole = 6; % 6 finds my peaks well, 5 is more conservative (ignores dimmest spots)
sigma_pole = 4; % 4 finds my peaks well
kernelgauss = fspecial('gaussian', kernelsize_pole, sigma_pole);

% Make a Gaussian filtered image
imtempgauss = imfilter(double(image_pole(:,:,:)), kernelgauss, 'symmetric', 'conv');
padded3D(pad+1+offsetBlueYdisplay:end-pad+offsetBlueYdisplay, pad+1+offsetBlueXdisplay:end-pad+offsetBlueXdisplay,:) = imtempgauss(:,:,:);
imgauss_pole = padded3D;

% Fix the blue image
finalblue = double(imgauss_pole(:,:,:));
finalblue = autogain_uint16(finalblue);
subtract = 20000; % adjust for each image?
finalblue = (finalblue - subtract);
gain = 3;
finalblue = finalblue.*gain;

% Create a maximum intensity projection of the spots image
padded(pad+1+offsetBlueYdisplay:end-pad+offsetBlueYdisplay, pad+1+offsetBlueXdisplay:end-pad+offsetBlueXdisplay) = max(double(image_pole),[],3);
MIP_pole = padded;
MIP_pole = autogain_uint16(MIP_pole);