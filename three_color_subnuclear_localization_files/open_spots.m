function [MIP_spot, spot_projection, imLOG_spot, imgauss_spot, finalgreen] = open_spots(imG, numplanes, offsetGreenYdisplay, offsetGreenXdisplay)

%% Nathaniel I. Krefman 2014-05-06

% Opens the green image, creates a max-intensity projection, 
% Gaussian- and LOG-filtered images, a blank image to store 
% a spot projection, and a contrast-adjusted image for 
% supervised analysis.

% Input variables:
%   imG:                    green image
%   numplanes:              number of planes in the image
%   offsetGreenYdisplay:    Y offset to match green to red
%   offsetGreenXdisplay:    X offset to match green to red

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   kernelsize_spot & sigma_spot (NOTE: Used twice -- for Gaussian filter
%   and LOG filter.)
%   threshold & im_max (NOTE: These variables may be partly redundant with
%   adjRangeAndBitDepth.m)

% Dependencies:
%   openTiffStack.m
%   autogain_uint16.m

%%

% Open the spots image
image_spot = openTiffStack(imG); % Need to have openTiffStack.m
% Convert the image to 8-bit
image_spot = uint8(image_spot);

% Pad the image before correcting for color offsets
pad = 50;

% Create a blank padded image
[sizeY, sizeX, sizeZ] = size(image_spot);
paddedY = sizeY+pad*2;
paddedX = sizeX+pad*2;
padded = zeros(paddedY,paddedX);

% Define the Laplacian of Gaussians (LOG) filter for spots
kernelsize_spot = 6; % 6 finds my peaks well, 5 is more conservative (ignores dimmest spots)
sigma_spot = 4; % 4 finds my peaks well
kernelspot = -fspecial('log', kernelsize_spot, sigma_spot);

% Define the Gaussian filter for spots
kernelsize_spot = 5; % 5 looks good
sigma_spot = 2; % 2 looks good
kernelgauss = fspecial('gaussian', kernelsize_spot, sigma_spot);

%%

% Create a maximum intensity projection of the spots image
padded(pad+1+offsetGreenYdisplay:end-pad+offsetGreenYdisplay, pad+1+offsetGreenXdisplay:end-pad+offsetGreenXdisplay) = max(double(image_spot),[],3);
MIP_spot = padded;
MIP_spot = autogain_uint16(MIP_spot);

% Create a matrix to store a projection of all the spots
spot_projection = zeros(size(padded));

% Filter then pad the stack using LOG & Gaussian filters

% Create blank padded image
padded3D = zeros(paddedY,paddedX,sizeZ);

% Make a LOG filtered image
imtemp = imfilter(double(image_spot(:,:,:)), kernelspot, 'symmetric', 'conv');
padded3D(pad+1+offsetGreenYdisplay:end-pad+offsetGreenYdisplay, pad+1+offsetGreenXdisplay:end-pad+offsetGreenXdisplay,:) = imtemp(:,:,:);
imLOG_spot(:,:,:) = padded3D;

% Make a Gaussian filtered image
imtempgauss = imfilter(double(image_spot(:,:,:)), kernelgauss, 'symmetric', 'conv');
padded3D(pad+1+offsetGreenYdisplay:end-pad+offsetGreenYdisplay, pad+1+offsetGreenXdisplay:end-pad+offsetGreenXdisplay,:) = imtempgauss(:,:,:);
imgauss_spot = padded3D;

% Make blank stack to store final spot images in
finalgreen = zeros(size(imgauss_spot));
