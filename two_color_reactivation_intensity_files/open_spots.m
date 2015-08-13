function [MIP_spot, spot_projection, imLOG_spot, imgauss_spot, finalgreen] = open_spots(imG, numplanes, offsetGreenYdisplay, offsetGreenXdisplay)

%% Nathaniel I. Krefman 2014-05-06

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   kernelsize_spot & sigma_spot (NOTE: Used twice -- for Gaussian filter
%   and LOG filter.)
%   threshold & im_max (NOTE: These variables may be partly redundant with
%   adjRangeAndBitDepth.m)

% Dependencies:
%   openTiffStack.m

%%

% Open the spots image
image_spot = openTiffStack(imG);
% Convert the image to 8-bit
image_spot = uint8(image_spot);

% Pad the image before correcting for color offsets
pad = 50;

% Create a blank padded image
[sizeY, sizeX, sizeZ] = size(image_spot);
paddedY = sizeY+pad*2;
paddedX = sizeX+pad*2;
padded = zeros(paddedY,paddedX);

% Define the Gaussian filter for spots
kernelsize_spot = 5; % 5 looks good
sigma_spot = 2; % 2 looks good
kernelgauss = fspecial('gaussian', kernelsize_spot, sigma_spot);

% Define the Laplacian of Gaussians (LOG) filter for spots
kernelsize_spot = 7; % 7 finds my peaks well, 5 is more conservative (ignores dimmest spots)
sigma_spot = 4; % 4 finds my peaks well
kernelspot = -fspecial('log', kernelsize_spot, sigma_spot);

%%

% Create a maximum intensity projection of the spots image
padded(pad+1+offsetGreenYdisplay:end-pad+offsetGreenYdisplay, pad+1+offsetGreenXdisplay:end-pad+offsetGreenXdisplay) = max(double(image_spot),[],3);
MIP_spot = padded;

% make sure the image is a double
MIP_spot = double(MIP_spot);

% set a threshold
im_min = 0;
MIP_spot = MIP_spot - im_min;
threshold = 0;
MIP_spot(find(MIP_spot < threshold)) = threshold;

% select a max instensity (im_max)
im_max = 70;

% adjust the display range for uint16 
% (i.e. 65535 values), setting any 
% pixel over im_max to saturated
MIP_spot = 65535*MIP_spot/im_max;

% turn the image into a uint16 array
MIP_spot = uint16(MIP_spot);

%%

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
