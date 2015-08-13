function [finalred, maxCircles, centers, radii] = fit_nuclei(nucfilename, numplanes)

%% Nathaniel I. Krefman 2014-05-06

% Loops through a stack identifying nuclei and finding their centroids

% Input variables:
%   nucfilename:      complete path to nucleus image
%   numplanes:        number of planes in the image

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   kernelsize_nuc & sigma_nuc
%   gausskernelsize_nuc & gausssigma_nuc
%   radius range and sensitivity of imfindcircles()

% Other file requirements:
%       medthresh.m
%       autogain_uint16.m

%%

% Open the nucleus image, filter it, then 
% calculate the centroids and approximate 
% radii of the nuclei.

% Open the nucleus image
image_nuc = openTiffStack(nucfilename);
% Convert the image to 8-bit
image_nuc = uint8(image_nuc);

% Pad the nucleus image to match the dimensions of the spot image
pad = 50;
[sizeY, sizeX, sizeZ] = size(image_nuc);
paddedY = sizeY+pad*2;
paddedX = sizeX+pad*2;

% Create blank padded image
finalred = zeros(paddedY,paddedX,sizeZ);

% Define the Laplacian of Gaussians (LOG) filter for 
% fitting the nuclei.
kernelsize_nuc = 20; % 20 seems good
sigma_nuc = 4.4; % 4.4 seems good
kernelnuc = -fspecial('log', kernelsize_nuc, sigma_nuc);

% Make a Gaussian filtered image
gausskernelsize_nuc = 9; % 9 finds my peaks well, 5 is more conservative (ignores dimmest spots)
gausssigma_nuc = 6; % 6 finds my peaks well
kernelgauss = fspecial('gaussian', gausskernelsize_nuc, gausssigma_nuc);
imtempgauss = imfilter(double(image_nuc(:,:,:)), kernelgauss, 'symmetric', 'conv');
finalred(pad+1:end-pad, pad+1:end-pad,:) = imtempgauss(:,:,:);

% Fit the nuclei using imfindcircles()

% Create containers to store centers and radii.
% maxCircles is to create a large enough container
% based on an overestimation of the number of cells
% expected per image.  (e.g. You'll probably never 
% see >400 nuclei with the Nikon C-MOS.)
maxCircles = 400; 
centers = zeros(maxCircles,2,numplanes) + 1;
radii = zeros(maxCircles,numplanes) + 1;

% Loop over each plane
for i=1:numplanes
    
    % Filter image with LOG (an edge-detection filter)
	imtemp = imfilter(double(image_nuc(:,:,i)), kernelnuc, 'symmetric', 'conv');
    
    % Resize the image after filtering to prevent edge artifacts
    padded = zeros(paddedY,paddedX);
    padded(pad+1:end-pad, pad+1:end-pad) = imtemp;
    log_nuc(:,:,i) = padded;
    
    % Find the circles in the image.
    % [15 20] pixels seems to be a good radius range.
    % 'Sensitivity' 0.955 seems good.
    [middles, radiuses] = imfindcircles(log_nuc(:,:,i), [15 20], 'Sensitivity', .955);
    
    % Figure out how many circles there were in the plane
    numCircles = size(middles);
    
    % Store the current circle data with the rest
    centers(1:numCircles(1),1:numCircles(2),i) = middles(:,:);
    % Need a correction factor for the radius, because 
    % imfindcircles() makes the circles too big
    window = 1.3; % probably need to optimize this
    radii(1:numCircles(1),i) = radiuses(:);
    radii(1:numCircles(1),i) = radii(1:numCircles(1),i)./window;

end

% Adjust the display range for the red image
finalred = double(finalred);
rawred = finalred;
% Threshold based on the median
finalred = medthresh(finalred);
% Convert to 16 bit with autogain
finalgreen = autogain_uint16(finalred);
gain = 2; % change gain for the display
finalred = finalred*gain;
