function final_results = main(dir, stack, imG, imR)

%% Nathaniel I. Krefman 2014-05-06

% Summary:

% Script for identifying detached centromeres and spindle peaks in
% cells following centromere reactivation, and then quantifying the 
% Mad1 intensity. Cells should have GFP-TUB1 and GAL-CEN-TetOs marked 
% with TetR-GFP and Mad1-3xmCherry.

% Adjust the min/max range of the images using adjRangeAndBitDepth.m 
% or ImageJ before running.

% We generally call this script using run_main.m, although it 
% can also be called directly from the Command Window. The 
% variables in run_main.m need to be updated every time you run 
% the code for a new image stack.

% Most of the code in this program does not need to be updated 
% at runtime, except to optimize for localization of a new spot 
% type, or when using new microscope or illumination settings, 
% for example.

% This script saves several output images, including data from every 
% cell. Be sure you have plenty of hard-drive space before running
% We delete the output images after we double check our 
% classifications due to digital storage constraints.

% Input variables:
%   dir:        Directory for storing the results.
%               The directory must already exist!
%   stack:      Reference number for the image being analyzed.
%               Will be stored in the final output table.
%               If you compile all your results from each timepoint
%               into one Excel spreadsheet, this is useful for sorting 
%               the data by image.
%   imG & imR:  green (CEN/MT) & red (MAD1) images
%               Give the complete paths.

% Dependencies (in order of first call):
%	open_spots.m
%	openTiffStack.m
%	open_SAC.m
%	autogain_uint16.m
%	find_centroids.m
%	link_centroids.m
%	find_peak_intensities.m

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   numplanes
%   offsetYprecise & offsetXprecise
%   pixelsize_nm
%   status

% Each function called from main.m also lists variables that 
% may require optimization for your set-up.

% Other parameters may require adjustment for different microscopes. 
% For instance, if the pixel size deviates substantially from 
% 65 nm, several distance parameters will need alteration.

%%

hold off
close all

% Input the number of z-planes in your stack.
numplanes = 28;

% Input the offsets calculated between the images.
% We recommend calculating the offset by aligning bead
% images (e.g. TetraSpeck beads (Life Technologies, Grand 
% Island, NY).
% Here, we shift the green image to match the red image.
% NOTE: This will require some playing to be sure you have 
% the signs correct, and you don't have Y and X backwards.
offsetYprecise = -4.3470;
offsetXprecise = -0.5360;

% Calculate the offset for displaying aligned images
% without interpolation (rounded to nearest pixel).
offsetYdisplay = round(offsetYprecise);
offsetXdisplay = round(offsetXprecise);

% Calculate the remaining subpixel offset.
% We use these values to correct our green:red distance 
% measurements much later!
offsetYremainder = offsetYprecise - offsetYdisplay;
offsetXremainder = offsetXprecise - offsetXdisplay;

%%

% Open the green image, create a
% max-intensity projection, Gaussian-
% and LOG-filtered images, a blank 
% image to store a spot projection,
% and a contrast-adjusted image for 
% supervised analysis.

[MIP_spot, spot_projection, imLOG_spot, imgauss_spot, finalgreen] = open_spots(imG, numplanes, offsetGreenYdisplay, offsetGreenXdisplay);

'Opened the spot image and filtered it.'

%%

% Open the red image, create a
% max-intensity projection, Gaussian-
% filtered images, and a contrast-
% adjusted image for supervised analysis.

[MIP_SAC, imgauss_SAC, finalred] = open_SAC(imR);

'Opened the SAC image and filtered it.'

%%

% Loop through stack identifying spots, 
% finding their centroids, measuring their
% intensities, and measuring the red intensity
% in the same region

[centroids, finalgreen, spot_projection] = find_centroids(numplanes, imLOG_spot, imgauss_spot, imgauss_SAC, spot_projection, finalgreen, dir);

'Found the spots in the green image and measured the intensities in the green and red channels.'

%%

% Link spots visible in adjacent z-planes.
% This part of the script will eliminate 
% any spot that is not present in at
% least two adjacent planes.

[centroids_unique, centroids_links, indexes, indexes_unique] = link_centroids(numplanes, centroids);

'Stitched spots together across z-stacks.'

%%

% Find the planes at which each spot is at its
% maximum intensity. Creates a list of non-
% redundant spots at their max intensities. 

spotmaxes = find_peak_intensities(centroids_unique, centroids_links, indexes, indexes_unique);

'Found the peak positions of the spots.'

figure;
imshow(MIP_spot,[]);
hold on, plot(spotmaxes(:,7) + offsetGreenYremainder, spotmaxes(:,8) + offsetGreenXremainder, 'b+', 'LineWidth', .25, 'MarkerSize',2); 
hold off
% Prevent it from making the background white
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, strcat(dir, 'spots_MIP.tif'));
close


%%

'Now plotting the spots for your review and saving the results.'

% Create a place to keep track of the final results
final_results = [];
finalcount = 0;
pixelsize_nm = 64.5;

numSpots = size(spotmaxes);

% Loop over the spots
for i=1:numSpots(1)
    
    % Determine which plane the spot was in
    plane = spotmaxes(i,2);
    
    % Go through the spots one by one to validate them
    
    % Recalculate the coords to correct the remaining red:green offset
    offsetspotcoords = [spotmaxes(i,8) + offsetGreenYremainder, spotmaxes(i,7) + offsetGreenXremainder];
    
    % Set figure display size based on screen size
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 scrsz(4) scrsz(3)*(2/3) scrsz(4)]);
    
    % Set the width of the image you want to see
    % Actual width will be 2 * dimensions
    dimensions = 48;
    
    % Get the coords of the spot and round it
    % for indexing, then determine the adjustment factor
    % for correcting the coords of the localizations for
    % displaying.
    midy = spotmaxes(i,8);
    midy_floor = floor(midy);
    adj_y = midy_floor - dimensions - 1;
    midx = spotmaxes(i,7);
    midx_floor = floor(midx);
    adj_x = midx_floor - dimensions - 1;
    
    % Calculate the display range for the nucleus
    horiz_range = midx_floor-dimensions:midx_floor+dimensions;
    vert_range = midy_floor-dimensions:midy_floor+dimensions;
    
    % Display all the spots from the same nucleus
    one = subplot(2,2,1);
    oneposition = get(one, 'pos');
    oneposition(1) = oneposition(1) + 0.05;
    oneposition(2) = oneposition(2) - 0.05;
    oneposition(3:4) = oneposition(3:4) + 0.05;
    set (one, 'pos', oneposition);
    hold on
    MIP_SACOfInterest = MIP_SAC(vert_range,horiz_range);
    MIP_spotOfInterest = MIP_spot(vert_range,horiz_range);
    edgetest = MIP_spotOfInterest;
    MIP_spotOfInterest = MIP_spotOfInterest - 12000;
    MIP_spotOfInterest(find(MIP_spotOfInterest < 0)) = 0;
    colorim = cat(3,MIP_SACOfInterest,MIP_spotOfInterest,MIP_SACOfInterest);
    imshow(colorim);
    hold on, plot(spotmaxes(i,7)-adj_x, spotmaxes(i,8)-adj_y,'b+', 'LineWidth', 2, 'MarkerSize',6);
    hold off
    
    % Display the images along with the localization data
    two = subplot(2,2,2);
    twoposition = get(two, 'pos');
    twoposition(1) = twoposition(1) - 0.05;
    twoposition(2) = twoposition(2) - 0.05;
    twoposition(3:4) = twoposition(3:4) + 0.05;
    set (two, 'pos', twoposition);
    hold on
    redgain = 2;
    MIP_SACOfInterest = (MIP_SAC(vert_range,horiz_range)).*redgain;
    greengain = 2;
    MIP_spotOfInterest = (MIP_spot(vert_range,horiz_range)).*greengain;
    colorim = cat(3,MIP_SACOfInterest,MIP_spotOfInterest,MIP_SACOfInterest);
    imshow(colorim);
    hold on, plot(spotmaxes(i,7)-adj_x, spotmaxes(i,8)-adj_y,'b+', 'LineWidth', 2, 'MarkerSize',6);
    hold off

    % Display the nucleus and spot along with the localization data
    three = subplot(2,2,3);
    threeposition = get(three, 'pos');
    threeposition(1) = threeposition(1) + 0.05;
    threeposition(2) = threeposition(2) + 0.05;
    threeposition(3:4) = threeposition(3:4) + 0.05;
    set (three, 'pos', threeposition);
    hold on
    
    % Extract the image of MAD1 and spots from the
    % reference image based on the calculations above
    
    % Sum the neighboring 2 planes to increase your 
    % ability to see a MAD1 peak
    try
        
        red = finalred(vert_range,horiz_range,plane-1) + finalblue(vert_range,horiz_range,plane) + finalblue(vert_range,horiz_range,plane+1);
    
    % But just use the single plane if there aren't
    % neighboring planes (at the top and bottom)
    catch error
        
        red = finalred(vert_range,horiz_range,plane);
        
    end
    
    green = finalgreen(vert_range,horiz_range,plane);
    
    % Display the color image
    colorim = cat(3,red*2,green*3,red*2);
    peakint = max(max(red));
    
    if peakint <= 0
        peakint = 1;
    end
    
    imshow(colorim, [0 peakint]);
    
    % Show the position of the spot as a blue plus
    hold on, plot(spotmaxes(i,7)-adj_x, spotmaxes(i,8)-adj_y,'b+', 'LineWidth', 2, 'MarkerSize',6);
    
    % Supervise spot inclusion
    % User should classify spots based on morphology into
    % the following categories. The final results will have a 
    % column where the number input here is stored. If you want
    % to add or eliminate categories, just modify the list below
    % to reflect your classification scheme.
    status = input('1: spindle\n2: MT\n3: CEN detached\n4: CEN lateral lateral\n5: CEN end-on\n6: reject\n');
    reject = 6; % Number of classifications
    
    % Eliminate the spot if the cell was on the
    % border of the image
    if (sum(edgetest(1,:)) * sum(edgetest(41,:)) * sum(edgetest(:,1)) * sum(edgetest(:,41))) == 0
        status = reject;
    end
    
    % If the user classified the spot, store the data
    if status < reject
        
        % Update the number of spots in the data
        finalcount = finalcount + 1;
        % Take the spot data
        finalspot = spotmaxes(i,:);
        % Update with the corrected coords
        finalspot(7:8) = [offsetspotcoords(2) offsetspotcoords(1)];
        % Replace the index of the spot within the plane
        % with the overall spot number
        finalspot(3) = finalcount;
        
        % Append the results
        final_results = [final_results
            stack finalspot status];
        
        % Save the figure
        % Prevent it from making the background white
        set(gcf, 'InvertHardCopy', 'off');
        saveas(gcf, strcat(dir, num2str(finalcount),'.tif'));
        hold off
        
    end
    
    close all

end

'Finished.'
