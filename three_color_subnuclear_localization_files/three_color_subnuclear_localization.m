function final_results = map_nuclear_spots_3colors(dir, stack, imG, imR, imB)
%% Nathaniel I. Krefman 2014-05-06

% Summary:

% Script for determining the subnuclear localization of spots
% using a three-color imaging set up. Usually spots are in green
% and Nsg1-mCherry is the a marker for the nuclear boundary.
% Spc42-CFP is a reference for spindle poles to determine if
% kinetochores are attached or detached.

% Adjust the min/max range of the images using adjRangeAndBitDepth.m 
% or ImageJ before running.

% We generally call this script using run_3csl.m, although it 
% can also be called directly from the Command Window. The 
% variables in run_3csl.m need to be updated every time you run 
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
%   imG, imR, & imB:  green (Mtw1-3xGFP), red (MAD1-3xmCherry), & 
%               blue (Spc42-CFP) images.
%               Give the complete paths.

% Dependencies (in order of first call):
%	open_spots.m
%	openTiffStack.m
%	autogain_uint16.m
%   find_centroids.m
% 	link_centroids.m
%	find_peak_intensities.m
%	fit_nuclei.m
%	medthresh.m
% 	match_spots_to_nuclei.m
%  	group_spots_within_cells.m
%  	autogain_uint16.m

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   numplanes
%   offsetGreenYprecise & offsetGreenXprecise
%   offsetBlueYprecise  &offsetBlueXprecise
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
% Here, we shift the green and blue image2 to match the red.
% NOTE: This will require some playing to be sure you have 
% the signs correct, and you don't have Y and X backwards.
offsetGreenYprecise = -4.3470;
offsetGreenXprecise = -0.5360;
offsetBlueYprecise = -2.0000;
offsetBlueXprecise = -1.0000;

% Calculate the correction factor for displaying
% the corrected colocalization without interpolation
offsetGreenYdisplay = round(offsetGreenYprecise);
offsetGreenXdisplay = round(offsetGreenXprecise);
offsetBlueYdisplay = round(offsetBlueYprecise);
offsetBlueXdisplay = round(offsetBlueXprecise);

% Calculate the remaining subpixel offset.
% We use these values to correct our green:red distance 
% measurements much later!
offsetGreenYremainder = offsetGreenYprecise - offsetGreenYdisplay;
offsetGreenXremainder = offsetGreenXprecise - offsetGreenXdisplay;

% We don't calculate the offset remainder for blue, because we 
% don't need that information for our purpose.

%%

% Open the green image, create a
% max-intensity projection, Gaussian-
% and LOG-filtered images, a blank 
% image to store a spot projection,
% and a contrast-adjusted image for 
% supervised analysis.

[MIP_spot, spot_projection, imLOG_spot, imgauss_spot, finalgreen] = open_spots(imG, numplanes, offsetGreenYdisplay, offsetGreenXdisplay);

'Opened the image and filtered it.'

%%

% Loop through stack identifying spots, 
% and finding their centroids

[centroids, finalgreen, spot_projection] = find_centroids(numplanes, imLOG_spot, imgauss_spot, spot_projection, finalgreen, dir);

'Found the centroids in the green image.'

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

%%

% Open the nucleus image, filter it, then 
% calculate the centroids and approximate 
% radii of the nuclei.

[finalred, maxCircles, centers, radii] = fit_nuclei(imR, numplanes);

'Fit the nuclei and calculated their centroids.'

%%

% Identify the nuclei to which each spot belongs.

spots_nuclei = match_spots_to_nuclei(numplanes, spotmaxes, maxCircles, centers, radii);

'Matched spots to their corresponding nuclei.'

%%

% Figure out which spots were in the same nucleus
% by finding out how close the nuclei are to each 
% other.

% Check to see that the spots have been grouped 
% appropriately

grouped_spots_nuclei = group_spots_within_cells(spots_nuclei)

'Matched spots to their corresponding nuclei.'

%%

% Uncomment to see the localizations of the spots
% at their max intensity planes superimposed over
% the maximum intensity projection.

figure;
imshow(MIP_spot,[]);
hold on, plot(grouped_spots_nuclei(:,6), grouped_spots_nuclei(:,7), 'b+', 'LineWidth', .25, 'MarkerSize',2); 
hold off
% wait = input('Hit enter/return to continue.'); % Uncomment to hold.
% Prevent it from making the background white
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, strcat(dir, 'spots_MIP.tif'));
close

figure;
imshow(MIP_spot,[], 'InitialMagnification', 500);
hold on, plot(grouped_spots_nuclei(:,6), grouped_spots_nuclei(:,7), 'b+', 'LineWidth', .25, 'MarkerSize', 2); 
hold on, plot(grouped_spots_nuclei(:,8), grouped_spots_nuclei(:,9),'ro', 'LineWidth', .25, 'MarkerSize', 2);
h = viscircles(grouped_spots_nuclei(:,8:9),grouped_spots_nuclei(:,10), 'LineWidth', .25, 'EdgeColor', 'w');
hold off
% wait = input('Hit enter/return to continue.'); % Uncomment to hold.
% Prevent it from making the background white
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, strcat(dir, 'spots_nuclei_MIP.tif'));
close

%%

% Open the blue image, create a
% max-intensity projection, Gaussian-
% -filtered images, and a contrast-
% adjusted image for supervised
% analysis.

[MIP_pole, imgauss_pole, finalblue] = open_poles(imB, offsetBlueYdisplay, offsetBlueXdisplay);

'Opened the image and filtered it.'

%%

figure;
imshow(MIP_pole,[], 'InitialMagnification', 500);
hold on, plot(grouped_spots_nuclei(:,6), grouped_spots_nuclei(:,7), 'b+', 'LineWidth', .25, 'MarkerSize', 2); 
hold on, plot(grouped_spots_nuclei(:,8), grouped_spots_nuclei(:,9),'ro', 'LineWidth', .25, 'MarkerSize', 2);
hold off
% wait = input('Hit enter/return to continue.'); % Uncomment to hold.
% Prevent it from making the background white
set(gcf, 'InvertHardCopy', 'off');
saveas(gcf, strcat(dir, 'poles_nuclei_MIP.tif'));
close

%%

% Create a place to keep track of the final results
final_results = [];
finalcount = 0;
pixelsize_nm = 64.5;

numSpots = size(grouped_spots_nuclei);

% Loop over the spots
for i=1:numSpots(1)
    
    % Get the xcoords and ycoords from the spot and nucleus
    xs = [grouped_spots_nuclei(i,6),grouped_spots_nuclei(i,8)];
    ys = [grouped_spots_nuclei(i,7),grouped_spots_nuclei(i,9)];
    
    % Fit a line to the points (y = mx + b)
    line = polyfit(xs,ys,1);
    slope = line(1);
    intercept = line(2);

    % Specify the length of the line to trace
    length = 20;
    % Determine how many x pixels the line must cross
    % to make the a line of the specified length
    % based on the slope of the line
    window = length/sqrt(1+(slope*slope));
    % Specify the number of points to compute on the line
    sample = 100;
    % Specify the x spacing for the desired number of points
    increment = window/sample;
    % Get the x coord of the nuclear center
    nucleusx = grouped_spots_nuclei(i,8);
    
    % Default is to assume the spot is to the left of the
    % the center of the nucleus (true 50% of the time)
    right = 0;
    
    % If the spot is actually to the right
    if grouped_spots_nuclei(i,6) < grouped_spots_nuclei(i,8)
        
        % Flag it
        right = 1;
        
    end
    
    % If the spot was on the left
    if right == 0
        
        % Sample from the center of the nucleus
        % to the right
        xrange = nucleusx:increment:nucleusx+window;

    else
        
        % Sample from the center of the nucleus 
        % to the left
        xrange = nucleusx-window:increment:nucleusx;

    end
    
    % Evaluate the equation with the specified range
    yfit = slope*xrange + intercept;
    
    % Determine which plane the spot was in
    plane = grouped_spots_nuclei(i,3);
    
    % Measure the intensities along the line
    prof_nuc = improfile(finalred(:,:,plane),xrange,yfit);
    prof_spot = improfile(imgauss_spot(:,:,plane),xrange,yfit); %xxx
        
    % Find the coordinate of maximum intensity.
    boundarypeak = find(prof_nuc == max(prof_nuc));
    boundarycoordX = mean(xrange(boundarypeak));
    boundarycoordY = slope*boundarycoordX + intercept;
    % Adjust the line profiles so all the points are above zero
    prof_nuc = prof_nuc - min(prof_nuc);
    prof_spot = prof_spot - min(prof_spot);
    % Adjust the line profile of the spots data so that it has
    % the same maximum as the line profile of the nucleus data
    ratio = max(prof_spot)/max(prof_nuc);
    prof_spot = prof_spot ./ ratio;
        
    % Go through the spots one by one to validate them
    
    % Recalculate the coords to correct the remaining red:green offset
    offsetspotcoords = [grouped_spots_nuclei(i,7) + offsetGreenYremainder, grouped_spots_nuclei(i,6) + offsetGreenXremainder];
    % Get the coords for the fitted nuclear boundary
    boundarycoords = [boundarycoordY, boundarycoordX];
    centercoords = [grouped_spots_nuclei(i,9), grouped_spots_nuclei(i,8)];
    % Calculate the distance between the spot and
    % the nuclear boundary
    boundarydistance = pdist2(offsetspotcoords,boundarycoords,'euclidean');
    % Calculate the distance between the spot and
    % the center of the nucleus
    centerdistance = pdist2(offsetspotcoords,centercoords,'euclidean');
    % Calculate the nuclear radius
    radius = pdist2(boundarycoords,centercoords,'euclidean');
    
    % Filter out all nuclei with aberrantly small radii
    % or that peak at the last pixel of the line profile
    
    if (radius > 3) & (radius <= 19)
    
        % Convert distances to nm
        boundarydistance_nm = boundarydistance * pixelsize_nm;
        centerdistance_nm = centerdistance * pixelsize_nm;
        radius_nm = radius * pixelsize_nm;
        
        % Determine if the spot is in zone 1, 2, or 3
        % If the spot is farther from the center than the fitted nuclear
        % boundary (correcting for the final offset, which is necessarily
        % less than 1 pixel in either axis (or the square root of 2
        % pixels in both axes)
        if (boundarydistance + centerdistance) > (radius + sqrt(2))
            
            zone = 1;
            % The spot is on the outside of the boundary
            % so switch its sign to negative
            boundarydistance_nm = -boundarydistance_nm;
            
        elseif boundarydistance < (0.184*radius)
            
            zone = 1;
            
        elseif boundarydistance < (0.422*radius)
            
            zone = 2;
            
        elseif boundarydistance >= (0.422*radius)
            
            zone = 3;
            
        end
        
        close all
        
        scrsz = get(0,'ScreenSize');
        figure('Position',[1 scrsz(4) scrsz(3)*(2/3) scrsz(4)]);        
        % Display the line profiles and Gaussian fit
        subplot(2,2,1);
        
        % Set the y-axis to adjust the range range
        axismax = max(prof_nuc)*1.5;
        
        if axismax <= 0
            axismax = 1;
        end
        
        hold on, axis([min(xrange) max(xrange) 0 axismax]);
        % Make the plot background black
        set(subplot(2,2,1),'Color','Black')
        % Show the position of the nuclear center as a white circle
        hold on, plot(grouped_spots_nuclei(i,8), 0, 'wo', 'LineWidth', 2, 'MarkerSize',12);
        % Plot the spot profile along the line as a green line
        hold on, plot(xrange,prof_spot,'g', 'LineWidth', 2);
        % Show the position of the spot as a blue x
        hold on, plot(grouped_spots_nuclei(i,6),0,'bx', 'LineWidth', 2, 'MarkerSize',12);
        % Plot the nuclear profile along the line as a magenta line
        hold on, plot(xrange,prof_nuc,'m', 'LineWidth', 2);
        % Show the position of the nuclear edge as a red x
        hold on, plot(boundarycoordX,0,'rx', 'LineWidth', 2, 'MarkerSize',12);
        hold off
        
        % Set the width of the image you want to see
        % Width will be 2 * dimensions
        dimensions = 20;
        
        % Get the coords of the nuclear center and round it
        % for indexing, then determine the adjustment factor
        % for correcting the coords of the localizations for
        % displaying.
        midy = grouped_spots_nuclei(i,9);
        midy_floor = floor(midy);
        adj_y = midy_floor - dimensions - 1;
        midx = grouped_spots_nuclei(i,8);
        midx_floor = floor(midx);
        adj_x = midx_floor - dimensions - 1;
        
        % Calculate the display range for the nucleus
        horiz_range = midx_floor-dimensions:midx_floor+dimensions;
        vert_range = midy_floor-dimensions:midy_floor+dimensions;
        
        % Display all the spots from the same nucleus
        subplot(2,2,2);
        hold on 
        greengain = 1.2;
        MIP_cellspots = MIP_spot(vert_range,horiz_range).*greengain;
        bluegain = 3;
        MIP_cellpoles = MIP_pole(vert_range,horiz_range).*bluegain;
        colorim = cat(3,MIP_cellpoles,MIP_cellspots,MIP_cellpoles);            
        imshow(colorim);
        samenucleus = find(grouped_spots_nuclei(:,1) == grouped_spots_nuclei(i,1));
        hold on, plot(grouped_spots_nuclei(samenucleus,6)-adj_x, grouped_spots_nuclei(samenucleus,7)-adj_y,'b+', 'LineWidth', 1, 'MarkerSize',4);
        hold on, plot(grouped_spots_nuclei(samenucleus,8)-adj_x, grouped_spots_nuclei(samenucleus,9)-adj_y,'ro', 'LineWidth', 1, 'MarkerSize',4);
        h = viscircles([grouped_spots_nuclei(samenucleus,8)-adj_x,grouped_spots_nuclei(samenucleus,9)-adj_y],grouped_spots_nuclei(samenucleus,10), 'LineWidth', .25, 'EdgeColor', 'w');
        hold off

        % Display the spot and the reference image (e.g. spindle pole 
        % marker) along with the localization data
        subplot(2,2,3);
        hold on
        
        % Extract the image of the reference and spots from the
        % reference image based on the calculations above
        
        try
       
            blue = finalblue(vert_range,horiz_range,plane-1) + finalblue(vert_range,horiz_range,plane) + finalblue(vert_range,horiz_range,plane+1);
        
        catch error
            
            blue = finalblue(vert_range,horiz_range,plane);
        
        end
        
        green = finalgreen(vert_range,horiz_range,plane);
        
        % Display the color image
        colorim = cat(3,blue,green,blue);            
        imshow(colorim);
     
        % Show the position of the spot as a blue +
        hold on, plot(grouped_spots_nuclei(i,6)-adj_x, grouped_spots_nuclei(i,7)-adj_y,'b+', 'LineWidth', 2, 'MarkerSize',8);
        % Show the position of the nuclear center as a white circle
        hold on, plot(grouped_spots_nuclei(i,8)-adj_x, grouped_spots_nuclei(i,9)-adj_y,'wo', 'LineWidth', 2, 'MarkerSize',8);
        % Show the axis of the line scan as a white line
        hold on, plot(xrange-adj_x,yfit-adj_y, 'w');
        % Show the position of the Gaussian-fit nuclear edge as a red +
        hold on, plot(boundarycoordX-adj_x, boundarycoordY-adj_y, 'r+', 'LineWidth', 2, 'MarkerSize',8);
        
        % Show the fit of imfindcirces as a white circle
        % NOTE: Spots might appear outside of their real zone
        % by a bit, because the image does not display the final
        % correction for the color offsets on the microscope.
        
        try
            
            h = viscircles([grouped_spots_nuclei(i,8)-adj_x,grouped_spots_nuclei(i,9)-adj_y],radius,'LineWidth',1,'EdgeColor', 'w');
            h = viscircles([grouped_spots_nuclei(i,8)-adj_x,grouped_spots_nuclei(i,9)-adj_y],radius - (0.184*radius),'LineWidth',1,'EdgeColor', 'w');
            h = viscircles([grouped_spots_nuclei(i,8)-adj_x,grouped_spots_nuclei(i,9)-adj_y],radius - (0.422*radius),'LineWidth',1,'EdgeColor', 'w');
        
        catch error
            
            'failed'
                
        end
        
        % Display the nucleus and spot along with the localization data
        % and the line along which the profile is taken
        subplot(2,2,4);
        hold on
        
        % Extract the image of the nucleus and spots from the
        % reference image based on the calculations above

        try

            red = finalred(vert_range,horiz_range,plane-1) + finalred(vert_range,horiz_range,plane) + finalred(vert_range,horiz_range,plane+1);
        
        catch error    
            
            red = finalred(vert_range,horiz_range,plane);
            
        end
        
        green = finalgreen(vert_range,horiz_range,plane);
        
        % Display the color image
        colorim = cat(3,red,green,red);
        peakint = max(prof_nuc);
        
        if peakint <= 0
            peakint = 1
        end
            
        imshow(colorim, [0 peakint]);
     
        % Show the position of the spot as a blue +
        hold on, plot(grouped_spots_nuclei(i,6)-adj_x, grouped_spots_nuclei(i,7)-adj_y,'b+', 'LineWidth', 2, 'MarkerSize',8);
        % Show the position of the nuclear center as a white circle
        hold on, plot(grouped_spots_nuclei(i,8)-adj_x, grouped_spots_nuclei(i,9)-adj_y,'wo', 'LineWidth', 2, 'MarkerSize',8);
        % Show the axis of the line scan as a white line
        hold on, plot(xrange-adj_x,yfit-adj_y, 'w');
        % Show the position of the Gaussian-fit nuclear edge as a red +
        hold on, plot(boundarycoordX-adj_x, boundarycoordY-adj_y, 'r+', 'LineWidth', 2, 'MarkerSize',8);
        
        % Show the fit of imfindcircles as a white circle
        % NOTE: Spots might appear outside of their real zone
        % by a bit, because the image does not display the final
        % correction for the color offsets on the microscope.
        
        try
            
            h = viscircles([grouped_spots_nuclei(i,8)-adj_x,grouped_spots_nuclei(i,9)-adj_y],radius,'LineWidth',1,'EdgeColor', 'w');
            h = viscircles([grouped_spots_nuclei(i,8)-adj_x,grouped_spots_nuclei(i,9)-adj_y],radius - (0.184*radius),'LineWidth',1,'EdgeColor', 'w');
            h = viscircles([grouped_spots_nuclei(i,8)-adj_x,grouped_spots_nuclei(i,9)-adj_y],radius - (0.422*radius),'LineWidth',1,'EdgeColor', 'w');
        
        catch error
            
            'failed'
                
        end
        
        % Supervise spot inclusion
        % User should classify spots based on morphology into
        % the following categories. The final results will have a
        % column where the number input here is stored. If you want
        % to add or eliminate categories, just modify the list below
        % to reflect your classification scheme.
        status = input('Type 1 for attached, 2 for detached, 3 for ambiguous, 4 for reject\n');
        reject = 6; % Number of classifications

        % Eliminate the spot if the cell was on the 
        % border of the image
        if (sum(MIP_cellspots(1,:)) * sum(MIP_cellspots(41,:)) * sum(MIP_cellspots(:,1)) * sum(MIP_cellspots(:,41))) == 0
        	status = reject;
        end
        
        % If the user says yes, store the data
        if status ~= reject
            
            % Update the number of spots in the data
            finalcount = finalcount + 1;
            
            % Group the spot data with the edgecoords and distance
            finalspot = [grouped_spots_nuclei(i,:) boundarycoordX boundarycoordY boundarydistance_nm radius_nm zone];
            % Update the corrected spot coords in the data
            finalspot(6:7) = [offsetspotcoords(2) offsetspotcoords(1)];
            % Replace the index of the spot within the plane
            % with the overall spot number
            finalspot(4) = finalcount;
                        
            % Append the results
            final_results = [final_results
                stack finalspot status];
            
            % Save the figure
            % Prevent it from making the background white
            set(gcf, 'InvertHardCopy', 'off');
            saveas(gcf, strcat(dir, num2str(finalcount),'.tif'));
            hold off
            
        end
    
    end
    
    close all
    
end

'Finished.'