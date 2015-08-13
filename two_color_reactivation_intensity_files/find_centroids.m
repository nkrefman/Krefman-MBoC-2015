function [centroids, finalgreen, spot_projection] = find_centroids(numplanes, imLOG_spot, imgauss_spot, imgauss_SAC, spot_projection, finalgreen, dir)

%% Nathaniel I. Krefman 2014-05-06

% Loops through a stack identifying spots and finding their centroids

% Input variables:
%	numplanes:              number of planes in the image
%   imLOG_spot:             LOG filtered stack for filtering 
%   imgauss_spot:           Gaussian filtered stack
%   imgauss_SAC:            Gaussian filtered stack
%   spot_projection:        an empty padded image for storing 
%                           thresholded spots
%   finalgreen:             an empty padded stack for storing
%                           background-subtracted spots
%   dir:                    Directory for storing the results.
%                           The directory must already exist!

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   level
%   intensity_threshold
%   spread_threshold & area_threshold
%   bg_mask & spot_mask
%   im_max

% Other file requirements:
%       None

%%

% Set background threshold.
% Threshold is used to identify spots.
% Value can be set arbitrarily, since this threshold is not 
% used to quantify fluorescence intensities.
% Should be the same value for all images from an experiment.
level = 0.035; 

% Set spot intensity threshold
intensity_threshold = 0; % 0 includes all spots

for i=1:numplanes
    
    % Convert to binary
    filtered(:,:,i) = im2bw(imLOG_spot(:,:,i),level);
    % Fill in holes in regions
    filtered(:,:,i) = imfill(filtered(:,:,i),'holes');
    % Eliminate all regions smaller than 20 pixels
    filtered(:,:,i) = bwareaopen(filtered(:,:,i),25);
    % Make a connectivity map of the regions
    regions = bwconncomp(filtered(:,:,i),8);
    
    % Calc the coords of the boxes containing the regions
    borders = regionprops(regions, 'BoundingBox');
    boxes = cat(1, borders.BoundingBox);
        
    % Filter ROIs whose peak intensity is below 
    % a threshold.
    
    numregs = size(boxes);  % Calculate number of regions
    ROI_idx = [];           % To store indices of ROIs
        
    % Loop through the regions to find spots
    for j=1:numregs(1)
        
        try
            
            % Extract coords of the region to plot the spot.
            % BoundingBox gives coordinates as half pixels
            % so this corrects for that.box = boxes(j,:);
            box = boxes(j,:);
            xmin = box(1,1);
            xmin = xmin + 0.5;
            xmax = xmin+box(1,3)-1;
            ymin = box(1,2);
            ymin = ymin + 0.5;
            ymax = ymin+box(1,4)-1;
            
            % Get spot image from the reference image
            spot_ref = imgauss_spot(ymin:ymax,xmin:xmax,i);
            
            % Calculate the max pixel intensity of the region
            max_intensity = max(max(spot_ref));
            
            % Store indices/intensities of ROIs above threshold
            if max_intensity >= intensity_threshold
                
                ROI_idx = [ROI_idx
                    j];
                
            end
            
%         % Uncomment to see the error message if this fails.
%         catch err
%             
%             rethrow(err)
            
        end
        
    end;
        
    % Create a map of the ROIs above intensity threshold
    intensity_thresh = ismember(labelmatrix(regions), ROI_idx);
    regions = bwconncomp(intensity_thresh,8);
    
    % Calculate region properties to eliminate
    % irregular spots (high perimeter:area ratio).
    
    perim_area = regionprops(regions, 'Perimeter', 'Area');
    perimeter = cat(1, perim_area.Perimeter);
    area = cat(1, perim_area.Area);
    
    % Takes ratio of perimeter/area, while correcting
    % for how circumference scales with the square root
    % of the area.  Otherwise, this eliminates small spots
    % that are still quite round.
    spread = (perimeter ./ area) .* ((area) ./ (2*sqrt(area)));
    
    % Set a threshold for spread (2.8 is good to start. 100 allows identification of MTs)
    % and area (16-25 is a good range for finding only spots)
    spread_threshold = 1000; % Early on, we were much more stringent (2.8), but later improvements to the code eliminated spurious spots through other means, and so we increased the threshold to 1000.
    area_threshold = 5; % Spot must be at least this number of pixels.
    
    % Create a map of the ROIs below spread threshold
    shape_idx = find((spread < spread_threshold) & (area > area_threshold));
    spots = ismember(labelmatrix(regions), shape_idx);

    % Add the good spots to the spots projection
    % Ultimately, we'll only use spots that were seen
    % in two or more adjacent planes.
    spot_projection = spot_projection + spots;
    
    % Make a connectivity map of the spots
    ROIs = bwconncomp(spots,8);
    ROImap(i) = ROIs;

    % Calc the coords of the boxes containing the regions
    boundaries = regionprops(ROIs, 'BoundingBox');
    containers = cat(1, boundaries.BoundingBox);
    % Calculate the intensities of spots
    % and their coordinates.
    
    numspots = size(containers);    % Calculate number of regions
    ROI_idx = [];                   % To store indices of ROIs
    ROI_data = [];                  % To store intensities
    
    % Loop through the ROIs to find spots for
    % localization and measuring integrated intensity.
        
    for j=1:numspots(1)
        
        % Extract coords of the region to plot the spot.
        % BoundingBox gives coordinates as half pixels
        % so this corrects for that.
        container = containers(j,:);
        xmin = container(1,1);
        xmin = xmin + 0.5;
        xwidth = container(1,3);
        xmax = xmin+xwidth-1;
        ymin = container(1,2);
        ymin = ymin + 0.5;
        ywidth = container(1,4);
        ymax = ymin+ywidth-1;

        % Measure the spot intensities
        
        % Calculate local background using local pixels below
        % the threshold set previously
        bg_ref = imgauss_spot(ymin-1:ymax+1,xmin-1:xmax+1,i);
        bg_mask = ~filtered(ymin-1:ymax+1,xmin-1:xmax+1,i);      
        filtered_bg = bg_ref.*bg_mask;
        bg_sum = sum(sum(filtered_bg));
        num_bg_pxls = sum(sum(bg_mask));
        bg_per_pxl = bg_sum/num_bg_pxls;
                
        % Calculate a refined intensity by subtracting
        % the local background
        spot_ref = imgauss_spot(ymin:ymax,xmin:xmax,i);
        spot_mask = filtered(ymin:ymax,xmin:xmax,i);
        filtered_spot = spot_ref.*spot_mask;
        spot_sum = sum(sum(filtered_spot));
        num_spot_pxls = sum(sum(spot_mask));
        spot_bg = bg_per_pxl * num_spot_pxls;
        integrated_intensity_spot = spot_sum - spot_bg;
        
        % Generate a background-subtracted image for
        % localization.
        bg_matrix = (spot_mask .* 0) + bg_per_pxl;
        subtracted_spot = (spot_ref - bg_per_pxl).*spot_mask;

        % Add the subtracted spot image to the final green image
        finalgreen(ymin:ymax,xmin:xmax,i) = finalgreen(ymin:ymax,xmin:xmax,i) + subtracted_spot;
        
        % Localization by peak intensity pixel.
        [Ycom, Xcom]=find(subtracted_spot == max(max(subtracted_spot)));
        
        % Correct for fitting within an ROI
        Xcoord = Xcom + xmin-1;
        Ycoord = Ycom + ymin-1;
        
        % If more than one pixel has the max value
        % choose one randomly
        numMaxPix = size(Xcoord);
        if numMaxPix(1) > 1
            choice = ceil(rand(1)*numMaxPix(1));
            Xcoord = Xcoord(choice);
            Ycoord = Ycoord(choice);
        end
        
        % Recalculate signal and background
        
        % Define new spot & background masks
        bg_mask = [0	0	1	1	1	1	1	1	1	0	0
            0	1	1	0	0	0	0	0	1	1	0
            1	1	0	0	0	0	0	0	0	1	1
            1	0	0	0	0	0	0	0	0	0	1
            1	0	0	0	0	0	0	0	0	0	1
            1	0	0	0	0	0	0	0	0	0	1
            1	0	0	0	0	0	0	0	0	0	1
            1	0	0	0	0	0	0	0	0	0	1
            1	1	0	0	0	0	0	0	0	1	1
            0	1	1	0	0	0	0	0	1	1	0
            0	0	1	1	1	1	1	1	1	0	0];
        num_bg_pxls = sum(sum(bg_mask));
        spot_mask = [0	0	0	0	0	0	0	0	0	0	0
            0	0	0	1	1	1	1	1	0	0	0
            0	0	1	1	1	1	1	1	1	0	0
            0	1	1	1	1	1	1	1	1	1	0
            0	1	1	1	1	1	1	1	1	1	0
            0	1	1	1	1	1	1	1	1	1	0
            0	1	1	1	1	1	1	1	1	1	0
            0	1	1	1	1	1	1	1	1	1	0
            0	0	1	1	1	1	1	1	1	0	0
            0	0	0	1	1	1	1	1	0	0	0
            0	0	0	0	0	0	0	0	0	0	0];
        num_spot_pxls = sum(sum(spot_mask));
        
        % Calculate local background using
        % the new background mask
        bg_ref = imgauss_spot(Ycoord-5:Ycoord+5,Xcoord-5:Xcoord+5,i);
        filtered_bg = bg_ref.*bg_mask;
        bg_sum = sum(sum(filtered_bg));
        bg_per_pxl = bg_sum/num_bg_pxls;
        
        % Calculate a refined intensity by subtracting
        % the local background
        spot_ref = imgauss_spot(Ycoord-5:Ycoord+5,Xcoord-5:Xcoord+5,i);
        filtered_spot = spot_ref.*spot_mask;
        spot_sum = sum(sum(filtered_spot));
        spot_bg = bg_per_pxl * num_spot_pxls;
        integrated_intensity_spot = spot_sum - spot_bg;
        
        % Regenerate a background-subtracted image
        % for localization.
        bg_matrix = (spot_mask .* 0) + bg_per_pxl;
        subtracted_spot = (spot_ref - bg_per_pxl).*spot_mask;
        
        % Measure the SAC intensities in the same place
        
        % Calculate local background using local pixels below
        % the threshold set previously
        bg_ref = imgauss_SAC(Ycoord-5:Ycoord+5,Xcoord-5:Xcoord+5,i);
        filtered_bg = bg_ref.*bg_mask;
        bg_sum = sum(sum(filtered_bg));
        bg_per_pxl = bg_sum/num_bg_pxls;
                
        % Calculate a refined intensity by subtracting
        % the local background
        spot_ref = imgauss_SAC(Ycoord-5:Ycoord+5,Xcoord-5:Xcoord+5,i);;
        filtered_spot = spot_ref.*spot_mask;
        spot_sum = sum(sum(filtered_spot));
        spot_bg = bg_per_pxl * num_spot_pxls;
        integrated_intensity_SAC = spot_sum - spot_bg;
        
%         % Uncomment to see spot, mask, & localizations
%         if i < numImages % Can ignore by setting above numImages
%             range = [-40 150] % optional display range
%             hold off
%             close all
%             imshow(subtracted_spot, [], 'InitialMagnification', 2000);
%             hold on, plot(6, 6, 'b+')
%             hold off
%             wait = input('Hit enter/return to continue.');
%             saveas(gcf, strcat(dir, 'COM.tif'));
%             close
%         end
        
        ROI_idx = [ROI_idx
            j];
        
        ROI_data = [ROI_data
            j integrated_intensity_spot integrated_intensity_SAC spot_sum Xcoord Ycoord];
        
    end;
    
    % Append the spot data for this plane to the spot
    % data for all the planes.
    numSpots = size(ROI_data);
    centroids(1:numSpots(1),1:numSpots(2),i) = ROI_data;

end

%%

% Set everything below zero to zero
finalgreen = double(finalgreen);
finalgreen(find(finalgreen < 0)) = 0;

%%

% Save a sum-projection image of all
% the background-subtracted spots
sumspots = zeros(size(finalgreen(:,:,1)));
for i=1:numplanes
    sumspots = sumspots + finalgreen(:,:,i);
end
maxint = max(max(sumspots));
imshow(sumspots, [0 maxint]);
colormap('Jet')
colorbar
saveas(gcf, strcat(dir, 'spot_sum.tif'));
close

% Save a projection of all the spots
close all
figure;
hold on
imshow(spot_projection,[]);
hold off
saveas(gcf, strcat(dir, 'spot_projection.tif'));
close

%%

% Convert to a uint16 image
% with a fixed display range

% make sure the image is a double
finalgreen = double(finalgreen);

% select a max instensity (im_max)
im_max = 170;

% adjust the display range for uint16 
% (i.e. 65535 values), setting any 
% pixel over im_max to saturated
finalgreen = 65535*finalgreen/im_max;

% turn the image into a uint16 array
finalgreen = uint16(finalgreen);