function final_results = two_color_subnuclear_localization(dir, stack, imG, imR) % Could add filenames as input arguments

%% Nathaniel I. Krefman 2014-05-05

% Summary:

% Script for determining the subnuclear localization of spots
% using a two-color imaging set up. Usually spots are in green
% and Nsg1-mCherry is the a marker for the nuclear boundary.

% Adjust the min/max range of the images using ImageJ or
% adjRangeAndBitDepth.m before running.

% We generally call this script using run_tcsn.m, although it 
% can also be called directly from the Command Window. The 
% variables in run_tcsn.m need to be updated every time you run 
% the code for a new image stack.

% Most of the code in this program does not need to be updated 
% at runtime, except to optimize for localization of a new spot 
% type, or when using new microscope or illumination settings, 
% for example.

% We suggest adjusting the following variables for your set-up:
%   numplanes
%   offsetYprecise & offsetXprecise
%   kernelsize_spot & sigma_spot
%   level
%   intensity_threshold
%   spread_threshold & area_threshold
%   distance_threshold (NOTE: used at least three different ways -- for spot-spot distances, for spot-nuclei distances, and for nuclei-nuclei distances)
%   dead_cell_threshold
%   kernelsize_nuc & sigma_nuc
%   gausskernelsize_nuc & gausssigma_nuc
%   radius range and sensitivity of imfindcircles()

% Other parameters may require adjustment for different microscopes. 
% For instance, if the pixel size deviates substantially from 
% 65 nm, several distance parameters will need alteration.

% This script saves several output images, including data from every 
% nucleus. Be sure you have plenty of hard-drive space before 
% running. We delete the output images after we double check our 
% classifications due to digital storage constraints.

% Dependencies:
%   openTiffStack.m
%   autogain_uint16.m

%% Open spot image and apply LOG and Gaussian filters

hold off
close all

% Input the number of z-planes in your stack.
numplanes = 27;

% Open the spots image.
image_spot = openTiffStack(imG); % Need openTiffStack.m in Matlab path

% Convert the image to 8-bit.
image_spot = uint8(image_spot);

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
offsetYremainder = (offsetYprecise) - offsetYdisplay;
offsetXremainder = (offsetXprecise) - offsetXdisplay;

% Create a blank padded image.
pad = 50; % Pad width is arbitrary, but 50 is good.
[sizeY, sizeX, sizeZ] = size(image_spot);
paddedY = sizeY+pad*2;
paddedX = sizeX+pad*2;
padded = zeros(paddedY,paddedX);
padded3D = zeros(paddedY,paddedX,sizeZ);

% Create a maximum intensity projection (MIP) of the spot image.
padded(pad+1+offsetYdisplay:end-pad+offsetYdisplay, pad+1+offsetXdisplay:end-pad+offsetXdisplay) = max(double(image_spot),[],3);
MIP_spot = padded;

% Define the Laplacian of Gaussians (LOG) & Gaussian filters for spots.
kernelsize_spot = 6; % 6 finds peaks well, 5 is more conservative (misses dimmest spots)
sigma_spot = 4; % 4 finds peaks well
kernelgauss = fspecial('gaussian', kernelsize_spot, sigma_spot); % for finding spots
kernelspot = -fspecial('log', kernelsize_spot, sigma_spot); % for localizing spots & measuring intensities

% Filter then pad the stack using LOG & Gaussian filters.

% Make a LOG filtered image stack.
imtemp = imfilter(double(image_spot(:,:,:)), kernelspot, 'symmetric', 'conv');
padded3D(pad+1+offsetYdisplay:end-pad+offsetYdisplay, pad+1+offsetXdisplay:end-pad+offsetXdisplay,:) = imtemp(:,:,:);
clear imtemp
log_spot(:,:,:) = padded3D;

% Threshold filtered image using Otsu's Method.
% You should use the same threshold for all images. So try
% Otsu's on representative images, then use a consensus value.
level = graythresh(log_spot(:,:,:)) % Otsu's method, comment to set arbitrary threshold
% level = 0.04; % Uncomment to set arbitrary threshold

% Optional: Correct for failed thresholding
if level == 0 level = 0.04;

    'Thresholding this image failed.'

end

% Make a Gaussian filtered image stack.
imtempgauss = imfilter(double(image_spot(:,:,:)), kernelgauss, 'symmetric', 'conv');
padded3D(pad+1+offsetYdisplay:end-pad+offsetYdisplay, pad+1+offsetXdisplay:end-pad+offsetXdisplay,:) = imtempgauss(:,:,:);
imgauss_spot = padded3D;

% Set spot intensity threshold.
intensity_threshold = 0; % 0 includes all spots

clear padded3D imtempgauss image_spot

'Opened the image and filtered it.'

%% Localize spots & measure intensities

% Create a matrix to store the final background subtracted spot image.
finalspot = zeros(size(imgauss_spot));
% Create a matrix to store a projection of all the spots.
spot_projection = zeros(size(padded));

% Loop through stack identifying spots,
% and finding their centroids.

for i=1:numplanes
    
    % Convert to binary.
    filtered(:,:,i) = im2bw(log_spot(:,:,i),level);
    % Fill in holes in regions.
    filtered(:,:,i) = imfill(filtered(:,:,i),'holes');
    % Eliminate all regions smaller than 8 pixels.
    filtered(:,:,i) = bwareaopen(filtered(:,:,i),8);
    % Make a connectivity map of the regions.
    regions = bwconncomp(filtered(:,:,i),8);
    
    % Calc the coords of the boxes containing the regions.
    borders = regionprops(regions, 'BoundingBox');
    boxes = cat(1, borders.BoundingBox);
    
    % Filter ROIs whose peak intensity is below 
    % a threshold.
    
    numregs = size(boxes);  % Calculate number of regions.
    ROI_idx = [];           % To store indices of ROIs.
    
    % Loop through the regions to find spots.
    for j=1:numregs(1)
        
        try
            
            % Extract coords of the region to plot the spot.
            % BoundingBox gives coordinates as half pixels
            % so this corrects for that.
            box = boxes(j,:);
            xmin = box(1,1);
            xmin = xmin + 0.5;
            xmax = xmin+box(1,3)-1;
            ymin = box(1,2);
            ymin = ymin + 0.5;
            ymax = ymin+box(1,4)-1;
            
            % Get spot image from the reference image.
            spot_ref = imgauss_spot(ymin:ymax,xmin:xmax,i);
            
            % Calculate the max pixel intensity of the region.
            max_intensity = max(max(spot_ref));
            
            % Store indices/intensities of ROIs above threshold.
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
    
    % Create a map of the ROIs above intensity threshold.
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
    
    % Set thresholds for spread and area.
    spread_threshold = 25; % Early on, we were much more stringent (2.8), but later improvements to the code eliminated spurious spots through other means, and so we increased the threshold to 25.
    area_threshold = 5; % Spot must be at least this number of pixels.
    
    % Create a map of the ROIs below spread threshold.
    shape_idx = find((spread < spread_threshold) & (area > area_threshold));
    spots = ismember(labelmatrix(regions), shape_idx);

    % Add the good spots to the spots projection.
    % Ultimately, we'll only use spots that were seen
    % in two or more adjacent planes.
    spot_projection = spot_projection + spots;
    
    % Make a connectivity map of the spots.
    ROIs = bwconncomp(spots,8);
    ROImap(i) = ROIs;

    % Calc the coords of the boxes containing the regions.
    boundaries = regionprops(ROIs, 'BoundingBox');
    containers = cat(1, boundaries.BoundingBox);
    
    % Calculate the intensities of spots and their coordinates.
    
    numspots = size(containers);    % Calculate number of regions
    ROI_idx = [];                   % To store indices of ROIs
    ROI_data = [];                  % To store intensities
    
    % Loop through the ROIs to find spots for crude
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
        

        % Localize the spots and measure their intensities
        
        % Calculate local background using local pixels below
        % the threshold set previously.
        bg_ref = imgauss_spot(ymin-1:ymax+1,xmin-1:xmax+1,i);
        bg_mask = ~filtered(ymin-1:ymax+1,xmin-1:xmax+1,i);      
        filtered_bg = bg_ref.*bg_mask;
        bg_sum = sum(sum(filtered_bg));
        num_bg_pxls = sum(sum(bg_mask));
        bg_per_pxl = bg_sum/num_bg_pxls;
        
        % Calculate a refined intensity by subtracting
        % the local background.
        spot_ref = imgauss_spot(ymin:ymax,xmin:xmax,i);
        spot_mask = filtered(ymin:ymax,xmin:xmax,i);
        filtered_spot = spot_ref.*spot_mask;
        spot_sum = sum(sum(filtered_spot));
        num_spot_pxls = sum(sum(spot_mask));
        spot_bg = bg_per_pxl * num_spot_pxls;
        integrated_intensity = spot_sum - spot_bg;
        
        % Generate a background-subtracted image for
        % localization.
        bg_matrix = (spot_mask .* 0) + bg_per_pxl;
        subtracted_spot = (spot_ref - bg_per_pxl).*spot_mask;
        
        % Add the subtracted spot image to the final green image.
        finalspot(ymin:ymax,xmin:xmax,i) = finalspot(ymin:ymax,xmin:xmax,i) + subtracted_spot;
        
        % Localization by center-of-mass.
        X_hist=sum(subtracted_spot,1);
        N = size(X_hist);
        Y_hist=sum(subtracted_spot,2);
        Y_hist = rot90(Y_hist);
        M = size(Y_hist);
        X=1:N(2);
        Y=1:M(2);
        Xcom = sum(X.*X_hist)/sum(X_hist);
        Ycom = sum(Y.*Y_hist)/sum(Y_hist);
        
        % Correct for fitting within an ROI.
        Xcoord = Xcom + xmin - 1;
        Ycoord = Ycom + ymin - 1;
        
%         % Uncomment to see spot, mask, & localizations
%         if i < numImages
%             range = [-40 150] % optional display range
%             hold off
%             close all
%             imshow(subtracted_spot, [], 'InitialMagnification', 2000);
%             hold on, plot(Xcom, Ycom, 'b+')
%             hold off
%             wait = input('Hit enter/return to continue.');
%             saveas(gcf, strcat(dir, 'COM.tif'));
%             close
%         end
        
        ROI_idx = [ROI_idx
            j];
        
        ROI_data = [ROI_data
            j integrated_intensity Xcoord Ycoord];
        
    end;
    
    % Append the spot data for this plane to the spot
    % data for all the planes.
    numSpots = size(ROI_data);
    centroids(1:numSpots(1),1:numSpots(2),i) = ROI_data;

end

'Localized spots and measured their intensities'

%% View and save images of spots

% Uncomment to see all background-subtracted spots in each plane.
sumspots = zeros(size(finalspot(:,:,1)));
for i=1:numplanes
    sumspots = sumspots + finalspot(:,:,i);
end
maxint = max(max(sumspots));

imshow(sumspots, [0 maxint]);
colormap('Jet')
colorbar
% wait = input('Hit enter/return to continue.'); % Uncomment to hold.
saveas(gcf, strcat(dir, 'spot_sum.tif'));
close

clear bg_ref bg_sum borders boundaries box boxes container containers
clear filtered_bg filtered_spot i imtemp integrated_intensity
clear j kernelsize kernelspot max_intensity numSpots num_bg_pixels
clear num_spot_pxls numregs numspots perim_area perimeter regions
clear shape_idx spot_bg spot_mask spot_ref spot_sum spots
clear  spread_threshold subtracted_spot water log_spot

% Uncomment to see a projection of all the spots.
close all
figure;
hold on
imshow(spot_projection,[]);
hold off
% wait = input('Hit enter/return to continue.'); % Uncomment to hold.
saveas(gcf, strcat(dir, 'spot_projection.tif'));
close

'Saved images of spots'

%% Link spots across z

% Link spots visible in adjacent z-planes and
% create a list of non-redundant spots at their 
% max intensities. This part of the script will
% eliminate any spot that is not present in at
% least two adjacent planes.

% Make lists to temporarily store centroids.
centroids_previous = [];
centroids_current = [];

% Make lists to permanently store links and unique centroids.
centroids_links = [];
centroids_unique = [];
indexes = [];
indexes_unique = [];
uniqueCounter = 0;

% Loop through the centroids lists (1 list per plane).
for i=1:numplanes
    
    % Determine the size of the current centroids list.
    numSpots = size(centroids(:,3:4,i));
    tempcentroids = centroids(:,3:4,i);
    
    % Make a list of centroids for the plane.
    centroids_current = [];
    
    for j=1:numSpots(1)
        
        % Get the current centroid.
        centroid = tempcentroids(j,:);
        
        % Add it to the list, if it's not a placeholder centroid.
        if centroid ~= [0 0]
            
            centroids_current = [centroids_current
                centroid];
            
        end
        
    end
    
    % Link centroids & identify unique centroids.
    % NOTE: This will ignore spots that cannot
    % be detected in at least two adjacent planes.
    
    % If there were centroids in the previous plane.
    if centroids_previous 
        
        % And there are centroids in the current plane.
        if centroids_current % Can't happen the first time through loop.
            
            % Measure all pairwise distances between the centroids
            % in centroids_previous & centroids_current.
            pwds = pdist2(centroids_previous,centroids_current, 'euclidean');
            
            % Identify the indices of centroids that are within a
            % maximum distance.
            distance_threshold = 4; % 4 seems to work well for most purposes.
            
            % Collect the indices of linked spots.
            % idx_previous will contain the indices for centroids_previous.
            % idx_current will contain the indices for centroids_current.
            [idx_previous, idx_current] = find(pwds<distance_threshold);
            numCurrentLinks = size(idx_previous);
            
            % If there are linked centroids
            if idx_previous
                
                % Create lists for storing links.
                counters = [];
                planes_previous = [];
                index_previous = [];
                coords_previous = [];
                intensities_previous = [];
                planes_current = [];
                index_current = [];
                coords_current = [];
                intensities_current = [];
                
                % Make the list as long as the number of 
                % linked centroids.
                for j=1:numCurrentLinks
                    
                    % Give each spot a counter.
                    uniqueCounter = uniqueCounter + 1;
                    counters = [counters
                        uniqueCounter];
                    
                    plane_previous = i-1;
                    planes_previous = [planes_previous
                        plane_previous];
                    index_previous = [index_previous
                        plane_previous centroids(idx_previous(j),1,planes_previous(j))];
                    coord_previous = centroids(idx_previous(j),3:4,planes_previous(j));
                    coords_previous = [coords_previous
                        coord_previous];
                    intensity_previous = centroids(idx_previous(j),2,planes_previous(j));
                    intensities_previous =[intensities_previous
                        intensity_previous];
                      
                    plane_current = i;
                    planes_current = [planes_current
                        plane_current];
                    index_current = [index_current
                        plane_current centroids(idx_current(j),1,planes_current(j))];
                    coord_current = centroids(idx_current(j),3:4,planes_current(j));
                    coords_current = [coords_current
                        coord_current]; 
                    intensity_current = centroids(idx_current(j),2,planes_current(j));
                    intensities_current =[intensities_current
                        intensity_current];
                
                end
                
                % Append the new links to the big list.
                current_links = [counters planes_previous idx_previous intensities_previous coords_previous planes_current idx_current intensities_current coords_current];
                centroids_links = [centroids_links
                    current_links];
                
                indexes = [indexes
                    index_previous index_current];
                
                % Calculate the current number of links
                % and the number of unique centroids.
                numCumulativeLinks = size(centroids_links);
                numUnique = size(centroids_unique);
                
                % Initialize a flag.
                % 1: add it to centroids_unique
                % 0: pass
                flag = 1;
                
                % Loop over all accumulated links.
                for j=1:numCurrentLinks
                    
                    flag = 1; % Reset flag each time through.
                  
                    % Loop over all accumulated links.
                    for m=1:numCumulativeLinks(1)
                                                
                        % If a new centroid is found
                        % in the third pair of columns of 
                        % centroids_links
                        if centroids_links(m,10:11) == current_links(j,5:6)
                            
                            flag = 0; % Don't add it.
                            
                            % Go back over all the accumulated links.
                            for n=1:numCumulativeLinks(1)
                                
                                % When you find the matching spot
                                if centroids_links(n,5:6) == current_links(j,5:6)
                                    
                                    % Renumber it so that it has a matching
                                    % counter.
                                    centroids_links(n,1) = centroids_links(m,1);
                                                                    
                                end
                            
                            end
                            
                        end
                        
                    end
                    
                    % If it hasn't been linked yet, add it.
                    if flag == 1;
                        
                        centroids_unique = [centroids_unique
                            current_links(j,1:6)];
                        indexes_unique = [indexes_unique
                            indexes(j,:)];
                                            
                    end
                    
                end
                
            end
            
        end
        
    end
    
    % Store the current centroid list as the
    % previous centroid list for next round.
    centroids_previous = centroids_current;

end

'Linked spots across z'

%% Link gapped spots

% Determine number of unique/non-unique centroids.
numUnique = size(centroids_unique);
numLinks = size(centroids_links);

% Loop over unique centroids & link anything that was
% overlooked because there was a gap across z.
for i=1:numUnique(1)
    
    % Loop over all of the links
    for j=1:numLinks
        
        % if the current unique centroid & the current
        % centroid in the big list list don't have
        % the same index
        current_unique_centroid_idx = centroids_unique(i,1);
        current_linked_centroid_idx = centroids_links(j,1);
        if current_unique_centroid_idx ~= current_linked_centroid_idx
            
            % But they are within the distance threshold
            pwd = pdist2(centroids_links(j,10:11),centroids_unique(i,5:6), 'euclidean');
            distance_threshold = 4;
            if pwd < distance_threshold
                
                % Re-index all the observations that share the second 
                % spot index so that they match the current spot.
                [rows, col] = find(centroids_links(:,1) == current_linked_centroid_idx);
                numrows = size(rows);
                for k=1:numrows(1)
                    
                    centroids_links(rows(k),1) = centroids_unique(i,1);
                    
                end
                
                % Re-index the spot in centroids_unique
                % so the spot isn't repeated.
                [row, col] = find(centroids_unique(:,1) == current_linked_centroid_idx);
                centroids_unique(row,1) = current_unique_centroid_idx;
                
            end
            
        end
        
    end
    
end

spots_unique = [centroids_unique(1,:)];

% Loop over unique centroids & get rid of
% known redundant spots.
for i=2:numUnique(1)
    
    % Find any redundant spots.
    [rows, cols] = find(spots_unique(:,1) == centroids_unique(i,1));
        
    if isempty(rows)
        
        spots_unique = [spots_unique
            centroids_unique(i,:)];
    end
    
end

% Variable swap/cleanup.
clear centroids_unique
centroids_unique = spots_unique;
clear spots_unique

'Linked gapped spots'

%% Find the planes at which each spot is at its maximum intensity

% Determine number of unique/non-unique centroids.
numUnique = size(centroids_unique);
numLinks = size(centroids_links);

% Make a list to store the data.
current_spot = [];
spotmaxes = [];

% Loop over unique centroids.
for i=1:numUnique(1)
    
    % Add the data for the first observation
    % to a temporary list of data that will
    % contain all observations of the current
    % spot.
    current_spot = [centroids_unique(i,:)];
    current_indexes = [indexes_unique(i,1:2)];
    
    % Loop over all of the links.
    for j=1:numLinks
        
        % If the current link & the first observation
        % have the same spot index
        if centroids_unique(i,1) == centroids_links(j,1)
            
            % Add the linked data to the list of
            % all observations of the current spot.
            current_spot = [current_spot
                centroids_unique(i,1) centroids_links(j,7:11)];
            current_indexes = [current_indexes
                indexes(j,3:4)];
          
        end
        
    end
    
    % Find the plane in which the spot was at 
    % its max intensity.
    spotmax = find(current_spot(:,4) == max(current_spot(:,4)));
    
    % Add the data from that plane to the list
    % of spots at their max intensity planes.
    spotmaxes = [spotmaxes
        current_spot(spotmax,:)];
    
%     % Uncomment to recalculate the coords by Gaussian fitting.
%     plane = current_spot(spotmax,2);
% 	  spot = current_indexes(spotmax,2);
%     [Ycoord, Xcoord, err] = gaussfit2d_for_map_nuclear_spots(plane, spot, ROImap, current_indexes, imgauss, level, imfilt, filtered);
%     spotmaxes(i,5:6) = [Xcoord, Ycoord];

%     % Uncomment to see spot, mask, & localizations.
%     hold off
%     close all
%     imshow(autogain_uint16(subtracted_spot), 'InitialMagnification', 2000);
%     hold on, plot(XgaussTemp, YgaussTemp, 'b+')
%     hold off
%     wait = input('Hit enter/return to continue.');
    
end

'Found the planes at which each spot is at its maximum intensity'

%% Open nuclear boundary image and apply LOG and Gaussian filters

% Open the nucleus image, filter it using a 
% Laplacian of Gaussians filter, then 
% calculate the centroids and approximate 
% radii of the nuclei.

% Open the nucleus image.
image_nuc = openTiffStack(imR);

% If you have really bright dead cells, set a 
% dead_cell_threshold above the maximum expected intensity 
% of NSG1. Otherwise, imfindcircles() will only find dead 
% cells!
dead_cell_threshold = 110; % This must be customized for your experimental set-up
image_nuc(find(image_nuc > dead_cell_threshold)) = dead_cell_threshold;

% Convert the image to 8-bit.
image_nuc = uint8(image_nuc);

% Pad the nucleus image to match the dimensions of the spot image.
pad = 50;
[sizeY, sizeX, sizeZ] = size(image_nuc);
paddedY = sizeY+pad*2;
paddedX = sizeX+pad*2;
padded = zeros(paddedY,paddedX);
padded3D = zeros(paddedY,paddedX,sizeZ);

% Define the Laplacian of Gaussians (LOG) filter for 
% fitting the nuclei.
kernelsize_nuc = 20; % 20 works well with Nsg1-mCherry
sigma_nuc = 4.4; % 4.4 works well with Nsg1-mCherry
kernelnuc = -fspecial('log', kernelsize_nuc, sigma_nuc);

% Make a Gaussian filtered image.
% NOTE: The larger sigma is for Gaussian kernel, the smaller the 
% apparent nuclear radius will be, and the more peripheral any 
% given spot will appear. However, if sigma is too small, the 
% calculated boundary of the nucleus is impacted by statistical
% errors and background signal. So there is a trade-off between 
% precision and accuracy. (Ideally, you would use a negative 
% control to optimize this parameter empirically. The correct 
% negative control would be a subnuclear feature that is randomly 
% positioned within the nucleus, and  optimization would require 
% determining the value of sigma that gives a localization result 
% most close to random. When we began our studies, we weren't 
% aware of any subnuclear features that had been described with 
% such a random distribution. However in light of our results,
% we suggest that kinetochores detached in nocodazole may be a 
% good control of this type for future studies.)
gausskernelsize_nuc = 9; % 9 works well for Nsg1-mCherry
gausssigma_nuc = 6; % 6 works well for Nsg1-mCherry
kernelgauss = fspecial('gaussian', gausskernelsize_nuc, gausssigma_nuc);
imtempgauss = imfilter(double(image_nuc(:,:,:)), kernelgauss, 'symmetric', 'conv');
padded3D(pad+1:end-pad, pad+1:end-pad,:) = imtempgauss(:,:,:);
imgauss_nuc = padded3D;

clear padded3D imtempgauss 

'Opened nuclear boundary image and applied LOG and Gaussian filters'

%% Fit nuclei using imfindcircles()

% Create containers to store centers and radii.
% maxCircles is to create an arbitrarily-sized container
% based on an overestimation of the most circles you
% expect per image.
maxCircles = 300;
centers = zeros(maxCircles,2,numplanes) + 1;
radii = zeros(maxCircles,numplanes) + 1;

% Loop over each plane
for i=1:numplanes
    
    % Filter image with LOG
	imtemp = imfilter(double(image_nuc(:,:,i)), kernelnuc, 'symmetric', 'conv');

    % Resize the image after filtering to prevent edge artifacts
    padded = zeros(paddedY,paddedX);
    padded(pad+1:end-pad, pad+1:end-pad) = imtemp;
    log_nuc(:,:,i) = padded;
    
    % Find the circles in the image
    % [15 20] seems to be a good radius range
    % 'Sensitivity' 0.955 seems good
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

clear log_nuc imtemp padded image_nuc

'Fit nuclei using imfindcircles()'

%% Identify the nuclei to which each spot belongs.

% Make a list to pool all the spots that can be 
% matched to nuclei.
spots_nuclei = []; 
numSpots = size(spotmaxes);
spotsCount = 0;

% Loop over the spots
for i=1:numSpots(1)
    
    % Make a temporary list to pool the nuclei
    % from the current plane.
    circus = [];
    circus_below = [];
    circus_above = [];
    plane = spotmaxes(i,2);
    plane_below = plane - 1;
    plane_above = plane + 1;
        
    % Loop over the data defining the circles
    % in the current plane and above and below.
    for j=1:maxCircles
        
        % If the spot was not in the bottom plane
        if plane ~= 1
            
            % If there is data for a circle in the plane below
            if radii(j,plane_below) ~= 1
                
                % Add it to the temporary list.
                circus_below = [circus_below
                    centers(j,:,plane_below) radii(j,plane_below)];
                
            end
            
        end
        
        % If there is data for a circle in the same plane
        if radii(j,plane) ~= 1
            
            % Add it to the temporary list.
            circus = [circus
                centers(j,:,plane) radii(j,plane)];
            
        end
        
        % If the spot was not in the top plane
        if plane ~= numplanes
            
            % If there is data for a circle in the plane above
            if radii(j,plane_above) ~= 1
                
                % Add it to the temporary list.
                circus_above = [circus_above
                    centers(j,:,plane_above) radii(j,plane_above)];
                
            end
        
        end
    
    end
    
    % Create a flag to monitor if matching the spot
    % to a nucleus has been successful.
    flag = 0;
    
    % If there weren't any circles, put an arbitrary point
    % in the list that won't be near any spot.
    % This prevents an error at the next step.
    if size(circus) == [0 0]
        
        circus = [circus
            9999999 9999999 10];
        
    end

    % Calculate the distance between the current spot
    % and the centers of all detected nuclei.
    pwds = pdist2(spotmaxes(i,5:6),circus(:,1:2), 'euclidean');
    
    % Find the nearest nucleus in the plane.
    [col, nearest] = find(pwds == min(pwds));
    
    % Define a distance threshold.
    distance_threshold = 17; % 17 works well (with a pixel width of 0.065 um)
    
    % If the nearest nucleus is near enough
    if pwds(col,nearest) < distance_threshold
        
        % Find how many other nuclei are within the distance.
        [cols, neighbors] = find(pwds < distance_threshold);
        numNeighbors = size(neighbors);
        
        % If there aren't any others
        if numNeighbors(2) == 1
            
            % Then assign the spot to the nearest nucleus.
            spotsCount = spotsCount + 1;
            circ = circus(nearest,:);
            spots_nuclei = [spots_nuclei
                spotsCount spotmaxes(i,:) circ];
                        
        end
        
        % Set the flag to 1 if the spot was matched or it
        % could not be matched because more than one
        % nucleus within the distance.
        flag = 1;
        
    end
    
    % If the spot could not be linked to a nucleus
    if flag == 0
        
        % Set flags for if the spot could be linked to a
        % nucleus above or below.
        below = 0;
        above = 0;
        
        % Try to find the nucleus in the neighboring planes.
        
        try
                        
        	% Calculate the distance between the current spot
            % and the centers of all detected nuclei
            % in the plane below.
        	pwds_below = pdist2(spotmaxes(i,5:6),circus_below(:,1:2), 'euclidean');
            
            % Find the nearest nucleus above and below
            [col_below, nearest_below] = find(pwds_below == min(pwds_below));
            
            % If the nearest nucleus below is near enough
            if pwds_below(col_below,nearest_below) < distance_threshold
                
                % Find how many others are within the distance
                [cols_below, neighbors_below] = find(pwds_below < distance_threshold);
                numNeighbors_below = size(neighbors_below);
                
                % If there aren't any others
                if numNeighbors_below(2) == 1
                    
                    % Change the flag to indicate 'found.'
                    below = 1;
                    flag = 1;
                    
                end
                
            end
            
        end
        
        try
            
            % Calculate the distance between the current spot
            % and the centers of all detected nuclei
            % in the plane above.
            pwds_above = pdist2(spotmaxes(i,5:6),circus_above(:,1:2), 'euclidean');
            
            % Find the nearest nucleus above and below.
            [col_above, nearest_above] = find(pwds_above == min(pwds_above));
        
            % If the nearest nucleus below is near enough
            if pwds_above(col_above,nearest_above) < distance_threshold
                
                % Find how many others are within the distance.
                [cols_above, neighbors_above] = find(pwds_above < distance_threshold);
                numNeighbors_above = size(neighbors_above);
                
                % If there aren't any others
                if numNeighbors_above(2) == 1
                    
                    % Change the flag to indicate 'found.'
                    above = 1;
                    flag = 1;
                    
                end
                
            end
        
        end
              
        % If the spot could be linked to a nucleus above or below
        if (below + above) > 0
                                    
            spotsCount = spotsCount + 1;

            % If the spot was near a nucleus below
            if below
                
                % And the spot was also near a nucleus above
                if above;
            
                    % Then assign the spot to the average of the
                    % nearest nuclei above and below.
                    xcircus = (circus_below(nearest_below,1)+circus_above(nearest_above,1))/2;
                    ycircus = (circus_below(nearest_below,2)+circus_above(nearest_above,2))/2;
                    radcircus = (circus_below(nearest_below,3)+circus_above(nearest_above,3))/2;
                    circ = [xcircus ycircus radcircus];
                    spots_nuclei = [spots_nuclei
                        spotsCount spotmaxes(i,:) circ];
            
                % If it was only found below
                else
                    
                    % Assign the spot to the nearest nucleus below.
                    circ = circus_below(nearest_below,:);
                    spots_nuclei = [spots_nuclei
                        spotsCount spotmaxes(i,:) circ];
                    
                end
                            
            % If it was only found above
            else
                
                % Assign the spot to the nearest nucleus above.
                circ = circus_above(nearest_above,:);
                spots_nuclei = [spots_nuclei
                    spotsCount spotmaxes(i,:) circ];

            end
                        
        end
                
    end
    
end

'Identified the nucleus to which each spot belongs.'

%% Group spots within nuclei

% Figure out which spots were in the same nucleus
% by finding out how close the nuclei are to each other.

% Determine how many spots there are.
numSpots = size(spots_nuclei);
% Extract the coords of all the nuclei.
nuclei = spots_nuclei(:,8:9);
% Decide how far apart the centers of nuclei can be.
distance_threshold = 8;

% Loop over the spots.
for i=1:numSpots(1)
    
    % Get the coords of nucleusA.
    nucleusA = nuclei(i,:);
    
    % Loop over the spots again.
    for j=1:numSpots
        
        % This assures you don't repeat comparisons
        % or compare a spot to itself.
        if j > i
            
            % Get the coords of nucleusB.
            nucleusB = nuclei(j,:);
            % Determine the distance between the nuclei.
            pwd = pdist2(nucleusA,nucleusB,'euclidean');
            
            % If the nuclei are nearby.
            if pwd < distance_threshold

                % Change the index of nucleusB to match 
                % the index of nucleusA.
                spots_nuclei(j,1) = spots_nuclei(i,1);                
   
            end
    
        end
        
    end
    
end

% Check to see that the spots have been grouped appropriately.

% Determine the number of unique nuclei.
uniqulii = unique(spots_nuclei(:,1));
numUniqulii = size(uniqulii);

% Create a temporary container to store all the spot
% and nucleus data that maps to the same cell.
uniqulum = [];

% Create a permanent container to store the spots grouped.
grouped_spots_nuclei = [];

% Set threshold for distance between unique spots.
distance_threshold = 4;

% Loop over the list of unique nuclei.
for i=1:numUniqulii(1)
    
    uniqulum = [];

    % Get the index of the current nucleus.
    uniqulus = uniqulii(i);
    
    % Retreive the indices for spots that map to 
    % the current nucleus.
    [rows,cols] = find(spots_nuclei(:,1) == uniqulus);

    % Create a list with all the data for the current
    % nucleus.
    uniqulum = spots_nuclei(rows,:);
    numspots = size(uniqulum);
    
    % Renumber spots within the current nucleus.
    for j=1:numspots(1)
        
                uniqulum(j,4) = j;
                
    end
    
    % Loop over the spots.
    for j=1:numspots(1)
        
        % Get the coords of the spot.
        coordsA = [uniqulum(j,6),uniqulum(j,7)];
        
        % Loop over the spots
        for k=1:numspots(1)
            
            % Skipping redundant and self comparisons.
            if k>j
                
                % Get the coords of the spot.
                coordsB = [uniqulum(k,6),uniqulum(k,7)];
                
                % Measure the pairwise distance
                pwd = pdist2(coordsA,coordsB,'euclidean');
                
                % If the spots are nearby
                if pwd < distance_threshold
                    
                    % Change the index of spotB to match
                    % the index of spotA
                    uniqulum(k,4) = j;
                    
                end
                    
            end
            
        end
                
    end
    
    % Create an empty list to store the final spot
    % data for this nucleus.
    uniqulus = [];
    
    % Find the list of unique spots
    uniqspots = unique(uniqulum(:,4));
    % And determine how many there are.
    numUniqSpots = size(uniqspots);
    
    % Loop over the spots again to eliminate redundant ones
    % by finding and keeping the max intensity one.
    for j=1:numUniqSpots(1)
        instances = find(uniqulum(:,4) == uniqspots(j));
        peak = find(uniqulum(:,5) == max(uniqulum(instances,5)));
        spot = uniqulum(peak,:);
        uniqulus = [uniqulus
            spot(1,:)];
        
    end
    
    % Generate indices for the nucleus number.
    number = cols(1:numUniqSpots) .* i;
    
    % Add the data to the list with the new numbers.
    grouped_spots_nuclei = [grouped_spots_nuclei
        number uniqulus(:,2:10)];

%     % Uncomment to see the localizations of all the
%     % spots within a nucleus at their max intensity
%     % planes superimposed over the MIP.
%     imshow(MIP_spot,[]);
%     hold on, plot(uniqulus(:,6) + 1, uniqulus(:,7) + 1,'b+');
%     hold on, plot(uniqulus(:,8) + 1, uniqulus(:,9) + 1,'ro');
%     h = viscircles(uniqulus(:,8:9) + 1,uniqulus(:,10),'LineWidth',1);
%     hold off
%     wait = input('Hit enter/return to continue.');
    
end

' Grouped spots within nuclei'

%% Save images of spots at max intensity planes

% Uncomment to see the localizations of the spots
% at their max intensity planes superimposed over
% the maximum intensity projection.

figure;
imshow(MIP_spot,[]);
hold on, plot(grouped_spots_nuclei(:,6), grouped_spots_nuclei(:,7), 'b+', 'LineWidth', .25, 'MarkerSize',2); 
hold off
% Prevent it from making the background white.
set(gcf, 'InvertHardCopy', 'off');
% wait = input('Hit enter/return to continue.'); % Uncomment to hold.
saveas(gcf, strcat(dir, 'spots_MIP.tif'));
close

figure;
imshow(MIP_spot,[], 'InitialMagnification', 500);
hold on, plot(grouped_spots_nuclei(:,6), grouped_spots_nuclei(:,7), 'b+', 'LineWidth', .25, 'MarkerSize', 2); 
hold on, plot(grouped_spots_nuclei(:,8), grouped_spots_nuclei(:,9),'ro', 'LineWidth', .25, 'MarkerSize', 2);
h = viscircles(grouped_spots_nuclei(:,8:9),grouped_spots_nuclei(:,10), 'LineWidth', .25, 'EdgeColor', 'w');
hold off
% Prevent it from making the background white.
set(gcf, 'InvertHardCopy', 'off');
% wait = input('Hit enter/return to continue.'); % Uncomment to hold.
saveas(gcf, strcat(dir, 'spots_nuclei_MIP.tif'));
close

'Saved images of spots at max intensity planes'

%% Adjust final display ranges.

% These settings will not affect any of the quantitative results,
% only the appearance during the accept/reject step.
% The display ranges you choose should probably be set the same 
% within a given experiment, but they can be adjusted for 
% different types of spots in order to allow you to best interpret 
% the images, since some spots have different intensity ranges.
% (For instance, TetO/TetR spots are usually much brighter than 
% detached kinetochores labeled with Mtw1-3xGFP.)

% NOTE: This may need be optimized for a particular set-up.
finalspot = double(finalspot);
mostoften = mode(mode(mode(finalspot)));
finalspot(find(finalspot < mostoften)) = mostoften;
finalspot = finalspot - mostoften;
finalspot = autogain_uint16(finalspot);

% NOTE: This needs to be optimized for a particular set-up,
% since subtract defines an arbitrary threshold.
finalnuc = double(imgauss_nuc(:,:,:));
finalnuc = autogain_uint16(finalnuc);
subtract = 5000; % Arbitrary. Depends on experimental set-up.
finalnuc = (finalnuc - subtract);
redmax = max(max(max(finalnuc)));
finalnuc = finalnuc.*(65535/redmax);

'Adjusted final display ranges.'

%% Supervise localization and save final results

% Create a place to keep track of the final results.
final_results = [];
finalcount = 0;
pixelsize_nm = 64.5;

numSpots = size(grouped_spots_nuclei);

% Loop over the spots.
for i=1:numSpots(1)
    
    % Get the xcoords and ycoords from the spot and nucleus.
    xs = [grouped_spots_nuclei(i,6),grouped_spots_nuclei(i,8)];
    ys = [grouped_spots_nuclei(i,7),grouped_spots_nuclei(i,9)];
    
    % Fit a line to the points (y = mx + b).
    line = polyfit(xs,ys,1);
    slope = line(1);
    intercept = line(2);

    % Specify the length of the line to trace.
    length = 20;
    % Calculate how many x pixels the line must cross
    % to make the a line of the specified length
    % based on the slope of the line.
    window = length/sqrt(1+(slope*slope));
    % Specify the number of points to compute on the line.
    sample = 100;
    % Specify the x spacing for the desired number of points.
    increment = window/sample;
    % Get the x coord of the nuclear center.
    nucleusx = grouped_spots_nuclei(i,8);
    
    % Default is to assume the spot is to the left of the
    % the center of the nucleus (true 50% of the time).
    right = 0;
    
    % If the spot is actually to the right.
    if grouped_spots_nuclei(i,6) < grouped_spots_nuclei(i,8)
        
        % Flag it.
        right = 1;
        
    end
    
    % If the spot was on the left.
    if right == 0
        
        % Sample from the center of the nucleus
        % to the right.
        xrange = nucleusx:increment:nucleusx+window;

    else
        
        % Sample from the center of the nucleus 
        % to the left.
        xrange = nucleusx-window:increment:nucleusx;

    end
    
    % Evaluate the equation with the specified range.
    yfit = slope*xrange + intercept;
    
    % Determine which plane the spot was in.
    plane = grouped_spots_nuclei(i,3);
    
    % Measure the intensities along the line.
    prof_nuc = improfile(imgauss_nuc(:,:,plane),xrange,yfit);
    prof_spot = improfile(imgauss_spot(:,:,plane),xrange,yfit);
        
    % Find the coordinate of maximum intensity.
    boundarypeak = find(prof_nuc == max(prof_nuc));
    boundarycoordX = mean(xrange(boundarypeak));
    boundarycoordY = slope*boundarycoordX + intercept;
    % Adjust the line profiles so all the points are above zero.
    prof_nuc = prof_nuc - min(prof_nuc);
    prof_spot = prof_spot - min(prof_spot);
    % Adjust the line profile of the spots data so that it has
    % the same maximum as the line profile of the nucleus data.
    ratio = max(prof_spot)/max(prof_nuc);
    prof_spot = prof_spot ./ ratio;
        
    % Go through the spots one by one to validate them.
    
    % Recalculate the coords to correct the remaining red:green offset.
    offsetspotcoords = [grouped_spots_nuclei(i,7) + offsetYremainder, grouped_spots_nuclei(i,6) + offsetXremainder];
    % Get the coords for the fitted nuclear boundary.
    boundarycoords = [boundarycoordY, boundarycoordX];
    centercoords = [grouped_spots_nuclei(i,9), grouped_spots_nuclei(i,8)];
    % Calculate the distance between the spot and
    % the nuclear boundary.
    boundarydistance = pdist2(offsetspotcoords,boundarycoords,'euclidean');
    % Calculate the distance between the spot and
    % the center of the nucleus.
    centerdistance = pdist2(offsetspotcoords,centercoords,'euclidean');
    % Calculate the nuclear radius.
    radius = pdist2(boundarycoords,centercoords,'euclidean');
    
    % Filter out all nuclei with aberrantly small radii
    % or that peak at the last pixel of the line profile.
    
    if (radius > 3) & (radius <= 19) % May need optimization.
    
        % Convert distances to nm.
        boundarydistance_nm = boundarydistance * pixelsize_nm;
        centerdistance_nm = centerdistance * pixelsize_nm;
        radius_nm = radius * pixelsize_nm;
        
        % Determine if the spot is in zone 1, 2, or 3.
        % If the spot is farther from the center than the fitted nuclear
        % boundary (correcting for the final offset, which is necessarily
        % less than 1 pixel in either axis (or the square root of 2
        % pixels in both axes).
        if (boundarydistance + centerdistance) > (radius + sqrt(2))
            
            zone = 1;
            % The spot is on the outside of the boundary
            % so switch its sign to negative.
            boundarydistance_nm = -boundarydistance_nm;
            
        elseif boundarydistance < (0.184*radius)
            
            zone = 1;
            
        elseif boundarydistance < (0.422*radius)
            
            zone = 2;
            
        elseif boundarydistance >= (0.422*radius)
            
            zone = 3;
            
        end
        
        close all
        
        figure;
        
        % Display the line profiles and Gaussian fit.
        subplot(1,2,1);
        
        % Set the y-axis to adjust the range range.
        axismax = max(prof_nuc)*1.5;
        
        if axismax <= 0
            axismax = 1;
        end
        
        hold on, axis([min(xrange) max(xrange) 0 axismax]);
        % Make the plot background black.
        set(subplot(2,2,[1 3]),'Color','Black')
        % Show the position of the nuclear center as a white circle.
        hold on, plot(grouped_spots_nuclei(i,8), 0, 'wo', 'LineWidth', 2, 'MarkerSize',12);
        % Plot the spot profile along the line as a green line.
        hold on, plot(xrange,prof_spot,'g', 'LineWidth', 2);
        % Show the position of the spot as a blue x.
        hold on, plot(grouped_spots_nuclei(i,6),0,'bx', 'LineWidth', 2, 'MarkerSize',12);
        % Plot the nuclear profile along the line as a magenta line.
        hold on, plot(xrange,prof_nuc,'m', 'LineWidth', 2);
        % Show the position of the nuclear edge as a red x.
        hold on, plot(boundarycoordX,0,'rx', 'LineWidth', 2, 'MarkerSize',12);
        hold off
        
        % Set the width of the image you want to see
        % Width will be 2 * dimensions.
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
        
        % Calculate the display range for the nucleus.
        horiz_range = midx_floor-dimensions:midx_floor+dimensions;
        vert_range = midy_floor-dimensions:midy_floor+dimensions;
        
        % Display all the spots from the same nucleus.
        subplot(2,2,2);
        hold on     
        imshow(MIP_spot(vert_range,horiz_range),[]);
        samenucleus = find(grouped_spots_nuclei(:,1) == grouped_spots_nuclei(i,1));
        hold on, plot(grouped_spots_nuclei(samenucleus,6)-adj_x, grouped_spots_nuclei(samenucleus,7)-adj_y,'b+', 'LineWidth', 1, 'MarkerSize',4);
        hold on, plot(grouped_spots_nuclei(samenucleus,8)-adj_x, grouped_spots_nuclei(samenucleus,9)-adj_y,'ro', 'LineWidth', 1, 'MarkerSize',4);
        h = viscircles([grouped_spots_nuclei(samenucleus,8)-adj_x,grouped_spots_nuclei(samenucleus,9)-adj_y],grouped_spots_nuclei(samenucleus,10), 'LineWidth', .25, 'EdgeColor', 'w');
        hold off

        % Display the nucleus and spot along with the localization data
        % and the line along which the profile is taken.
        subplot(2,2,4);
        hold on
        
        % Extract the image of the nucleus and spots from the
        % reference image based on the calculations above.
        red = finalnuc(vert_range,horiz_range,plane);
        green = finalspot(vert_range,horiz_range,plane);
        
        % Display the color image.
        colorim = cat(3,red,green,red);
        peakint = max(prof_nuc);
        
        if peakint <= 0
            peakint = 1
        end
            
        imshow(colorim, [0 peakint]);
     
        % Show the position of the spot as a blue plus.
        hold on, plot(grouped_spots_nuclei(i,6)-adj_x, grouped_spots_nuclei(i,7)-adj_y,'b+', 'LineWidth', 2, 'MarkerSize',8);
        % Show the position of the nuclear center as a white circle.
        hold on, plot(grouped_spots_nuclei(i,8)-adj_x, grouped_spots_nuclei(i,9)-adj_y,'wo', 'LineWidth', 2, 'MarkerSize',8);
        % Show the axis of the line scan as a white line
        hold on, plot(xrange-adj_x,yfit-adj_y, 'w');
        % Show the position of the Gaussian-fit nuclear edge as a red plus.
        hold on, plot(boundarycoordX-adj_x, boundarycoordY-adj_y, 'r+', 'LineWidth', 2, 'MarkerSize',8);
        
        % Show the fit of imfindcircles as a white circle.
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
        
        % Allow supervision of spot inclusion.
        reject = input('Hit enter/return to keep.  Type any number to reject.\n');
%         reject = []
        keep = size(reject);
        
        % Eliminate the spot if the cell was on the 
        % border of the image.
        if (sum(red(1,:)) * sum(red(41,:)) * sum(red(:,1)) * sum(red(:,41))) == 0
            keep(1) = 1
        end
        
        % If the user says yes, store the data.
        if keep(1) == 0
            % Update the number of spots in the data.
            finalcount = finalcount + 1;
            
            % Group the spot data with the edgecoords and distance.
            finalspot = [grouped_spots_nuclei(i,:) boundarycoordX boundarycoordY boundarydistance_nm radius_nm zone];
            % Update the corrected spot coords in the data.
            finalspot(6:7) = [offsetspotcoords(2) offsetspotcoords(1)];
            % Replace the index of the spot within the plane
            % with the overall spot number.
            finalspot(4) = finalcount;
            
            % Append the results.
            final_results = [final_results
                stack finalspot];
            
            % Save the figure.
            % Prevent it from making the background white.
            set(gcf, 'InvertHardCopy', 'off');
            saveas(gcf, strcat(dir, num2str(finalcount),'.tif'));
            hold off
            
        end
    
    end

    
    close all
    
end

'Finished.'