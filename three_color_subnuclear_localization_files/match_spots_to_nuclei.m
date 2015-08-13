function spots_nuclei = match_spots_to_nuclei(numplanes, spotmaxes, maxCircles, centers, radii)

%% Nathaniel I. Krefman 2014-05-06

% Identifies the nuclei to which each spot belongs.

% Input variables:
%   numplanes:        number of planes in the image
%   spotmaxes:        output from find_peak_intensities.m
%   maxCircles:       output from fit_nuclei.m
%   centers:          output from fit_nuclei.m
%   radii:            output from fit_nuclei.m

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   distance_threshold

% Other file requirements:
%       None

%%

% Make a list to pool all the spots that can be 
% matched to nuclei.
spots_nuclei = []; 
numSpots = size(spotmaxes);
spotsCount = 0;

% Loop over the spots
for i=1:numSpots(1)
    
    % Make a temporary list to pool the nuclei
    % from the current plane
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
                
                % Add it to the temporary list
                circus_below = [circus_below
                    centers(j,:,plane_below) radii(j,plane_below)];
                
            end
            
        end
        
        % If there is data for a circle in the same plane
        if radii(j,plane) ~= 1
            
            % Add it to the temporary list
            circus = [circus
                centers(j,:,plane) radii(j,plane)];
            
        end
        
        % If the spot was not in the top plane
        if plane ~= numplanes
            
            % If there is data for a circle in the plane above
            if radii(j,plane_above) ~= 1
                
                % Add it to the temporary list
                circus_above = [circus_above
                    centers(j,:,plane_above) radii(j,plane_above)];
                
            end
        
        end
    
    end
    
    % Create a flag to monitor if matching the spot
    % to a nucleus has been successful
    flag = 0;
    
    % If there weren't any circles, put an arbitrary point
    % in the list that won't be near any spot.
    % This prevents an error at the next step.
    if size(circus) == [0 0]
        
        circus = [circus
            9999999 9999999 10];
        
    end

    % Calculate the distance between the current spot
    % and the centers of all detected nuclei
    pwds = pdist2(spotmaxes(i,5:6),circus(:,1:2), 'euclidean');
    
    % Find the nearest nucleus in the plane
    [col, nearest] = find(pwds == min(pwds));
    
    % Define a distance threshold.
    % 17 seems to be good, but this could be based on data
    % once I have more.
    distance_threshold = 17;
    
    % If the nearest nucleus is near enough
    if pwds(col,nearest) < distance_threshold
        
        % Find how many other nuclei are within the distance
        [cols, neighbors] = find(pwds < distance_threshold);
        numNeighbors = size(neighbors);
        
        % If there aren't any others
        if numNeighbors(2) == 1
            
            % Then assign the spot to the nearest nucleus
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
        % nucleus above or below
        below = 0;
        above = 0;
        
        % Try to find the nucleus in the neighboring planes
        
        try
                        
        	% Calculate the distance between the current spot
            % and the centers of all detected nuclei
            % in the plane below
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
                    
                    % Change the flag to indicate 'found'
                    below = 1;
                    flag = 1;
                    
                end
                
            end
            
        end
        
        try
            
            % Calculate the distance between the current spot
            % and the centers of all detected nuclei
            % in the plane above
            pwds_above = pdist2(spotmaxes(i,5:6),circus_above(:,1:2), 'euclidean');
            
            % Find the nearest nucleus above and below
            [col_above, nearest_above] = find(pwds_above == min(pwds_above));
        
            % If the nearest nucleus below is near enough
            if pwds_above(col_above,nearest_above) < distance_threshold
                
                % Find how many others are within the distance
                [cols_above, neighbors_above] = find(pwds_above < distance_threshold);
                numNeighbors_above = size(neighbors_above);
                
                % If there aren't any others
                if numNeighbors_above(2) == 1
                    
                    % Change the flag to indicate 'found'
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
                    % nearest nuclei above and below
                    xcircus = (circus_below(nearest_below,1)+circus_above(nearest_above,1))/2;
                    ycircus = (circus_below(nearest_below,2)+circus_above(nearest_above,2))/2;
                    radcircus = (circus_below(nearest_below,3)+circus_above(nearest_above,3))/2;
                    circ = [xcircus ycircus radcircus];
                    spots_nuclei = [spots_nuclei
                        spotsCount spotmaxes(i,:) circ];
            
                % If it was only found below
                else
                    
                    % Assign the spot to the nearest nucleus below
                    circ = circus_below(nearest_below,:);
                    spots_nuclei = [spots_nuclei
                        spotsCount spotmaxes(i,:) circ];
                    
                end
                            
            % If it was only found above
            else
                
                % Assign the spot to the nearest nucleus above
                circ = circus_above(nearest_above,:);
                spots_nuclei = [spots_nuclei
                    spotsCount spotmaxes(i,:) circ];

            end
                        
        end
                
    end
    
end