function [centroids_unique, centroids_links, indexes, indexes_unique] = link_centroids2(numplanes, centroids)

%% Nathaniel I. Krefman 2014-05-06

% Links spots visible in adjacent z-planes. Will eliminate any 
% spot that is not present in at  least two adjacent planes.

% Input variables:
%	numplanes:        number of planes in the image
%	centroids:        output from find_centroids.m

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   distance_threshold

% Other file requirements:
%       None

%%

% Make lists to temporarily store centroids
centroids_previous = [];
centroids_current = [];

% Make lists to permanently store links and unique centroids
centroids_links = [];
centroids_unique = [];
indexes = [];
indexes_unique = [];
uniqueCounter = 0;

% Loop through the centroids lists (1 list per plane)
for i=1:numplanes
    
    % Determine the size of the current centroids list
    numSpots = size(centroids(:,5:6,i));
    tempcentroids = centroids(:,5:6,i);
    
    % Make a list of centroids for the plane
    centroids_current = [];
    
    for j=1:numSpots(1)
        
        % Get the current centroid
        centroid = tempcentroids(j,:);
        
        % Add it to the list, if it's
        % not a placeholder centroid
        if centroid ~= [0 0]
            
            centroids_current = [centroids_current
                centroid];
            
        end
        
    end
    
    % Link centroids & identify unique centroids
    % NOTE: This will ignore spots that cannot
    % be detected in at least two adjacent planes.

    % If there were centroids in the previous plane
    if centroids_previous 
        
        % And there are centroids in the current plane
        if centroids_current % Can't happen the first time through loop
            
            % Measure all pairwise distances between
            % the centroids in centroids_previous & centroids_current.
            pwds = pdist2(centroids_previous,centroids_current, 'euclidean');
            
            % Identify the indices of centroids that are 
            % within a maximum distance.
            % NOTE: This might fail if a spot in one plane is
            % near 2 spots in an adjacent plane.
            distance_threshold = 4;
            
            % Collect the indices of linked spots
            % idx_previous will contain the indices for centroids_previous
            % idx_current will contain the indices for centroids_current
            [idx_previous, idx_current] = find(pwds<distance_threshold);
            numCurrentLinks = size(idx_previous);
            
            % If there are linked centroids
            if idx_previous
                
                % Create lists for storing links
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
                % linked centroids
                for j=1:numCurrentLinks
                    
                    % Give each spot a counter
                    uniqueCounter = uniqueCounter + 1;
                    counters = [counters
                        uniqueCounter];
                    plane_previous = i-1;
                    planes_previous = [planes_previous
                        plane_previous];
                    index_previous = [index_previous
                        plane_previous centroids(idx_previous(j),1,planes_previous(j))];
                    coord_previous = centroids(idx_previous(j),5:6,planes_previous(j));
                    coords_previous = [coords_previous
                        coord_previous];
                    intensity_previous = centroids(idx_previous(j),2:4,planes_previous(j));
                    intensities_previous =[intensities_previous
                        intensity_previous];
                      
                    plane_current = i;
                    planes_current = [planes_current
                        plane_current];
                    index_current = [index_current
                        plane_current centroids(idx_current(j),1,planes_current(j))];
                    coord_current = centroids(idx_current(j),5:6,planes_current(j));
                    coords_current = [coords_current
                        coord_current]; 
                    intensity_current = centroids(idx_current(j),2:4,planes_current(j));
                    intensities_current =[intensities_current
                        intensity_current];
                
                end
            
                % Append the new links to the big list
                current_links = [counters planes_previous idx_previous intensities_previous coords_previous planes_current idx_current intensities_current coords_current];
                centroids_links = [centroids_links
                    current_links];
                indexes = [indexes
                    index_previous index_current];
                
                % Calculate the current number of links
                % and the number of unique centroids
                numCumulativeLinks = size(centroids_links);
                numUnique = size(centroids_unique);
                
                % Initialize a flag
                % 1: add it to centroids_unique
                % 0: pass
                flag = 1;
                
                % Loop over all accumulated links
                for j=1:numCurrentLinks
                    
                    flag = 1; % reset flag each time through
                  
                    % Loop over all accumulated links
                    for m=1:numCumulativeLinks(1)
                                                
                        % If a new centroid is found
                        % in the third pair of columns of 
                        % centroids_links
                        if centroids_links(m,14:15) == current_links(j,7:8)
                            
                            flag = 0; % Don't add it
                            
                            % Go back over all the accumulated links
                            for n=1:numCumulativeLinks(1)
                                
                                % When you find the matching spot
                                if centroids_links(n,7:8) == current_links(j,7:8)
                                    
                                    % Renumber it so that it has a matching
                                    % counter
                                    centroids_links(n,1) = centroids_links(m,1);
                                                                    
                                end
                            
                            end
                            
                        end
                        
                    end
                    
                    % If it hasn't been linked yet, add it
                    if flag == 1;
                        
                        centroids_unique = [centroids_unique
                            current_links(j,1:8)];
                        indexes_unique = [indexes_unique
                            indexes(j,:)];
                                            
                    end
                    
                end
                
            end
                                        
        end
        
    end
    
    % Store the current centroid list as the
    % previous centroid list for next round
    centroids_previous = centroids_current;

end

%%

% Determine number of unique/non-unique centroids
numUnique = size(centroids_unique);
numLinks = size(centroids_links);

% Loop over unique centroids & link anything that was
% overlooked because there was a gap across z
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
            pwd = pdist2(centroids_links(j,14:15),centroids_unique(i,7:8), 'euclidean');
            distance_threshold = 4;
            if pwd < distance_threshold
                
                % Re-index all the observations that share the second 
                % spot index so that they match the current spot
                [rows, col] = find(centroids_links(:,1) == current_linked_centroid_idx);
                numrows = size(rows);
                for k=1:numrows(1)
                    
                    centroids_links(rows(k),1) = centroids_unique(i,1);
                    
                end
                
                % Re-index the spot in centroids_unique
                % so the spot isn't repeated
                [row, col] = find(centroids_unique(:,1) == current_linked_centroid_idx);
                centroids_unique(row,1) = current_unique_centroid_idx;
                
            end
            
        end
        
    end
    
end

spots_unique = [centroids_unique(1,:)];

% Loop over unique centroids & get rid of
% known redundant spots
for i=2:numUnique(1)
    
    % Find any redundant spots
    [rows, cols] = find(spots_unique(:,1) == centroids_unique(i,1));
        
    if isempty(rows)
        
        spots_unique = [spots_unique
            centroids_unique(i,:)];
    end
    
end

% Variable swap/cleanup
clear centroids_unique
centroids_unique = spots_unique;
clear spots_unique