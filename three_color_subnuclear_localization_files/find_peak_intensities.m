function spotmaxes = find_peak_intensities(centroids_unique, centroids_links, indexes, indexes_unique)

%% Nathaniel I. Krefman 2014-05-06

% Finds the planes at which each spot is at its maximum intensity.

% Input variables:
%	centroids_unique:     output from find_centroids.m
%	centroids_links:      output from find_centroids.m
%	indexes:              output from find_centroids.m
%	indexes_unique:       output from find_centroids.m

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   None

% Other file requirements:
%       None

%%

% Determine number of unique/non-unique centroids
numUnique = size(centroids_unique);
numLinks = size(centroids_links);

% Make a list to store the data
current_spot = [];
spotmaxes = [];

% Loop over unique centroids
for i=1:numUnique(1)
    
    % Add the data for the first observation
    % to a temporary list of data that will
    % contain all observations of the current
    % spot.
    current_spot = [centroids_unique(i,:)];
    current_indexes = [indexes_unique(i,1:2)];
    
    % Loop over all of the links
    for j=1:numLinks
        
        % If the current link & the first observation
        % have the same spot index
        if centroids_unique(i,1) == centroids_links(j,1)
            
            % Add the linked data to the list of
            % all observations of the current spot
            current_spot = [current_spot
                centroids_unique(i,1) centroids_links(j,7:11)];
            current_indexes = [current_indexes
                indexes(j,3:4)];
          
        end
        
    end
    
    % Find the plane in which the spot was at 
    % its max intensity
    spotmax = find(current_spot(:,4) == max(current_spot(:,4)));
    
    % Add the data from that plane to the list
    % of spots at their max intensity planes
    spotmaxes = [spotmaxes
        current_spot(spotmax,:)];
    
end
