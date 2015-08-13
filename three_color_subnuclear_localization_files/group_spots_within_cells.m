function grouped_spots_nuclei = group_spots_within_cells(spots_nuclei)

%% Nathaniel I. Krefman 2014-05-06

% Figures out which spots were in the same nucleus
% by finding out how close the nuclei are to each 
% other. Checks to see that the spots have been 
% grouped appropriately.

% Input variables:
%   spots_nuclei:     output from match_spots_to_nuclei.m

% We suggest adjusting the following variables in this function 
% for your particular set-up:
%   distance_threshold (NOTE: Used twice in different ways.)

% Other file requirements:
%       None

%%

% Figure out which spots were in the same nucleus
% by finding out how close the nuclei are to each other.

% Determine how many spots there are
numSpots = size(spots_nuclei);
% Extract the coords of all the nuclei
nuclei = spots_nuclei(:,8:9);
% Decide how far apart the centers of nuclei can be
distance_threshold = 8;

% Loop over the spots
for i=1:numSpots(1)
    
    % Get the coords of nucleusA
    nucleusA = nuclei(i,:);
    
    % Loop over the spots again
    for j=1:numSpots
        
        % This assures you don't repeat comparisons
        % or compare a spot to itself
        if j > i
            
            % Get the coords of nucleusB
            nucleusB = nuclei(j,:);
            % Determine the distance between the nuclei
            pwd = pdist2(nucleusA,nucleusB,'euclidean');
            
            % If the nuclei are nearby
            if pwd < distance_threshold

                % Change the index of nucleusB to match 
                % the index of nucleusA
                spots_nuclei(j,1) = spots_nuclei(i,1);                
   
            end
    
        end
        
    end
    
end

% Check to see that the spots have been grouped appropriately

% Determine the number of unique nuclei
uniqulii = unique(spots_nuclei(:,1));
numUniqulii = size(uniqulii);

% Create a temporary container to store all the spot 
% and nucleus data that maps to the same cell
uniqulum = [];
% Create a permanent container to store the spots grouped
grouped_spots_nuclei = [];

% Set threshold for distance between unique spots
distance_threshold = 4;

% Loop over the list of unique nuclei
for i=1:numUniqulii(1)
    
    uniqulum = [];

    % Get the index of the current nucleus
    uniqulus = uniqulii(i);
    
    % Retreive the indices for spots that map to 
    % the current nucleus
    [rows,cols] = find(spots_nuclei(:,1) == uniqulus);

    % Create a list with all the data for the current
    % nucleus
    uniqulum = spots_nuclei(rows,:);
    numspots = size(uniqulum);
    
    % Renumber spots within the current nucleus
    for j=1:numspots(1)
        
                uniqulum(j,4) = j;
                
    end
    
    % Loop over the spots
    for j=1:numspots(1)
        
        % Get the coords of the spot
        coordsA = [uniqulum(j,6),uniqulum(j,7)];
        
        % Loop over the spots
        for k=1:numspots(1)
            
            % Skipping redundant and self comparisons
            if k>j
                
                % Get the coords of the spot
                coordsB = [uniqulum(k,6),uniqulum(k,7)];
                
                % Measure the PWD
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
    % data for this nucleus
    uniqulus = [];
    
    % Find the list of unique spots 
    uniqspots = unique(uniqulum(:,4));
    % And determine how many there are
    numUniqSpots = size(uniqspots);
    
    % Loop over the spots again to eliminate redundant ones
    % by finding and keeping the max intensity one
    for j=1:numUniqSpots(1)
        instances = find(uniqulum(:,4) == uniqspots(j));
        peak = find(uniqulum(:,5) == max(uniqulum(instances,5)));
        spot = uniqulum(peak,:);
        uniqulus = [uniqulus
            spot(1,:)];
        
    end
    
    % Generate indices for the nucleus number
    number = cols(1:numUniqSpots) .* i;
    
    % Add the data to the list with the new numbers
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