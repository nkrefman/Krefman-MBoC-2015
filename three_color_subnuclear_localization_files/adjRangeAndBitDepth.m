%% Nathaniel I. Krefman 2014-05-06

% Summary:

% Sets intensity ranges for a group of 16-bit tiff stacks stored 
% in the same folder. Reduces bit depth to 8-bit (reduces file
% size by 50%, but also loses reduces precision). Splits images
% into left half and right half (necessary to reduce file size 
% and processing time with cameras that aquire very large fields
% of view).

% GFP images should say 'GFP' in the filename, etc
% Directories should be organized like the following example:
%   '/Users/Username/infolder/25C/01 - WT/00min/01_GFP.tif'
%   '/Users/Username/infolder/25C/01 - WT/00min/01_mCherry.tif'
%   '/Users/Username/infolder/25C/01 - WT/00min/01_CFP.tif'

% Dependencies:
%   openTiffStack.m
%   writeTiffStack.m

%%

% Specify the path to a folder containing all the data
folder = '/Users/Username/infolder';
% Specify the path to an empty folder where the results will be 
% stored
savefolder = '/Users/Username/outfolder';

% Set max intensity thresholds.
% Must be higher than the max expected intensity of any real signal
% fluorescence. It's good to take a look at the raw images in ImageJ
% first to see what the intensities of the peaks are like. The final 
% threshold may be lower than the intensities of pixels near dead
% cells, since autofluorescence can be much brighter than signal from
% fluorescent proteins.
maxGFP = 9000;
maxmCherry = 600;
maxCFP = 600;

% Number of unique strains or conditions
numConditions = 1;
% Number of timepoints sampled
numTimepoints = 6;

for k=1:numConditions
    
    % tempstrains is a list of folders specific to each strain.
    % Each contains folders for each timepoint (listed below).
    tempsstrains = {'/25C/01 - WT'};
    tempstrain = tempsstrains(k);
    
    for j=1:numTimepoints
        
        % times is a list of subfolders (under tempstrains) specific 
        % to each timepoint. These folders should only contain image 
        % stacks.
        times = {'/000min','/030min','/060min','/090min','/120min','/150min'};
        tempstraintime = char(strcat(tempstrain, times(j)))
        
        % Get the list of images in the folder. Should work even if 
        % there are different numbers of images for each strain and 
        % timepoint.
        list = dir(strcat(folder,tempstraintime,'/*.TIF'));
        
        for i=1:size(list,1)
            fn = list(i).name;
            imNumber = fn(1:2);
            if findstr(fn,'GFP') > 0
                
                % Open the GFP image
                im = openTiffStack(strcat(folder, tempstraintime,'/',fn));
                
                % Scale the intensities,
                im = double(im);
                im = im/(maxGFP/255);
                %im = ceil(im); % This is closer to what ImageJ does.
                im = uint8(im);
                
                width = size(im,2)/2;
                
                %Save left half
                imleft = im(:,1:width,:);
                savefnleft = strcat(savefolder, tempstraintime, '/', imNumber,'L_GFP.tif');
                writeTiffStack(imleft, savefnleft);
                clear imleft
                
                % Save right half
                imright = im(:,width+1:end,:);
                savefnright = strcat(savefolder, tempstraintime, '/', imNumber,'R_GFP.tif');
                writeTiffStack(imright, savefnright);
                clear imright
                
                clear im
                
            end
            
            if findstr(fn,'mCherry') > 0
                
                % Open the mCherry image
                im = openTiffStack(strcat(folder, tempstraintime,'/',fn));
                
                % Scale the intensities.
                im = double(im);
                im = im/(maxmCherry/255);
                %im = ceil(im); % This is closer to what ImageJ does.
                im = uint8(im);
                
                width = size(im,2)/2;
                
                % Save left half.
                imleft = im(:,1:width,:);
                savefnleft = strcat(savefolder, tempstraintime, '/', imNumber,'L_mCherry.tif');
                writeTiffStack(imleft, savefnleft);
                clear imleft
                
                % Save right half.
                imright = im(:,width+1:end,:);
                savefnright = strcat(savefolder, tempstraintime, '/', imNumber,'R_mCherry.tif');
                writeTiffStack(imright, savefnright);
                clear imright
                
                clear im
                
            end
                
            if findstr(fn,'CFP') > 0
                
                % Open the CFP image
                im = openTiffStack(strcat(folder, tempstraintime,'/',fn));
                
                % Scale the intensities.
                im = double(im);
                im = im/(maxCFP/255);
                %im = ceil(im); % This is closer to what ImageJ does.
                im = uint8(im);
                
                width = size(im,2)/2;
                
                % Save left half.
                imleft = im(:,1:width,:);
                savefnleft = strcat(savefolder, tempstraintime, '/', imNumber,'L_CFP.tif');
                writeTiffStack(imleft, savefnleft);
                clear imleft
                
                % Save right half.
                imright = im(:,width+1:end,:);
                savefnright = strcat(savefolder, tempstraintime, '/', imNumber,'R_CFP.tif');
                writeTiffStack(imright, savefnright);
                clear imright
                
                clear im
                                
            end
        end
    end
end
