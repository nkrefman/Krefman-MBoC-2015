%% Nathaniel I. Krefman 2014-05-05

% Summary:

% Sets intensity ranges for a group of 16-bit tiff stacks stored 
% in the same folder. Reduces bit depth to 8-bit (reduces file
% size by 50%, but also loses reduces precision). Splits images
% into left half and right half (necessary to reduce file size 
% and processing time with cameras that aquire very large fields
% of view).

% Dependencies:
%   openTiffStack.m
%   writeTiffStack.m

%%

% Specify the path to the folder containing only stacks of images
folder = '/Users/Username/infolder';

% Get the list of images in the folder
list = dir(strcat(folder, '/*.TIF'));

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
maxSpot = 1400;
maxNucleus = 600;

for i=1:size(list,1)
    fn = list(i).name;
    imNumber = fn(1:2);
    if findstr(fn,'spot') > 0
        % Open the spot image.
        im = openTiffStack(strcat(folder,'/',fn));
        
        % Scale the intensities.
        im = double(im);
        im = im/(maxSpot/255);
        %im = ceil(im); % This is closer to what ImageJ does.
        im = uint8(im);
        
        width = size(im,2)/2;

        % Save left half.
        imleft = im(:,1:width,:);
        savefnleft = strcat(savefolder, '/', imNumber,'L_spot.tif');
        writeTiffStack(imleft, savefnleft);
        clear imleft
        
        % Save right half.
        imright = im(:,width+1:end,:);
        savefnright = strcat(savefolder, '/', imNumber,'R_spot.tif');
        writeTiffStack(imright, savefnright);
        clear imright
        
        clear im
        
    end
    
    if findstr(fn,'nucleus') > 0
        % Open the mCherry image.
        im = openTiffStack(strcat(folder,'/',fn));
        
        % Scale the intensities.
        im = double(im);
        im = im/(maxNucleus/255);
        %im = ceil(im); % This is closer to what ImageJ does.
        im = uint8(im);
        
        width = size(im,2)/2;

        % Save left half.
        imleft = im(:,1:width,:);
        savefnleft = strcat(savefolder, '/', imNumber,'L_nucleus.tif');
        writeTiffStack(imleft, savefnleft);
        clear imleft
        
        % Save right half.
        imright = im(:,width+1:end,:);
        savefnright = strcat(savefolder, '/', imNumber,'R_nucleus.tif');
        writeTiffStack(imright, savefnright);
        clear imright
        
        clear im

    end
end

'Finished saving the images'
