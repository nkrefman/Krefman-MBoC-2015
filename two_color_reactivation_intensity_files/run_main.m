function run_main()  

%% Nathaniel I. Krefman 2014-05-05

% Summary:

% Script for calling main.m

% The variables in this script need to be updated every time 
% you run the code for a new image. We save a unique copy of this 
% script for each timepoint.

% Refer to main.m for dependencies.

% Most of the code in main.m does not need to be updated at 
% runtime, except to optimize for localization of a new spot 
% type, or when using new microscope or illumination settings, 
% for example.

% Set image range prior to running code using adjRangeAndBitDepth.m
% or ImageJ.

%%

% Specify the path to the folder containing your current image stack
source_directory = '/Users/Username/folderA';

% Provide complete paths to the image stacks for the green 
% channel (GFP-TUB1 and GAL-CEN-TetOs marked with TetR-GFP) and 
% Mad1-3xmCherry
imG = strcat(source_directory, '/00min/01L_GFP.tif');
imR = strcat(source_directory, '/00min/01L_mCherry.tif');

% Specify the path to a folder for output images
% Program will overwrite any existing files with the same name.
% Name should begin and end with a slash.
% The depth of the path is arbitrary.
save_directory = '/Users/Username/folderB';

% Specify subfolder for saving images. 
% Directory must already exist. This will not create it.
% '01L' refers to the first image stack, left half. '01R' would
% be the right half. These name's come from the output of the 
% script adjRangeAndBitDepth.m. For each new image, the number
% should be incremented.
dir = strcat(save_directory,'/01 - WT/00min/01L/');

% Give the stack a number in the series to distinguish it
% (assuming you have multiple stacks per timepoint).
stack = 1;

results = main(dir,stack,imG,imR); 
