function [ans] = run_2csl() 

%% Nathaniel I. Krefman 2014-05-05

% Summary:

% Script for calling two_color_subnuclear_localization.m

% The variables in this script need to be updated every time 
% you run the code for a new image stack.

% Refer to the main code (two_color_subnuclear_localization.m) 
% for dependencies.

% Most of the code in two_color_subnuclear_localization.m does 
% not need to be updated at runtime, except to optimize for
% localization of a new spot type, or when using new microscope 
% or illumination settings, for example.

% Set image range prior to running code using MATLAB 
% (adjRangeAndBitDepth.m) or ImageJ.

%% Code

% Specify the path to the folder containing your current image stack
source_directory = '/Users/Username/folderA/';

% Provide complete paths to the image stacks for the spots (e.g. Mtw1-3xGFP,
% GAL locus tagged with GFP, Spc42-GFP, etc) and the nuclear boundary 
% (usually NSG1-mCherry).
imG = strcat(source_directory, '01L_spot.tif');
imR = strcat(source_directory, '01L_nucleus.tif');

% Specify the path to a folder for output images
% Program will overwrite any existing files with the same name.
% Below is an example, if the full path were
% '/Users/Username/folderB/subfolder1/subfolder2/', for example.
% Name should begin and end with a slash.
% The depth of the path is arbitrary.
save_directory = '/Users/Username/folderB/';

% Specify subfolder for saving images.
% Directory must already exist. This will not create it.
% This part of the path changes for each image stack.
% Subfolder1 might correspond to a strain or condition and 
% subfolder2 might be data from one image stack of several.
dir = strcat(save_directory,'subfolder1/subfolder2/');

% Give the stack a number in the series to distinguish it
% (assuming you have multiple stacks per strain/condition).
stack = 1;

% Execute the code
ans = two_color_subnuclear_localization(dir,stack,imG,imR);

