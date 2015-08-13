function writeTiffStack(tiffstack, outfilename)

% Based on suggestion at http://www.mathworks.com/matlabcentral/newsreader/view_thread/76992

% example: outfilename = 'TestTiffStack.tif'
% you can save a slice of tiffstack with
%   example: tiffstack = tiffstack(:,:,1:30)

imwrite(tiffstack(:,:,1), outfilename)
for k = 2:size(tiffstack,3)
    imwrite(tiffstack(:,:,k), outfilename, 'writemode', 'append');
end