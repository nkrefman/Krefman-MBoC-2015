function im = autogain_uint16(im)

% make sure the image is a double to begin with
im = double(im);

% set the minimum pixel to zero
im_min = min(im(:));
im = im - im_min;

% scale the intensities so the brightest pixel is 65535
im_max = max(im(:));
im = 65535*im/im_max;

% turn the image into a uint16 array
im = uint16(im);

end