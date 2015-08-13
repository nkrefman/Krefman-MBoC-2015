function im = medthresh(im)

%% Nathaniel I. Krefman 2014-05-05

% Thresholds using the median.

% Sets every pixel below the median pixel
% intensity to the median pixel intensity.

%%

im(find(im < (median(median(median(im)))))) = median(median(median(im)));
