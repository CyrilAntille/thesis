function [ bf_im ] = normalizeBFImage( bf_im, radius )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    minr = find(radius >= 43, 1);
    maxr = find(radius >= 47, 1);
    azl = size(bf_im, 2);
    az_range = floor(azl/3):ceil(2*azl/3); % To avoid edges samples
    bgd_zone = bf_im(minr:maxr, az_range);
    bgd_av = mean(bgd_zone(:));
    bf_im = bf_im / bgd_av;
end

