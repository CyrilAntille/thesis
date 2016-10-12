function [ bf_im ] = normalizeBFImage( bf_im, radius )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    minr = find(radius >= 43, 1);
    maxr = find(radius >= 47, 1);
    azl = size(bf_im, 2);
%     az_range = floor(azl/3):ceil(2*azl/3); 
    az_r1 = ceil(azl/8):floor(3*azl/8); % To avoid edges and center samples
    az_r2 = ceil(5*azl/8):floor(7*azl/8);
    bgd_z1 = bf_im(minr:maxr, az_r1);
    bgd_z2 = bf_im(minr:maxr, az_r2);
    bgd_av = mean([mean(bgd_z1(:)) mean(bgd_z2(:))]);
    bf_im = bf_im / bgd_av;
end

