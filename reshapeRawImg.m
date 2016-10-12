function [ reshaped_im ] = reshapeRawImg( mainP, raw_im )
%RESHAPERAWIMG Reshapes raw_im according to mainP attributes (medium_range,
%P.c, P.fs,...). Returns a new structure containing the modified image.

samples_per_mm = 2*mainP.P.fs / mainP.P.DesFactor / mainP.P.c / 1000; % S/mm
med_range_mm = mainP.medium_range(1):1/samples_per_mm:mainP.medium_range(2);
med_range_m = med_range_mm / 1000;

range_sample_start = med_range_mm(1) * samples_per_mm;
range_time_start = range_sample_start / (mainP.P.fs / mainP.P.DesFactor);
% Factor 2 because of two-way beamforming (-> 1mm range == 2mm prop. dist.)

zero_img = zeros([length(med_range_m), ...
    mainP.P.Tx.no_elements * mainP.P.Rx.no_elements]);

raw_sample_start = floor(raw_im.times * mainP.P.fs / mainP.P.DesFactor);
start_diff = raw_sample_start - range_sample_start;
if start_diff < 0
    % -> raw_im start out of boundary
    raw_length = size(raw_im.image,1) + start_diff;
    raw_length = min(raw_length, length(med_range_m));
    if raw_length > 0
        zero_img(1:raw_length,:) = ...
            raw_im.image(-start_diff+1:-start_diff+raw_length,:);
    end
else
    raw_length = min(size(raw_im.image,1), length(med_range_m) - start_diff);
    if raw_length > 0
        zero_img(start_diff+1:start_diff+raw_length,:) = ...
            raw_im.image(1:raw_length,:);
    end
end

reshaped_im = struct; reshaped_im.image = zero_img;
reshaped_im.times = range_time_start;
reshaped_im.LConvs = raw_im.LConvs;
end

