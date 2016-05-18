function [ ref_shift, all_beams, all_phantoms, all_DA, all_BF ] ...
    = mergeData( files, save_to_file )
%MERGEDATA Summary of this function goes here
%   Detailed explanation goes here

all_beams = [];
all_phantoms = [];
all_DA = [];
all_BF = [];
ref_shift = 0;
for f=1:length(files)
    load(files{f})
    if ref_shift == 0
        all_beams = num_beams;
        all_phantoms = data_phantoms;
        all_DA = data_DA;
        if exist('data_BF', 'var')
            all_BF = data_BF;
        end
        ref_shift = shift_type;
    else
        if shift_type ~= ref_shift || shift_type.shift ~= ref_shift.shift ...
                || shift_type.num_shifts ~= ref_shift.num_shifts
            error('The variable shift_type differs between files!')
        end
        if num_beams(end) < all_beams(1)
            all_beams = [num_beams all_beams];
            all_phantoms = [data_phantoms all_phantoms];
            all_DA = [data_DA all_DA];
            if exist('data_BF', 'var')
                for m=1:length(all_BF)
                    all_BF{m} = [data_BF{m} all_BF{m}];
                end
            end
        elseif num_beams(1) > all_beams(end)
            all_beams = [all_beams num_beams];
            all_phantoms = [all_phantoms data_phantoms];
            all_DA = [all_DA data_DA];
            if exist('data_BF', 'var')
                for m=1:length(all_BF)
                    all_BF{m} = [all_BF{m} data_BF{m}];
                end
            end
        else
            error('num_beams overlap between files!')
        end
    end
    clearvars -except all_beams all_phantoms ...
        all_DA all_BF ref_shift files save_to_file
end
if save_to_file
    output_file = ['2_1_' num2str(all_beams(1)) '_' ...
        num2str(all_beams(end)) '.mat'];
    fprintf('\nNSaving DA(S) data to file.')
    save(output_file, 'shift_type', 'num_beams', 'data_phantoms', ...
        'data_DA', '-v7.3')
end
end

