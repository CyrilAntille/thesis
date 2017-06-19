classdef MainParameters
    %MAINPARAMETERS Parameters for all experiments (main_*.m scripts)
    
    properties
        % Array and medium parameters
        num_beams = 65;
        P = Parameters(65); % Argument must match num_beams
        % Note: If num_beams changed. must call copyP()
        medium_range = [35, 60] % mm
        % medium_range should start P.MinRadImage (35mm)
        NoMLA = 1
        
        % Speckle raw image creation parameters
        speckle_create = false;
        speckle_seed = 2;
        speckle_numPts = 10^6;
        speckle_save = true;
        % Speckle creation can take a long time. It is therefore
        % recommended to save it in a file.
        speckle_load = true;
%         speckle_file = '..\data\2_1_speckle_2_10-6.mat';
        speckle_file = '..\data\2_1_speckle_42_10-6.mat';
        % Speckle noise is added to the beamformed image if the speckle raw
        % image (speckle_raw) exists in the workspace (so most of the times
        % if speckle_create or speckle_load are true)
        
        % raw and DA(S) image creation parameters
        % Add a range and azimuth (in mm) for each point
        pts_range = [40];
        pts_azimuth = [0]; % If non-zero, no guaranty that at least
        % on beam hits the point perfectly (-> risk of scalloping loss!)
        pts_gain = 30; % dB over background average
        shift = Shift(ShiftType.RadialVar, 1/2, 1, 0); % Ref Shift.m
        shift_per_beam = false;
        % If true, the scatterer points are shifted after each beam.
        % This results in a single image per num_beam value.
        % If false, the scatterer points are shifted after each frame.
        % This results in shift.num_shifts images per num_beam value.
        
        % Beamforming parameters
        % Beamformer options: DAS, MV, MV-MB, IAA-MBSB, IAA-MBMB,
        % IAA-MBMB-Upsampled, MV-L, where L = subarray length ratio
%         methods_set = {'DAS','MV','IAA-MBSB','IAA-MBMB','IAA-MBMB-Upsampled'};
        methods_set = {'DAS','MV','IAA-MBSB','IAA-MBMB'};
        dl = 5/100; % Diagonal loading. MV only.
        sbl = 1/2; % Subarray length ratio (-> 1/2 = L/2). MV only.
        dofb = true; % Do forward-backward averaging. MV only.
        tavg = 2; % Number of samples used for temporal averaging
        % nsd = ceil(rad2deg(P.Tx.Theta(end))); % Number of subdimensions.
        % All but DAS and MV.
        nsd = 31; % Number of subdimensions for IAA subspace span
        % If nsd = 0, optimal nsd estimated from eigenvalues decomposition
        upsample_number = 500; % Number of received beams for interpolation
        % by phase shifting. Only IAA MB/MB Upsampled.
        
        interp_upsample = 0; % 0 = Disabled
        interp_method = 'spline';
        normalize_bfim = false;
        norm_variant = 1; % 1 = Standalone normalization, 2 = DAS-Based normalization,
        % 3 = First-shift based normalization (only if shift_per_beam = false)
        
        save_all_data = false; % Can take multiple GB of memory
        save_plots = false; % If false, will display them instead
        save_folder = '../output/'; % Ref createOutputDir()
        files_prefix = '';
        % Multiprocessing might become an issue if little RAM space available.
        % For disable to work, must disable automatic creation of parallel pool
        % for parfor (in parallel preferences).
        disable_multiprocess = false;
    end
    
    methods
        function P = copyP(obj, NTheta, NoMLA)
            if nargin == 2
                NoMLA = 1;
            end
            P = copyStruct(obj.P);
            P.Tx.NTheta = NTheta;
            P.Tx.SinTheta = linspace(-P.Tx.SinThMax,P.Tx.SinThMax,P.Tx.NTheta);
            P.Tx.Theta = asin(P.Tx.SinTheta);
            P.Rx.NoMLA = NoMLA;
            P.Rx.DeltaSinTh = 0;
            nm = floor(P.Rx.NoMLA/2);
%             diffSinTh = (P.Tx.SinTheta(2) - P.Tx.SinTheta(1)) / (nm+2);
            diffSinTh = (P.Tx.SinTheta(2) - P.Tx.SinTheta(1)) / P.Rx.NoMLA;
            for n=1:nm
                P.Rx.DeltaSinTh = horzcat(P.Rx.DeltaSinTh, n*diffSinTh);
                if mod(P.Rx.NoMLA, 2) == 1
                    P.Rx.DeltaSinTh = horzcat(-n*diffSinTh, P.Rx.DeltaSinTh);
                end
            end
            P.Rx.SinTheta = zeros([1 P.Rx.NoMLA * P.Tx.NTheta]);
            for tx=1:P.Tx.NTheta
                for rx=1:P.Rx.NoMLA
                    P.Rx.SinTheta(P.Rx.NoMLA*(tx-1)+rx) = ...
                        P.Tx.SinTheta(tx) + P.Rx.DeltaSinTh(rx);
                end
            end
            nm_even = mod(P.Rx.NoMLA + 1, 2);% = 1 if NoMLA is even
            P.Rx.SinTheta = P.Rx.SinTheta(1+nm-nm_even:end-nm); % to keep within Tx range
%             P.Rx.SinTheta = linspace(-P.Tx.SinThMax, P.Tx.SinThMax, P.Rx.NoMLA*P.Tx.NTheta);
            P.Rx.Theta = asin(P.Rx.SinTheta);
        end
        function scanGridTheta = getScanGrid(obj, bf_method)
            scanGridTheta = obj.P.Rx.Theta;
            if strfind(bf_method, 'Upsampled')
                NTheta = obj.upsample_number;
                scanGridSin = linspace(obj.P.Rx.SinTheta(1),obj.P.Rx.SinTheta(end),NTheta);
                scanGridTheta = asin(scanGridSin);
            elseif strfind(bf_method, 'IAA-MBMB-') | strfind(bf_method, 'IAA-MBSB-')
                upsample = strsplit(bf_method, '-');
                upsample = str2double(upsample(end));
                NTheta = length(obj.P.Rx.SinTheta);
                scanGridSin = interp1(1:NTheta, obj.P.Rx.SinTheta, 1:1/upsample:NTheta);
                scanGridTheta = asin(scanGridSin);
            end
        end
        function output_file = outputFileName(obj, file_type)
            % Output is .png if is_plot=true, else .mat
            % Output file name: <prefix>_<shift_per_beam(0=false,1=true)>_
            % <num_beams>_<shift.type>_<shift.val>.mat
            % Example: 1_61_RadialVar_0.5.mat
            output_file = [obj.files_prefix num2str(obj.shift_per_beam,3)...
                '_' int2str(obj.num_beams) '_' char(obj.shift.type) '_' ...
                num2str(obj.shift.val)];
            if strcmp(file_type, 'png')
                output_file = strcat(obj.save_folder, 'png/', ...
                    output_file, '.png');
            elseif strcmp(file_type, 'fig')
                output_file = strcat(obj.save_folder, 'fig/', ...
                    output_file, '.fig');
            else
                output_file = strcat(obj.save_folder, output_file, '.mat');
            end
        end
        function obj = createOutputDir(obj)
            if ~(obj.save_all_data || obj.save_plots)
                return
            end
            output_folder = datestr(datetime('now'),...
                'yyyy-mm-dd_HH.MM.SS');
            mkdir(obj.save_folder, output_folder)
            obj.save_folder = strcat(obj.save_folder, output_folder, '/');
            save(strcat(obj.save_folder, 'MainP.mat'), 'obj')
            if obj.save_plots
                mkdir(obj.save_folder, 'png')
                mkdir(obj.save_folder, 'fig')
            end

            fid = fopen(strcat(obj.save_folder, 'MainP.txt'),'wt');
            fprintf(fid, '==============================\n');
            fields = fieldnames(obj);
            for f=1:length(fields)
                fname = fields{f};
                fval = obj.(fname);
                if isa(fval,'numeric')|| isa(fval,'logical')
                    fval = num2str(fval);
                    fprintf(fid, '%s\t\t%s\n', fname, fval);
                elseif isa(fval,'cell')
                    fprintf(fid, '%s\t\t%s', fname, fval{1});
                    for v=2:length(fval)
                        fprintf(fid, ', %s', fval{v});
                    end
                    fprintf(fid, ' \n');
                elseif isa(fval,'Shift') || isa(fval,'struct')
                    fprintf(fid, '%s\n', fname);
                    s_fields = fieldnames(fval);
                    for s=1:length(s_fields)
                        sname = s_fields{s};
                        sval = fval.(sname);
                        if isa(sval,'numeric')|| isa(sval,'logical')
                            sval = num2str(sval);
                        elseif isa(sval,'char')
                            sval = char(sval);
                        else
                            continue % complex fields not displayed
                        end
                        fprintf(fid, '\t%s\t\t%s\n', sname, sval);
                    end
                else
                    fprintf(fid, '%s\t\t%s\n', fname, fval);
                end
            end
            fprintf(fid, '==============================\n');
            fclose(fid);
        end
    end
end

