classdef MainParameters
    %MAINPARAMETERS Parameters for all experiments (main_*.m scripts)
    
    properties
        % Array and medium parameters
        P = Parameters();
        
        % Speckle raw image creation parameters
        speckle_create = false;
        speckle_seed = 2;
        speckle_numPts = 10^6;
        speckle_save = true;
        % Speckle creation can take a long time. It is therefore
        % recommended to save it in a file.
        speckle_load = true;
        speckle_file = '..\data\2_1_speckle_2_10-6.mat';
        % Speckle noise is added to the beamformed image if the speckle raw
        % image (speckle_raw) exists in the workspace (so most of the times
        % if speckle_create or speckle_load are true
        
        % raw and DA(S) image creation parameters
        pts_range = [40]; % Add a range (in mm) for each point
        shift = Shift(ShiftType.RadialVar, 1/2, 3, 0); % Ref Shift.m
        num_beams = 101; % can be a single value or list of values
        shift_per_beam = false;
        % If true, the scatterer points are shifted after each beam.
        % This results in a single image per num_beam value.
        % If false, the scatterer points are shifted after each frame.
        % This results in shift.num_shifts images per num_beam value.
        
        % Beamforming parameters
        % Beamformer options: DAS, MV, MV-MB, IAA-MBSB, IAA-MBMB, IAA-MBMB-Upsampled
        methods_set = {'DAS','MV','IAA-MBSB','IAA-MBMB','IAA-MBMB-Upsampled'};
        dl = 5/100; % Diagonal loading. All but DAS.
        sbl = 1/2; % Subarray length ratio (-> 1/2 = L/2). MV only.
        dofb = true; % Do forward-backward averaging. MV only.
        % nsd = ceil(rad2deg(P.Tx.Theta(end))); % Number of subdimensions.
        % All but DAS and MV.
        nsd = 30; % Number of subdimensions. All but DAS and MV.
        upsample_number = 500; % Number of received beams for interpolation
        % by phase shifting. Only IAA MB/MB Upsampled.
        
        save_all_data = false; % Can take multiple GB of memory
        save_plots = false; % If false, will display them instead
        save_folder = '../output/'; % Ref createOutputDir()
        % Multiprocessing might become an issue if little RAM space available.
        % For disable to work, must disable automatic creation of parallel pool
        % for parfor (in parallel preferences).
        disable_multiprocess = false;
    end
    
    methods
        function P = copyP(obj, NTheta)
            P = copyStruct(obj.P);
            P.Tx.NTheta = NTheta;
            P.Tx.SinTheta = linspace(-P.Tx.SinThMax,P.Tx.SinThMax,P.Tx.NTheta);
            P.Tx.Theta = asin(P.Tx.SinTheta);
        end
        function output_file = outputFileName(obj, used_speckle)
            % Output file name: <lowest_num_beams>_<highest_num_beams>_
            % <shift_per_beam(0=false,1=true)>_<speckle_name>.mat
            % Example: 61_211_0_speckle2.mat
            speckle_name = '_noSpeckle';
            if used_speckle
                speckle_name = ['_speckle' num2str(obj.P.Seed)];
            end
            output_file = [num2str(obj.num_beams(1)) '_' ...
                num2str(obj.num_beams(end)) '_' ...
                num2str(obj.shift_per_beam) speckle_name '.mat'];
            output_file = strcat(obj.save_folder, output_file);
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
                elseif isa(fval,'Shift')
                    fprintf(fid, '%s\n', fname);
                    s_fields = fieldnames(fval);
                    for s=1:length(s_fields)
                        sname = s_fields{s};
                        sval = fval.(sname);
                        if isa(sval,'numeric')|| isa(sval,'logical')
                            sval = num2str(sval);
                        else
                            sval = char(sval);
                        end
                        fprintf(fid, '\t%s\t\t%s\n', sname, sval);
                    end
                end
            end
            fprintf(fid, '==============================\n');
            fclose(fid);
        end
    end
end
