function source_TFR = WANDER_beamformer(isubject,rootpath,force)

if rootpath == 1
    fname_source = ['w:\WANDER\data\beamformer\TFR_source_s' num2str(isubject) '.mat'];
else
    fname_source = ['/shared/projects/project_wander/WANDER/data/beamformer/TFR_source_s' num2str(isubject) '.mat'];
end

% This is to create filter only - using artefact rejected data
if exist(fname_source,'file') && force == 0
    load(fname_source);
%     source_TFR = source;
else
    disp('No beamformer found - calculating');
    
     [common_filter leadfield] = WANDER_common_filter(isubject,rootpath,0);

    % load data
    data       = WANDER_ICA(isubject,0,rootpath,0);
    for ipart = 1 : 4
        % select only gradiometers
        cfg                     = [];
%         cfg.channel             = 'MEGGRAD';
        data{ipart}             = ft_selectdata(cfg,data{ipart});
    end
    
    % concatinate data over blocks
    data_append             = ft_appenddata([],data{:});
    clear data
     
    % initialize source-TFR data structure
    source_TFR              = [];
    source_TFR.pos          = leadfield.pos;
    source_TFR.dim          = leadfield.dim;
    source_TFR.inside       = leadfield.inside;
    source_TFR.trialdimord  = '{pos}_rpt_freq_time';

    % now project data through filter   
    for ivoxel = 1:size(common_filter,1)
        if ~isempty(common_filter{ivoxel})
            disp(['progress: ' round(num2str(ivoxel/size(common_filter,1)*100)) '%'])
            clear source_timecourse
            
            % trick fieldtrip is assuming regular timecourse datastructure
            source_timecourse           = data_append;
            source_timecourse           = rmfield(source_timecourse,'trial');
            source_timecourse           = rmfield(source_timecourse,'label');
            source_timecourse           = rmfield(source_timecourse,'cfg');
            source_timecourse.label{1}  = 'voxel';

            % get source timecourse, and do singular vector decomposition
            for itrial = 1 : size(data_append.trial,2) %numl 
                [source_timecourse.trial{itrial}, ~]= svdfft(common_filter{ivoxel} * data_append.trial{itrial}, 1);
            end
            
            % do mtmconvol on faux data structure
            cfg                 = [];
            cfg.pad             = 'nextpow2';
            cfg.channel         = 'all';
            cfg.method          = 'mtmconvol';
            cfg.keeptrials      = 'yes';
            cfg.foi             = [10,11];
            cfg.taper           = 'hanning';
            cfg.toi             = 0:0.05:30;
            cfg.t_ftimwin       = ones(size(cfg.foi));
            TFR                 = ft_freqanalysis(cfg, source_timecourse);
            
            % add to source TFR
            source_TFR.trial{ivoxel} = squeeze(permute(TFR.powspctrm,[1,3,4,2])); % reshape would be better
        end
    end
    save(fname_source,'source_TFR','-v7.3');

end





