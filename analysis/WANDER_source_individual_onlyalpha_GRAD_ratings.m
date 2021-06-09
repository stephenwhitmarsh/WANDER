function [FFT, source] = WANDER_source_individual_onlyalpha_GRAD_ratings(isubject,force,rootpath)

disp(['Working on Subjectnr. ' num2str(isubject)]);

if rootpath == 1
    fname_source_alpha = ['w:\WANDER\data\beamformer\s' num2str(isubject) '_source_individual_alpha_GRAD_ratings.mat'];
else
    fname_source_alpha = ['/shared/projects/project_wander/WANDER/data/beamformer/s' num2str(isubject) '_source_individual_alpha_GRAD_ratings.mat'];
end

if exist(fname_source_alpha,'file') && force ~= 1
    fprintf('Returning TFR-EGG phaselocking \n');
    load(fname_source_alpha);
else
    fprintf('TFR-EGG phaselocking not found, creating it now! \n');
    
    % load data
    data_MEG            = WANDER_ICA(isubject,0,rootpath,0);

    for iblock = 1 : 4
        
        % select all channels, only correct rejections
        cfg                             = [];
        cfg.channel                     = {'MEG*2','MEG*3'};
        cfg.trials                      = find(data_MEG{iblock}.trialinfo(:,3) == 4);
        data_MEG{iblock}                = ft_selectdata(cfg,data_MEG{iblock});
        data_MEG{iblock}                = rmfield(data_MEG{iblock},'cfg');
      
        % remove artefacts
        artdef                          = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
        cfg                             = [];
        cfg.reject                      = 'partial'; 
        cfg.artfctdef                   = artdef{iblock};
        cfg.artfctdef.minaccepttim      = 3;
        data_MEG_clean{iblock}          = ft_rejectartifact(cfg,data_MEG{iblock});
    end
    
    clear data_MEG
    data_MEG_append = ft_appenddata([],data_MEG_clean{:});
    clear data_MEG_clean 
    
    % load headmodel
    [~, headmodel, ~, ~, ~] = WANDER_grid(isubject,rootpath,0);
    
    % load sensor information
    [dataset_task, ~] = WANDER_subjectinfo(rootpath);
    hdr = ft_read_header(dataset_task{isubject,1});
    
    % load common filter
    [common_filter, leadfield] = WANDER_common_filter_DICS_individual_onlyalpha(isubject,rootpath,0);

    % get individual alpha
    if isubject ~= 13
        [FFT_MEG] = WANDER_FFT_MEG_alpha(isubject,0,'cue',rootpath,0);
        foi = FFT_MEG.max_freq;
    else
        foi = 11;
    end
    
    % recode ratings in set nr. of bins
    nrbins = 4;
    [ratings_binned, ~, ratings_stat] = bin_data(data_MEG_append.trialinfo(:,2),nrbins);
 
    for rating = 1:nrbins
        
        cfg = [];
        cfg.foi             = foi;
        cfg.output          = 'powandcsd';
        cfg.method          = 'mtmfft';
        cfg.channel         = {'MEG*2', 'MEG*3'};
        cfg.keeptrials      = 'no';
        cfg.tapsmofrq       = 0.5;
        cfg.pad             = 'nextpow2';
        cfg.trials          = find(ratings_binned == rating);
        FFT{rating}         = ft_freqanalysis(cfg, data_MEG_append);
        FFT{rating}         = rmfield(FFT{rating},'cfg');

        cfg                 = [];
        cfg.channel         = {'MEG*2', 'MEG*3'};
        cfg.grad            = hdr.grad;
        cfg.headmodel       = headmodel;
        cfg.keeptrials      = 'no';
        cfg.rawtrial        = 'no';
        cfg.method          = 'dics';
        cfg.grid            = leadfield;
        cfg.grid.filter     = common_filter;
        source{rating}      = ft_sourceanalysis(cfg, FFT{rating});
        source{rating}      = rmfield(source{rating},'cfg');
        source{rating}.ratings_stat = ratings_stat;
    end
    
    % save data
    fprintf('Saving data in %s. This can take a while.\n',fname_source_alpha);
    save(fname_source_alpha,'FFT','source','-v7.3');
end