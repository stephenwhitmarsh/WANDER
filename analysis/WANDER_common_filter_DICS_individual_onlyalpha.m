function [common_filter, leadfield] = WANDER_common_filter_DICS_individual_onlyalpha(isubject,rootpath,force)

if rootpath == 1
    fname_source = ['z:\WANDER\data\beamformer\common_filter_s' num2str(isubject) '_DICS_individual_onlyalpha.mat'];
else
    fname_source = ['/shared/projects/project_wander/WANDER/data/beamformer/common_filter_s' num2str(isubject) '_DICS_individual_onlyalpha.mat'];
end

% This is to create filter only - using artefact rejected data
if ~exist(fname_source,'file') || force == 1

    disp('No common filter found - calculating');
    
    % load headmodel
    [mri_segmented, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(isubject,rootpath,0);
    
    % load data
    data       = WANDER_ICA(isubject,0,rootpath,0);
    
    % load artefact definition
    artdef_MEG = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
    
    for ipart = 1 : 4
        % select only gradiometers
        cfg                     = [];
        cfg.channel             = 'MEGGRAD';
        data{ipart}             = ft_selectdata(cfg,data{ipart});
        
        % cut trials to size
        for itrial = 1 : size(data{ipart}.trial,2)
            
            fprintf('Cutting off beginning of trial %d block %d \n',itrial,ipart);
            data{ipart}.trial{itrial}               = data{ipart}.trial{itrial}(:,1501:end);
            data{ipart}.time{itrial}                = data{ipart}.time{itrial}(1501:end);
            data{ipart}.sampleinfo(itrial,1)        = data{ipart}.sampleinfo(itrial,1) + 1500;   
 
            fprintf('Cutting off end of trial %d block %d \n',itrial,ipart);
            data{ipart}.trial{itrial}               = data{ipart}.trial{itrial}(:,1:end-1000);
            data{ipart}.time{itrial}                = data{ipart}.time{itrial}(1:end-1000);
            data{ipart}.sampleinfo(itrial,2)        = data{ipart}.sampleinfo(itrial,2) - 1000;   
            
        end
        
        % partial artefact rejection for common filter
        cfg = [];
        cfg.artfctdef               = artdef_MEG{ipart};
        cfg.artfctdef.reject        = 'partial';
        cfg.artfctdef.minaccepttim  = 3;
        data{ipart}                 = ft_rejectartifact2(cfg,data{ipart});
    end
    
    % concatinate data over blocks
    data_append             = ft_appenddata([],data{:});
    clear data

    [dataset_task, dataset_rs] = WANDER_subjectinfo(rootpath);
    hdr = ft_read_header(dataset_task{isubject,1});
%     
%     cfg = [];
%     cfg.channel = 'MEGGRAD';
%     hdr = ft_selectdata(cfg,hdr);   % here is where it breaks
%     hdr.nChans = 204;
%     
    % compute leadfield once
    cfg             = [];
    cfg.grad        = hdr.grad;
    cfg.headmodel   = headmodel;
    cfg.grid        = subject_grid;
    cfg.channel     = {'MEGGRAD'};
    cfg.normalize   = 'yes'; % might want to check effect without ADDED AGAIN AFTER BLOB IN CENTER - 08-03-2019
    cfg.reducerank  = 2;
    leadfield       = ft_prepare_leadfield(cfg);
    
    % get individual alpha 
    %     if isubject ~= 13 % removed clause on 21-9-2018=
    %    [FFT_MEG] = WANDER_FFT_MEG_alpha(isubject,0,'cue',rootpath,0);
    %     cfg                 = [];
    %     cfg.foi             = FFT_MEG.max_freq;
    %     else
    %         cfg.foi = 11;
    %     end

    cfg                     = [];
    cfg.foi                 = 11; % NOTICED THIS WAS NOT DONE AND STILL STUCK ON SUBJECT SPECIFIC FREQ, SO CHANGED 08-03-2019
    cfg.output              = 'powandcsd';
%   cfg.channel             = 'MEGGRAD';    
    cfg.method              = 'mtmfft';
    cfg.keeptrials          = 'no';
    cfg.tapsmofrq           = 3;
    cfg.pad                 = 'nextpow2';
    common_FFT              = ft_freqanalysis(cfg, data_append);
    
    % make filter over all (artefact-free) trials
    cfg = [];
    cfg.method              = 'dics';
    cfg.grad                = hdr.grad;
    cfg.grid                = leadfield;
    cfg.channel             = 'all';
    cfg.headmodel           = headmodel;
    cfg.dics.keepfilter     = 'yes';
    cfg.dics.lambda         = '10%'; % try without regularizatin?
    cfg.dics.kappa          = rank(squeeze(data_append.trial{1}) * squeeze(data_append.trial{1})');
    cfg.dics.reducerank     = 2;
    source_common           = ft_sourceanalysis(cfg, common_FFT);
    
    % only save filter and leadfield
    common_filter           = source_common.avg.filter;
    save(fname_source,'common_filter','leadfield');  
else
    disp('Loading common filter');
    load(fname_source);
end

