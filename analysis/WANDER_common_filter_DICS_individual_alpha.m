function [common_filter, leadfield] = WANDER_common_filter_DICS_individual_alpha(isubject,rootpath,force)

if rootpath == 1
    fname_source = ['w:\WANDER\data\beamformer\common_filter_s' num2str(isubject) '_DICS_individual.mat'];
else
    fname_source = ['/shared/projects/project_wander/WANDER/data/beamformer/common_filter_s' num2str(isubject) '_DICS_individual.mat'];
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
        
%         % cut trials to size
%         for itrial = 1 : size(data{ipart}.trial,2)
%             fprintf('Cutting off beginning of trial %d block %d \n',itrial,ipart);
%             data{ipart}.trial{itrial} = data{ipart}.trial{itrial}(:,1501:end);
%             data{ipart}.time{itrial}  = data{ipart}.time{itrial}(1501:end);
%             fprintf('Cutting off end of trial %d block %d \n',itrial,ipart);
%             data{ipart}.trial{itrial} = data{ipart}.trial{itrial}(:,1:end-1000);
%             data{ipart}.time{itrial}  = data{ipart}.time{itrial}(1:end-1000);
%         end
%         
        % partial artefact rejection for common filter
        cfg = [];
        cfg.artfctdef               = artdef_MEG{ipart};
        cfg.artfctdef.reject        = 'partial';
        cfg.artfctdef.minaccepttim  = 3;
        data{ipart}                 = ft_rejectartifact(cfg,data{ipart});
    end
    
    % concatinate data over blocks
    data_append             = ft_appenddata([],data{:});
    clear data

    
%     % bandpass data
%     cfg = [];
%     cfg.bpfilter    = 'yes';
%     cfg.bpfreq      = [10 11];
%     cfg.bpfiltdir 	= 'twopass';
%     data_bpfilt     = ft_preprocessing(cfg,data_append);
%     
    % % bandpass data - using Craig Richter's (2017) code
    %
    % Nyq = data_append.fsample/2;
    % trans = .15;
    % data_bpfilt = data_append;
    % freq = [10,11];
    %
    % f = [0 (freq(1)/Nyq)-(trans*freq(1)/Nyq) freq(1)/Nyq freq(2)/Nyq freq(2)/Nyq+(trans*freq(2)/Nyq) 1]; m = [0 0 1 1 0 0];
    % b = fir2(150,f,m);
    %
    % for ii=1:size(data_bpfilt.trial{1},1)
    %     data_bpfilt.trial{1}(ii,:) = filtfilt(b,1,data_bpfilt.trial{1}(ii,:));
    % end
    %
    

    
    [dataset_task, dataset_rs] = WANDER_subjectinfo;
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
    cfg.normalize   = 'yes'; % might want to check effect without
    cfg.reducerank  = 2;
    leadfield       = ft_prepare_leadfield(cfg);
    
    % create CSD matrix for filter, using variable trial lenghts, but
    % controlled frequency smoothing
    
    cfg                 = [];

    % get individual alpha
    if isubject ~= 13
        [FFT_MEG_mag, ~] = WANDER_FFT_MEG_alpha(isubject,0,'cue',rootpath,1);
        cfg.foi = FFT_MEG_mag.max_freq;
    else
        cfg.foi = 11;
    end
    
    cfg.output          = 'powandcsd';
%     cfg.channel         = 'MEGGRAD';
    cfg.method          = 'mtmfft';
    cfg.keeptrials      = 'no';
    cfg.tapsmofrq       = 0.5;
    cfg.pad             = 'nextpow2';
%     cfg.foi             = 10.5;
    common_FFT          = ft_freqanalysis(cfg, data_append);
    
    % make filter over all (artefact-free) trials
    cfg = [];
    cfg.method              = 'dics';
    cfg.grad                = hdr.grad;
    cfg.grid                = leadfield;
    cfg.channel             = 'all';
    cfg.headmodel           = headmodel;
    cfg.dics.keepfilter     = 'yes';
    cfg.dics.lambda         = '5%';
    cfg.dics.reducerank     = 2;
    source_common           = ft_sourceanalysis(cfg, common_FFT);
    
    % only save filter and leadfield
    common_filter           = source_common.avg.filter;
    save(fname_source,'common_filter','leadfield');  
else
    disp('Loading common filter');
    load(fname_source);
end

