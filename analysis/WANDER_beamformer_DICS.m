function source = WANDER_beamformer_DICS(isubject,rootpath,force)

if rootpath == 1
    fname_source = ['w:\WANDER\data\beamformer\TFR_source_s' num2str(isubject) '_DICS.mat'];
else
    fname_source = ['/shared/projects/project_wander/WANDER/data/beamformer/TFR_source_s' num2str(isubject) '_DICS.mat'];
end

% This is to create filter only - using artefact rejected data
if exist(fname_source,'file') && force == 0
    load(fname_source);
else
    disp('No beamformer found - calculating');
    
    [common_filter leadfield] = WANDER_common_filter_DICS(isubject,rootpath,0);
    
    % load data
    data_MEG    = WANDER_ICA(isubject,0,rootpath,0);
    data_EGG    = WANDER_epoch_EGG(isubject,rootpath,0);
    
    for iblock = 1 : 4
        % select only gradiometers & correct rejections
        cfg                 = [];
        cfg.channel         = 'MEGGRAD';
        cfg.trials          = find(data_MEG{iblock}.trialinfo(:,3) == 4);
        data_MEG{iblock}    = ft_selectdata(cfg,data_MEG{iblock});
        
        cfg.channel         = {'phase'};
%         cfg.channel         = 'BIO004_phase'; % only for testing on computer
        data_EGG{iblock}    = ft_selectdata(cfg,data_EGG{iblock});
        
        % segment trial in second epochs
        start_sample = 0;
        end_sample   = 0;
        trl          = [];
        Fs = 1000;
        for itrial = 1 : size(data_MEG{iblock}.trial,2)
            disp(num2str(itrial));
            itime = 1; % start of segments, in seconds
            while (data_MEG{iblock}.sampleinfo(itrial,1) + (itime)*Fs) < (data_MEG{iblock}.sampleinfo(itrial,2) - 1*Fs) % remove last second
                start_sample    = round(data_MEG{iblock}.sampleinfo(itrial,1) + itime*Fs);
                end_sample      = round(data_MEG{iblock}.sampleinfo(itrial,1) + itime*Fs + 1*Fs);
                trl             = [trl; start_sample end_sample 0 data_MEG{iblock}.trialinfo(itrial,:) iblock itrial itime];
                itime           = itime + 1;
            end
        end
        
        cfg                     = [];
        cfg.trl                 = trl;
        data_MEG_timebinned{iblock}   = ft_redefinetrial(cfg,data_MEG{iblock});
        data_EGG_timebinned{iblock}   = ft_redefinetrial(cfg,data_EGG{iblock});
        
        % remove epochs with artefacts
        artdef = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
        
        cfg = [];
        cfg.reject                      = 'complete'; % make bugreport, does not work with cfg.artfctdef.reject
        cfg.artfctdef.minaccepttim      = 0;
        cfg.artfctdef                   = artdef{iblock};
        data_MEG_timebinned{iblock}     = ft_rejectartifact(cfg,data_MEG_timebinned{iblock});
        data_EGG_timebinned{iblock}     = ft_rejectartifact(cfg,data_EGG_timebinned{iblock});        

        % EGG artefact removal
        temp = WANDER_artefact_detection_EGG(isubject,0,rootpath,0); %%% bug in old code
        
        cfg = [];
        cfg.reject                      = 'complete'; % make bugreport, does not work with cfg.artfctdef.reject
        cfg.artfctdef.minaccepttim      = 0;
        cfg.artfctdef.EGG.artifact      = temp{iblock}; %%% bug in old code
        data_MEG_timebinned{iblock}     = ft_rejectartifact(cfg,data_MEG_timebinned{iblock});
        data_EGG_timebinned{iblock}     = ft_rejectartifact(cfg,data_EGG_timebinned{iblock});             

        % add phase of EGG to trialinfo (phase at middle of epoch)
        for itrial = 1 : size(data_MEG_timebinned{iblock}.trial,2)
            data_MEG_timebinned{iblock}.trialinfo(itrial,10) = data_EGG_timebinned{iblock}.trial{itrial}(500);
        end
        
        % prepare EGG phase bins
        phase_nrbins    = 18;                                              % The number of phase bins
        phase_bins      = -pi:2*pi/phase_nrbins:pi;                        % The extreme phases of each bin
        phase_axis      = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;   % The midpoint of each phase bin
        [~,~,data_MEG_timebinned{iblock}.trialinfo(:,11)] = histcounts(data_MEG_timebinned{iblock}.trialinfo(:,10),phase_bins);                                         % bin phase, check usage e.g. with: [count, edges, index] = histcounts([1 2 3 1 2 3 6 7 8 4 2 1],[1 2 3]);
        
    end
    data_MEG_timebinned_append = ft_appenddata([],data_MEG_timebinned{:});

    clear data_MEG_timebinned data_MEG data_EGG data_MEG
    
    %     do FFT for beamformer
    cfg                 = [];
    cfg.output          = 'powandcsd';
%     cfg.channel         = 'MEGGRAD';
    cfg.method          = 'mtmfft';
    cfg.keeptrials      = 'yes';
    cfg.taper           = 'hanning';
    cfg.pad             = 'nextpow2';
    cfg.foi             = [10 11];
    FFT                 = ft_freqanalysis(cfg, data_MEG_timebinned_append);
    
    % load headmodel
    [mri_segmented, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(isubject,rootpath,0);
    
    % load sensor information
    [dataset_task, ~] = WANDER_subjectinfo;
    hdr = ft_read_header(dataset_task{isubject,1});
    
    % beamforming
    for itrial = 1 : size(FFT.trialinfo,1)
        
        cfg = [];
        cfg.trials = itrial;
        FFT_sel = ft_selectdata(cfg,FFT);
        
        cfg              = [];
        cfg.channel      = {'MEG*2', 'MEG*3'};
        cfg.channel      = 'all';
        cfg.grad         = hdr.grad;
        
        cfg.method       = 'dics';
        cfg.grid         = leadfield;
        cfg.grid.filter  = common_filter;
        cfg.headmodel    = headmodel;
        
        cfg.keeptrials   = 'yes';
        cfg.rawtrial     = 'no';
        
        cfg.frequency    = 10;
        source10{itrial} = ft_sourceanalysis(cfg, FFT_sel);
        
        
        cfg.frequency    = 11;
        source11{itrial} = ft_sourceanalysis(cfg, FFT_sel);
    end
    
    source = source10{1};
    source = rmfield(source,{'avg','method','cfg'});
    source.dimord = 'pos_rpt';
    source.freq = (source10{1}.freq + source11{1}.freq) / 2;
    
    for itrial = 1 : size(source10,2)
        source.pow(:,itrial) = (source10{itrial}.avg.pow + source11{itrial}.avg.pow) ./ 2;
        source.trialinfo(itrial,:) = source10{itrial}.trialinfo;
    end
    
    save(fname_source,'source','-v7.3');

end





