function [FFT] = WANDER_FFT_fourier_MAG(isubject,force,rootpath)
% based on WANDER_TFR_phaselocking_balanced_V2b 

if rootpath == 1
    fname_FFT_fourier = ['w:\WANDER\data\PFR\s' num2str(isubject) '_FFT_fourier_MAG.mat'];
else
    fname_FFT_fourier = ['/shared/projects/project_wander/WANDER/data/PFR/s' num2str(isubject) '_FFT_fourier_MAG.mat'];
end

if exist(fname_FFT_fourier,'file') && force ~= 1
    fprintf('Returning FFT fourier \n');
    load(fname_FFT_fourier);
else
    fprintf('FFT fourier not found, creating it now! \n');
    
    % load FFT EGG
    FFT_EGG = WANDER_FFT_EGG(isubject,0,rootpath,0);
    
    % load data
    data_MEG    = WANDER_ICA(isubject,0,rootpath,0);
    data_EGG    = WANDER_epoch_EGG(isubject,rootpath,0);
    
    for iblock = 1 : 4
        
        %% select all channels, only correct rejections
        cfg                 = [];
        cfg.channel         = 'MEG*1';
        cfg.trials          = find(data_MEG{iblock}.trialinfo(:,3) == 4);
        data_MEG{iblock}    = ft_selectdata(cfg,data_MEG{iblock});
        data_MEG{iblock}    = rmfield(data_MEG{iblock},'cfg');
        
        cfg.channel         = {'phase'};
        data_EGG{iblock}    = ft_selectdata(cfg,data_EGG{iblock});
        data_EGG{iblock}    = rmfield(data_EGG{iblock},'cfg');
        
        %% segment trial in second epochs
        start_sample = 0;
        end_sample   = 0;
        trl          = [];
        Fs           = data_MEG{iblock}.fsample;
        for itrial = 1 : size(data_MEG{iblock}.trial,2)
            disp(num2str(itrial));
            itime = 1; % start of segments, in seconds
            while (data_MEG{iblock}.sampleinfo(itrial,1) + (itime)*Fs) < (data_MEG{iblock}.sampleinfo(itrial,2) - 1*Fs) % remove last second
                start_sample    = round(data_MEG{iblock}.sampleinfo(itrial,1) + itime*Fs);
                end_sample      = round(data_MEG{iblock}.sampleinfo(itrial,1) + itime*Fs + 1*Fs)-1;
                trl             = [trl; start_sample end_sample 0 data_MEG{iblock}.trialinfo(itrial,:) iblock itrial itime];
                itime           = itime + 0.5;
            end
        end
        
        cfg                     = [];
        cfg.trl                 = trl;
        data_MEG_timebinned{iblock}   = ft_redefinetrial(cfg,data_MEG{iblock});
        data_EGG_timebinned{iblock}   = ft_redefinetrial(cfg,data_EGG{iblock});
        
        %% remove artefacts, as well as trials with > 20% rejection

        artdef      = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
        artdef_EGG  = WANDER_artefact_detection_EGG(isubject,0,rootpath,0); %%% bug in old code
        
        % adding trialnr to original data
        data_orig                       = data_MEG_timebinned{iblock};
        data_orig.trialinfo(:,1)        = 1:size(data_orig.trialinfo,1);
        
        % remove artefactual segments
        cfg                             = [];
        cfg.reject                      = 'complete'; 
        cfg.artfctdef                   = artdef{iblock};
        cfg.artfctdef.EGG.artifact      = artdef_EGG{iblock};
        data_MEG_timebinned{iblock}     = ft_rejectartifact(cfg,data_MEG_timebinned{iblock});
        data_EGG_timebinned{iblock}     = ft_rejectartifact(cfg,data_EGG_timebinned{iblock});

        % find trials that have < 20% data removed que to artefacts
        trials_to_keep = [];
        for itrial = unique(data_EGG_timebinned{iblock}.trialinfo(:,8))'
            if any(data_orig.trialinfo(:,8) == itrial)
                nr_segments_orig = sum(data_orig.trialinfo(:,8) == itrial);
                nr_segments_art  = sum(data_MEG_timebinned{iblock}.trialinfo(:,8) == itrial);
                if nr_segments_orig ~= nr_segments_art
                    fprintf('%0.2f percent removed from trialnr %d', 100-nr_segments_art/nr_segments_orig*100, itrial);
                    if 100-nr_segments_art/nr_segments_orig*100 > 20
                        fprintf(': removing whole trial\n');
                    else
                        fprintf('\n');
                        trials_to_keep = [trials_to_keep itrial];
                    end
                else
                    trials_to_keep = [trials_to_keep itrial];
                end
            end
        end
        clear data_orig
        
        % only keep trials that have < 20% data removed que to artefacts
        cfg = [];
        cfg.trials = find(ismember(data_MEG_timebinned{iblock}.trialinfo(:,8),trials_to_keep));
        data_MEG_timebinned{iblock} = ft_selectdata(cfg,data_MEG_timebinned{iblock});
        data_EGG_timebinned{iblock} = ft_selectdata(cfg,data_EGG_timebinned{iblock});

        %% EGG phase binning
        
        % add phase of EGG to trialinfo (phase at middle of epoch)
        for itrial = 1 : size(data_MEG_timebinned{iblock}.trial,2)
            data_MEG_timebinned{iblock}.trialinfo(itrial,10) = data_EGG_timebinned{iblock}.trial{itrial}(1000);
        end
        
        % convert and enter EGG phase bins
        phase_nrbins    = 18;                                              % The number of phase bins
        phase_bins      = -pi:2*pi/phase_nrbins:pi;                        % The extreme phases of each bin
        phase_axis      = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;   % The midpoint of each phase bin
        
        % bin phase, check usage with e.g.: [count, edges, index] = histcounts([1 2 3 1 2 3 6 7 8 4 2 1],[1 2 3]);
        [~,~,data_MEG_timebinned{iblock}.trialinfo(:,11)] = histcounts(data_MEG_timebinned{iblock}.trialinfo(:,10),phase_bins);                                         
        
    end
    
    data_MEG_timebinned_append = ft_appenddata([],data_MEG_timebinned{:});
    clear data_MEG_timebinned data_MEG data_EGG data_MEG
    
    % give unique trialnumber
    itrial = 1;
    data_MEG_timebinned_append.trialinfo(1,12) = itrial;  
    for i = 2 : size(data_MEG_timebinned_append.trialinfo,1)
        if data_MEG_timebinned_append.trialinfo(i,6) ~= data_MEG_timebinned_append.trialinfo(i-1,6)
            itrial = itrial + 1;
        end
        data_MEG_timebinned_append.trialinfo(i,12) = itrial;
    end
    
    % do FFT on second windows
    cfg                 = [];
    cfg.output          = 'fourier';
    cfg.method          = 'mtmfft';
    cfg.keeptrials      = 'yes';
    cfg.taper           = 'hanning';
    cfg.foi             = [1:30];
    FFT                 = ft_freqanalysis(cfg, data_MEG_timebinned_append);
    clear data_MEG_timebinned_append

    % save data
    fprintf('Saving data. This can take a while.\n');
  
    save(fname_FFT_fourier,'FFT','-v7.3');
end