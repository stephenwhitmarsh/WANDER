function [data_MEG_stim_avg, data_MEG_stim_avg_high, data_MEG_stim_avg_low, data_MEG_stim_bin, data_MEG_stim_bin_high, data_MEG_stim_bin_low, data_MEG_stim_MI, data_MEG_stim_MI_high, data_MEG_stim_MI_low] = WANDER_timelock_SSEF(isubject,rootpath,force)
          
if rootpath == 1
    fname_SSEF_phase = ['w:\WANDER\data\SSEF\SSEF_phase_' num2str(isubject) '.mat'];
else
    fname_SSEF_phase = ['/shared/projects/project_wander/WANDER/data/SSEF/SSEF_phase_' num2str(isubject) '.mat'];
end

if exist(fname_SSEF_phase,'file') && force ~= 1
    fprintf('Returning SSEF phasebinned data\n');
    load(fname_SSEF_phase);
else
    fprintf('SSEF phasebinned data not found, creating it now! \n');
    addpath('D:/analysis/WANDER/scripts/');
    
    % load data
    data_MEG = WANDER_epoch_MEG(isubject,0,rootpath,0);
    data_EGG = WANDER_epoch_EGG(isubject,rootpath,0);    
    
    for iblock = 1 : 4
        
        % select all channels, only correct rejections
        cfg                 = [];
        cfg.channel         = 'MEG';
        cfg.trials          = find(data_MEG{iblock}.trialinfo(:,3) == 4);
        data_MEG{iblock}    = ft_selectdata(cfg,data_MEG{iblock});
        data_MEG{iblock}    = rmfield(data_MEG{iblock},'cfg');

        cfg.channel         = {'phase'};        
%         cfg.channel         = {'phase','BIO004_phase'};
        data_EGG{iblock}    = ft_selectdata(cfg,data_EGG{iblock});
        data_EGG{iblock}    = rmfield(data_EGG{iblock},'cfg');
        
        % change timecourse according to latency of stimulation
        for itrial = size(data_MEG{iblock}.time,2) 
            data_MEG{iblock}.time{itrial} = data_MEG{iblock}.time{itrial} + 0.033;
        end
        
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = [11 21];
        data_MEG{iblock} = ft_preprocessing(cfg,data_MEG{iblock});
        
        % segment trial in stimulus epochs
        start_sample = 0;
        end_sample   = 0;
        trl          = [];
        for itrial = 1 : size(data_MEG{iblock}.trial,2)
            disp(num2str(itrial));
            onsetsample = find(data_MEG{iblock}.time{itrial} > 1,1,'first'); % start after one second
%             40; %2500/(1000/16) - could also just base it on time, i.e. 1500+1000ms
            istim = 1;
            while round(data_MEG{iblock}.sampleinfo(itrial,1) + onsetsample + (istim-1)*1000/16 + 63) < (data_MEG{iblock}.sampleinfo(itrial,2) - 1000) % stop before last second
                start_sample    = round(data_MEG{iblock}.sampleinfo(itrial,1) + onsetsample + (istim-1)*1000/16);
                end_sample      = round(data_MEG{iblock}.sampleinfo(itrial,1) + onsetsample + (istim-1)*1000/16 + 63);
                trl             = [trl; start_sample end_sample 0 data_MEG{iblock}.trialinfo(itrial,:) iblock itrial istim];
                istim           = istim + 1;
            end
        end
        
        cfg                     = [];
        cfg.trl                 = trl;
        data_MEG_stim{iblock}   = ft_redefinetrial(cfg,data_MEG{iblock});
        data_EGG_stim{iblock}   = ft_redefinetrial(cfg,data_EGG{iblock});
        
        % add phase of EGG to trialinfo (phase at middle of epoch)
        for itrial = 1 : size(data_MEG_stim{iblock}.trial,2)
            data_MEG_stim{iblock}.trialinfo(itrial,10) = data_EGG_stim{iblock}.trial{itrial}(31);
        end
        
        % prepare EGG phase bins
        phase_nrbins    = 18;                                              % The number of phase bins
        phase_bins      = -pi:2*pi/phase_nrbins:pi;                        % The extreme phases of each bin
        phase_axis      = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;   % The midpoint of each phase bin
        [~,~,data_MEG_stim{iblock}.trialinfo(:,11)] = histcounts(data_MEG_stim{iblock}.trialinfo(:,10),phase_bins);                                         % bin phase, check usage e.g. with: [count, edges, index] = histcounts([1 2 3 1 2 3 6 7 8 4 2 1],[1 2 3]);
        
        % remove epochs with artefacts
        artdef = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
        
        cfg = [];
        cfg.reject                      = 'complete'; % make bugreport, does not work with cfg.artfctdef.reject
        cfg.artfctdef.minaccepttim      = 0;
        cfg.artfctdef                   = artdef{iblock};
        data_MEG_stim{iblock}           = ft_rejectartifact(cfg,data_MEG_stim{iblock});
        
        % EGG artefact removal
        temp = WANDER_artefact_detection_EGG(isubject,0,rootpath,0); %%% bug in old code
        
        cfg = [];
        cfg.reject                      = 'complete'; % make bugreport, does not work with cfg.artfctdef.reject
        cfg.artfctdef.minaccepttim      = 0;
        cfg.artfctdef.EGG.artifact      = temp{iblock}; %%% bug in old code
        data_MEG_stim{iblock}           = ft_rejectartifact(cfg,data_MEG_stim{iblock});

    end
    
    % concatinate trials over blocks
    data_MEG_stim_append = ft_appenddata([],data_MEG_stim{:});
    clear data_MEG_stim
    
    % demean
    cfg = [];
    cfg.demean = 'yes';
    data_MEG_stim_append = ft_preprocessing(cfg,data_MEG_stim_append);
    
    % combine planar gradiometers
%     data_MEG_stim_append = ft_combineplanar([],data_MEG_stim_append);
    
    % calculate RMS per epoch
    disp('Calculating RMS...');
    clear trial_rms
    for itrial = 1 : size(data_MEG_stim_append.trial,2)
        trial_rms(itrial) = sqrt(mean(mean(data_MEG_stim_append.trial{itrial}(:,:) .^ 2)));
    end
    
    % retain epochs with correct rejections and RMS < 3SD
    thresh = mean(trial_rms) + std(trial_rms)*3;
    fprintf('removing %2.2f percent of data \n',size(find(trial_rms > thresh)) / size(trial_rms) * 100);
    
    cfg = [];
    cfg.trials = find(trial_rms < thresh & (data_MEG_stim_append.trialinfo(:,3) == 4)') ;
    data_MEG_stim_append = ft_selectdata(cfg,data_MEG_stim_append);
    
    % average over all correct rejections
    cfg                                 = [];
    data_MEG_stim_avg                   = ft_timelockanalysis([],data_MEG_stim_append);
    
    F                                   = ceil(2 * tiedrank(data_MEG_stim_append.trialinfo(:,2)) / length(data_MEG_stim_append.trialinfo(:,2)));
    rating_split                        = ones(size(F));
    rating_split(F==2)                  = 2;
    
    split_diff = size(find(rating_split==1),1) - size(find(rating_split==2),1);
    if split_diff < 0
        diff_bar = min(data_MEG_stim_append.trialinfo(rating_split == 2,2)); % rating bin between high/low split
        searchspace      = find(rating_split == 2 & data_MEG_stim_append.trialinfo(:,2) == diff_bar);
        searchspace_size = size(searchspace,1);
        start_pos        = randsample(1:searchspace_size-abs(split_diff),1);
        to_remove        = searchspace(start_pos:start_pos+abs(split_diff)-1);
        rating_split(to_remove) = 0;
    elseif split_diff > 0
        diff_bar = max(data_MEG_stim_append.trialinfo(rating_split == 1,2)); % rating bin between high/low split
        searchspace      = find(rating_split == 1 & data_MEG_stim_append.trialinfo(:,2) == diff_bar);
        searchspace_size = size(searchspace,1);
        start_pos        = randsample(1:searchspace_size-abs(split_diff),1);
        to_remove        = searchspace(start_pos:start_pos+abs(split_diff)-1);
        rating_split(to_remove) = 0;
    else
        disp('CRAZY! equal amount of stimuli per median split!');
    end
    
    cfg                     = [];
    cfg.trials              = find(rating_split == 1);
    data_MEG_stim_avg_low   = ft_timelockanalysis(cfg,data_MEG_stim_append);
    cfg.trials              = find(rating_split == 2);
    data_MEG_stim_avg_high  = ft_timelockanalysis(cfg,data_MEG_stim_append);
    
    % bins x rating
    clear data_rms*
    for ibin = 1 : phase_nrbins
        cfg                              = [];
        cfg.trials                       = find(data_MEG_stim_append.trialinfo(:,11) == ibin);
        data_MEG_stim_bin{ibin}          = ft_timelockanalysis(cfg,data_MEG_stim_append);
        cfg.trials                       = find(data_MEG_stim_append.trialinfo(:,11) == ibin & rating_split == 1);
        data_MEG_stim_bin_low{ibin}      = ft_timelockanalysis(cfg,data_MEG_stim_append);
        cfg.trials                       = find(data_MEG_stim_append.trialinfo(:,11) == ibin & rating_split == 2);
        data_MEG_stim_bin_high{ibin}     = ft_timelockanalysis(cfg,data_MEG_stim_append);
        
        data_rms(:,ibin)      = sqrt(mean(data_MEG_stim_bin{ibin}.avg .^2,2));
        data_rms_low(:,ibin)  = sqrt(mean(data_MEG_stim_bin_low{ibin}.avg .^2,2));
        data_rms_high(:,ibin) = sqrt(mean(data_MEG_stim_bin_high{ibin}.avg .^2,2));
    end
    
    % prepare MI datastructure
    data_MEG_stim_MI = [];
    data_MEG_stim_MI.label                  = data_MEG_stim_append.label;
    data_MEG_stim_MI.dimord                 = 'chan';
    
    data_MEG_stim_MI_low = [];
    data_MEG_stim_MI_low.label              = data_MEG_stim_append.label;
    data_MEG_stim_MI_low.dimord             = 'chan';
    
    data_MEG_stim_MI_high = [];
    data_MEG_stim_MI_high.label             = data_MEG_stim_append.label;
    data_MEG_stim_MI_high.dimord            = 'chan';
    
    % modulation index
    for ichan = 1 : size(data_MEG_stim_bin{ibin}.label,1)
        data_MEG_stim_MI.avg(ichan,1)       = (log(phase_nrbins)-(-sum((data_rms(ichan,:)       / sum(data_rms(ichan,:)))       .* log((data_rms(ichan,:)      /sum(data_rms(ichan,:)))))))      /log(phase_nrbins);
        data_MEG_stim_MI_high.avg(ichan,1)  = (log(phase_nrbins)-(-sum((data_rms_high(ichan,:)  / sum(data_rms_high(ichan,:)))  .* log((data_rms_high(ichan,:) /sum(data_rms_high(ichan,:))))))) /log(phase_nrbins);
        data_MEG_stim_MI_low.avg(ichan,1)   = (log(phase_nrbins)-(-sum((data_rms_low(ichan,:)   / sum(data_rms_low(ichan,:)))   .* log((data_rms_low(ichan,:)  /sum(data_rms_low(ichan,:)))))))  /log(phase_nrbins);
    end
    
%     save(fname_SSEF_phase, 'data_MEG_stim_avg', 'data_MEG_stim_avg_high', 'data_MEG_stim_avg_low', 'data_MEG_stim_bin', 'data_MEG_stim_bin_high', 'data_MEG_stim_bin_low', 'data_MEG_stim_MI', 'data_MEG_stim_MI_high', 'data_MEG_stim_MI_low','-v7.3');
    save(fname_SSEF_phase, 'data_MEG_stim_avg', 'data_MEG_stim_avg_high', 'data_MEG_stim_avg_low', 'data_MEG_stim_MI', 'data_MEG_stim_MI_high', 'data_MEG_stim_MI_low','data_MEG_stim_bin','data_MEG_stim_bin_high','data_MEG_stim_bin_low','-v7.3');

end
% 
% 
% cfg         = [];
% cfg.layout  = 'neuromag306mag';
% cfg.channel = {'MEG*1'};
% cfg.zlim    = 'absmax';
% figure; ft_topoplotER(cfg,data_MEG_stim_MI);
% 
% temp = data_MEG_stim_MI;
% data_MEG_stim_MI.rms = mean(data_rms,2);
% cfg.parameter = 'rms';
% figure; ft_topoplotER(cfg,data_MEG_stim_MI);
