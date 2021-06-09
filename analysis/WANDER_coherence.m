% t = 1:2000;
% ti = 2000-t;
% p = mod(ti,160);
% 
% figure; plot(p)

function ERF = WANDER_coherence(isubject,force,timing,rootpath)

if ~(strcmp(timing,'cue') || strcmp(timing,'probe'))
    fprintf('Use "cue" or "probe" as third argument.\n');
    return
end

if rootpath == 1
    fname_SSPAC = ['i:\analysis\WANDER\data\SSPAC\s' num2str(isubject) '.mat'];
else
    fname_SSPAC = ['/home/swhitmarsh/WANDER/data/SSPAC/s' num2str(isubject) '.mat'];
end

if exist(fname_SSPAC,'file') && force ~= 1
    fprintf('Steady-state PAC \n');
    load(fname_SSPAC);
else
    fprintf('Steady-state PAC not found, creating it now! \n');
    data_epoch_MEG = WANDER_ICA(isubject,0,rootpath,0);
    artdef = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
    
    % reject artifacts
    for ipart = 1 : 4
        cfg = [];
        cfg.reject                      = 'complete'; % make bugreport, does not work with cfg.artfctdef.reject
        %         cfg.artfctdef.minaccepttim      = 0;
        cfg.artfctdef                   = artdef{ipart};
        data_epoch_MEG{ipart}           = ft_rejectartifact(cfg,data_epoch_MEG{ipart});
        for itrial = 1 : size(data_epoch_MEG{ipart}.trial,2)
            fprintf('Cutting off end of trial %d block %d \n',itrial,ipart);
            data_epoch_MEG{ipart}.trial{itrial} = data_epoch_MEG{ipart}.trial{itrial}(:,1:end-1000);
            data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial}(1:end-1000);
            fprintf('Cutting off beginning of trial %d block %d \n',itrial,ipart);
            data_epoch_MEG{ipart}.trial{itrial} = data_epoch_MEG{ipart}.trial{itrial}(:,1501:end);
            data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial}(1501:end);
            fprintf('Shifting timecourse of trial %d block %d \n',itrial,ipart);
            data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial} - data_epoch_MEG{ipart}.time{itrial}(end) + 1.000;
        end
    end
    
    % concatinate trials
    data_epoch_MEG                      = ft_appenddata([],data_epoch_MEG{:});
    
    % reject all but correct rejections
    cfg                                 = [];
    cfg.trials                          = find(data_epoch_MEG.trialinfo(:,3) == 4);
    data_epoch_MEG                      = ft_selectdata(cfg,data_epoch_MEG);
    
    % maximum length trial
    maxlength = 0;
    minlength = 30000;
    
    for itrial = 1 : size(data_epoch_MEG.trial,2)
        if size(data_epoch_MEG.trial{itrial},2) > maxlength
            maxlength = size(data_epoch_MEG.trial{itrial},2);
        end
        if size(data_epoch_MEG.trial{itrial},2) < minlength
            minlength = size(data_epoch_MEG.trial{itrial},2);
        end
    end
    
    clear power_bincount power_binmean power_binmedian

    % high frequency TFR
    for itrial = 1 : size(data_epoch_MEG.trial,2)
        cfg             = [];
        cfg.pad         = 'nextpow2';
        cfg.channel     = 'all';
        cfg.method      = 'mtmconvol';
        cfg.keeptrials  = 'no';
        cfg.channel     = {'MEG'};
        cfg.foi         = 30:100;
        cfg.taper       = 'hanning';
        cfg.toi         = -(minlength/1000):0.01:0;
        cfg.t_ftimwin   = (1./cfg.foi)*4;
        cfg.trials      = itrial;
        
        trialTFR        = ft_freqanalysis(cfg, data_epoch_MEG);
        trialTFR        = ft_combineplanar([],trialTFR);
        
        % prepare phase bins
        phase_nrbins    = 18;                                              % The number of phase bins
        phase_bins      = -pi:2*pi/phase_nrbins:pi;                        % The extreme phases of each bin
        phase_axis      = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;   % The midpoint of each phase bin
        time            = trialTFR.time;
        phase           = mod(trialTFR.time,1/16)*16*2*pi-pi;
        [~,~,phase_indx] = histcounts(phase,phase_bins);                                         % bin phase, check usage e.g. with: [count, edges, index] = histcounts([1 2 3 1 2 3 6 7 8 4 2 1],[1 2 3]);
        phase_indx(phase_indx == 0) = NaN;
        for ifreq = 1 : size(trialTFR.powspctrm,2)                                                                  % loop through all frequencies
            fprintf('working on freq %d of trial %d \n',ifreq,itrial);
            for ichan_MEG = 1 : size(trialTFR.powspctrm,1)                                                          % loop through all sensors
                for ibin = unique(phase_indx(~isnan(phase_indx)))
                    power_bincount(itrial,ichan_MEG,ifreq,ibin)  = sum(~isnan(trialTFR.powspctrm(ichan_MEG,ifreq,phase_indx == ibin)),3);
                    power_binmean(itrial,ichan_MEG,ifreq,ibin)   = nanmean(trialTFR.powspctrm(ichan_MEG,ifreq,phase_indx == ibin),3);
                    power_binmedian(itrial,ichan_MEG,ifreq,ibin) = nanmedian(trialTFR.powspctrm(ichan_MEG,ifreq,phase_indx == ibin),3);
                end
            end
        end
    end
    
%     save('power','power*');
    
    % manually average over all trials, excluding nans
    PAC_binmean         = squeeze(nanmean(power_binmean,1));
    
    % make index according to median split based on ratings
    F = ceil(2 * tiedrank(data_epoch_MEG.trialinfo(:,2)) / length(data_epoch_MEG.trialinfo(:,2)));
    
    % before equalizing high/low bins
    rating_split        = ones(size(F));
    rating_split(F==2)  = 2;
    low_cnt             = size(find(rating_split == 1),1);
    high_cnt            = size(find(rating_split == 2),1);
    
    for irand = 1 : 24
        rating_split       = ones(size(F));
        rating_split(F==2) = 2;
        if high_cnt > low_cnt
            diff_bar                = min(data_epoch_MEG.trialinfo(rating_split == 2,2)); % rating bin between high/low split
            nr_to_remove            = high_cnt-low_cnt;
            discarted               = randsample(find(data_epoch_MEG.trialinfo(:,2) == diff_bar & rating_split == 2),nr_to_remove);
            rating_split(discarted) = 0;
        elseif high_cnt < low_cnt
            diff_bar                = max(data_epoch_MEG.trialinfo(rating_split == 1,2)); % rating bin between high/low split
            nr_to_remove            = low_cnt-high_cnt;
            discarted               = randsample(find(data_epoch_MEG.trialinfo(:,2) == diff_bar & rating_split == 1),nr_to_remove);
            rating_split(discarted) = 0;
        else
            discarted = [];
        end
        fprintf('I made %d low ratings, and %d high ratings in randomization %d \n',sum(rating_split==1),sum(rating_split==2),irand);
        
        temp_low(irand,:,:,:)  = squeeze(nanmean(power_binmean(rating_split==1,:,:,:),1));
        temp_high(irand,:,:,:) = squeeze(nanmean(power_binmean(rating_split==2,:,:,:),1));
    end
    
    % average over randomizations
    PAC_binmean_high    = squeeze(nanmean(temp_high,1));
    PAC_binmean_low     = squeeze(nanmean(temp_low,1));
    clear temp_low temp_high
    
    % calculate modulation index per frequency and channel and put it in a
    % FieldTrip timelock data structure
    MI        = [];
    MI.label  = trialTFR.label;
    MI.dimord = 'chan_time';
    MI.time   = 30:100;
    MI_high   = MI;
    MI_low    = MI;
    
    for ichan = 1 : size(PAC_binmean,1)
        fprintf('Calculating MI, channel %d of %d \n',ichan, size(PAC_binmean,1) );
        for ifreq = 1 : size(PAC_binmean,2)
            MI.avg(ichan,ifreq,:)      = (log(phase_nrbins)-(-sum((PAC_binmean(ichan,ifreq,:)      / sum(PAC_binmean(ichan,ifreq,:)))      .* log((PAC_binmean(ichan,ifreq,:)      /sum(PAC_binmean(ichan,ifreq,:)))))))      /log(phase_nrbins);
            MI_high.avg(ichan,ifreq,:) = (log(phase_nrbins)-(-sum((PAC_binmean_high(ichan,ifreq,:) / sum(PAC_binmean_high(ichan,ifreq,:))) .* log((PAC_binmean_high(ichan,ifreq,:) /sum(PAC_binmean_high(ichan,ifreq,:))))))) /log(phase_nrbins);
            MI_low.avg(ichan,ifreq,:)  = (log(phase_nrbins)-(-sum((PAC_binmean_low(ichan,ifreq,:)  / sum(PAC_binmean_low(ichan,ifreq,:)))  .* log((PAC_binmean_low(ichan,ifreq,:)  /sum(PAC_binmean_low(ichan,ifreq,:)))))))  /log(phase_nrbins);
        end
    end
    
    %normalize median split with average over all trials
    PAC_binmean_norm      = nan(size(PAC_binmean));
    PAC_binmean_high_norm = nan(size(PAC_binmean_high));
    PAC_binmean_low_norm  = nan(size(PAC_binmean_low));
    for ichan = 1:size(PAC_binmean,1)
        for ifreq = 1 : size(PAC_binmean,2)
            PAC_binmean_norm(ichan,ifreq,:)      = PAC_binmean(ichan,ifreq,:)      ./ nanmean(PAC_binmean(ichan,ifreq,:));
            PAC_binmean_high_norm(ichan,ifreq,:) = PAC_binmean_high(ichan,ifreq,:) ./ nanmean(PAC_binmean(ichan,ifreq,:));
            PAC_binmean_low_norm(ichan,ifreq,:)  = PAC_binmean_low(ichan,ifreq,:)  ./ nanmean(PAC_binmean(ichan,ifreq,:));
        end
    end
    
    MI_diff = MI;
    MI_diff.avg = (MI_low.avg - MI_high.avg) ./ (MI_low.avg + MI_high.avg) ;
    
    % put it in a FieldTrip frequency data format
    PAC                = [];
    PAC.label           = trialTFR.label;
    PAC.dimord          = 'chan_freq_time';
    PAC.freq            = 30:100;
    PAC.time            = 1:18;
    PAC_high            = PAC;
    PAC_low             = PAC;
    PAC.powspctrm       = PAC_binmean_norm;
    PAC_high            = PAC;
    PAC_low             = PAC;
    PAC_diff            = PAC;
    PAC_high.powspctrm  = PAC_binmean_high_norm;
    PAC_low.powspctrm   = PAC_binmean_low_norm;
    PAC_diff.powspctrm  = PAC_binmean_low_norm - PAC_binmean_high_norm;

    cfg                 = [];
    cfg.layout          = 'neuromag306cmb.lay';
    cfg.channel         = 'MEG*3';
    cfg.xlim            = [1 18];
    cfg.ylim            = [30 100];
    figure; ft_singleplotTFR(cfg, PAC)
    figure; ft_singleplotTFR(cfg, PAC_high)
    figure; ft_singleplotTFR(cfg, PAC_low)
    figure; ft_singleplotTFR(cfg, PAC_diff)
     
    cfg                 = [];
    cfg.layout          = 'neuromag306cmb.lay';
    cfg.channel         = 'MEG*3';
    cfg.xlim            = [30 100];
    cfg.ylim            = 'maxabs';
    figure; ft_singleplotER(cfg, MI)
    figure; ft_singleplotER(cfg, MI_diff)
    figure; ft_singleplotER(cfg, MI_low, MI_high)
    
    cfg.layout          = 'neuromag306mag.lay';
    cfg.channel         = 'MEG*1';
    cfg.xlim            = [30 100];
    cfg.ylim            = 'maxabs';
    figure; ft_singleplotER(cfg, MI)
    figure; ft_singleplotER(cfg, MI_diff)
    figure; ft_singleplotER(cfg, MI_low, MI_high)
    
save(fname_SSPAC,'PAC*','MI*');
%     
%     
%     %%
%     
%     
%     
%     
%     
%     
% end
% 
% 
% 
% cfg                  = [];
% cfg.layout           = 'neuromag306all.lay';
% figure; ft_singleplotER(cfg, TFR)
% 
% 
% 
% % cut up to shortest trial
% cfg                 = [];
% cfg.latency         = [-round(minlength/1000) 0];
% data_epoch_MEG_sel  = ft_selectdata(cfg,data_epoch_MEG);
% 
% % get max sensors at SSSEF freq
% cfg             = [];
% cfg.pad         = 'nextpow2';
% cfg.channel     = 'all';
% cfg.method      = 'mtmfft';
% cfg.foi         = 16;
% cfg.taper       = 'hanning';
% cfg.channel     = {'MEG*2','MEG*3'};
% cfg.t_ftimwin	= 1;
% %     cfg.toilim      = [-round(minlength/1000) 0];
% FFT             = ft_freqanalysis(cfg, data_epoch_MEG_sel);
% 
% [~,indx]        = sort(FFT.powspctrm,'descend');
% for isens = 1 : 9
%     maxsens{isens} = FFT.label{indx(isens)};
% end
% %
% 
% cfg            = [];
% cfg.output     = 'fourier';
% cfg.method     = 'mtmfft';
% cfg.pad        = 'nextpow2';
% cfg.foilim     = [5 20];
% cfg.tapsmofrq  = 0.1;
% cfg.keeptrials = 'yes';
% cfg.channel    = {'MEG*',maxsens{1}; 'MEG*',maxsens{2}; 'MEG*',maxsens{3}; 'MEG*',maxsens{4}; 'MEG*',maxsens{5}; 'MEG*',maxsens{6}};
% freqfourier    = ft_freqanalysis(cfg, data_epoch_MEG_sel);
% 
% cfg            = [];
% cfg.method     = 'coh';
% cfg.channel    = {'all',maxsens{1}};
% fdfourier      = ft_connectivityanalysis(cfg, freqfourier);
% 
% %     fdfourier_cmb = ft_combineplanar([],fdfourier);
% 
% cfg                  = [];
% cfg.parameter        = 'cohspctrm';
% % cfg.xlim             = [5 80];
% % cfg.refchannel       = maxsens{:};
% cfg.refchannel       = maxsens{1};
% cfg.layout           = 'neuromag306all.lay';
% cfg.showlabels       = 'yes';
% figure; ft_multiplotER(cfg, fdfourier)
% 
% cfg.xlim = [16 16];
% figure; ft_topoplotER(cfg, fdfourier)
% 
% 
% freq_cmb = ft_combineplanar([],freqfourier);
% 
% cfg = [];
% cfg.frequency = [16 16];
% temp = ft_selectdata(
% 
% [~,indx]        = sort(FFT.powspctrm,'descend');
% for isens = 1 : 9
%     maxsens{isens} = FFT.label{indx(isens)};
% end
% 
% 
% 
% 
% 
% 
% freq_cmb = ft_combineplanar([],freq);
% 
% 
% 
% cfg = [];
% cfg.channel = maxsens;
% cfg.avgoverchan = 'yes';
% temp = ft_selectdata(cfg,data_epoch_MEG);
% temp.label{1} = 'max';
% 
% test = ft_appenddata([],data_epoch_MEG,temp);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % put data in matrix
% dat = NaN(size(data_epoch_MEG.trial,2),size(data_epoch_MEG.label,1),maxlength);
% for itrial = 1 : size(data_epoch_MEG.trial,2)
%     fprintf('Adding trial %d to data matrix \n',itrial);
%     triallength = size(data_epoch_MEG.trial{itrial},2);
%     dat(itrial,:,maxlength-triallength+1:maxlength) = data_epoch_MEG.trial{itrial};
% end
% 
% % make index according to median split based on ratings
% F = ceil(2 * tiedrank(data_epoch_MEG.trialinfo(:,2)) / length(data_epoch_MEG.trialinfo(:,2)));
% 
% % before equalizing high/low bins
% rating_split        = ones(size(F));
% rating_split(F==2)  = 2;
% low_cnt             = size(find(rating_split == 1),1);
% high_cnt            = size(find(rating_split == 2),1);
% 
% % manually average over all trials, excluding nans
% cfg = [];
% cfg.vartrllength    = 2;
% ERF_probe           = ft_timelockanalysis(cfg,data_epoch_MEG);
% ERF_probe.avg       = squeeze(nanmean(dat,1));
% ERF_probe           = rmfield(ERF_probe,'var');
% 
% for irand = 1 : 24
%     rating_split       = ones(size(F));
%     rating_split(F==2) = 2;
%     if high_cnt > low_cnt
%         diff_bar                = min(data_epoch_MEG.trialinfo(rating_split == 2,2)); % rating bin between high/low split
%         nr_to_remove            = high_cnt-low_cnt;
%         discarted               = randsample(find(data_epoch_MEG.trialinfo(:,2) == diff_bar & rating_split == 2),nr_to_remove);
%         rating_split(discarted) = 0;
%     elseif high_cnt < low_cnt
%         diff_bar                = max(data_epoch_MEG.trialinfo(rating_split == 1,2)); % rating bin between high/low split
%         nr_to_remove            = low_cnt-high_cnt;
%         discarted               = randsample(find(data_epoch_MEG.trialinfo(:,2) == diff_bar & rating_split == 1),nr_to_remove);
%         rating_split(discarted) = 0;
%     else
%         discarted = [];
%     end
%     fprintf('I made %d low ratings, and %d high ratings in randomization %d \n',sum(rating_split==1),sum(rating_split==2),irand);
%     
%     temp_low(irand,:,:)  = squeeze(nanmean(dat(rating_split==1,:,:),1));
%     temp_high(irand,:,:) = squeeze(nanmean(dat(rating_split==2,:,:),1));
% end
% 
% % average over randomizations
% ERF_probe_high        = ERF_probe;
% ERF_probe_high.avg    = squeeze(nanmean(temp_high,1));
% ERF_probe_low         = ERF_probe;
% ERF_probe_low.avg     = squeeze(nanmean(temp_low,1));
% clear temp_low temp_high
% 
% % save nr of ratings per high/low bin
% for irating = 1 : 7
%     rating_cnt_low(irating)       = size(find(data_epoch_MEG.trialinfo(:,2) == irating & rating_split == 1),1);
%     rating_cnt_high(irating)      = size(find(data_epoch_MEG.trialinfo(:,2) == irating & rating_split == 2),1);
%     rating_cnt_discarted(irating) = size(find(data_epoch_MEG.trialinfo(discarted,2) == irating),1);
% end
% 
% % nr of observations per timepoint
% clear dof*
% for itime = 1 : maxlength
%     dof(itime)      = sum(~isnan(dat(:,1,itime)));
%     dof_high(itime) = sum(~isnan(dat(rating_split==2,1,itime)));
%     dof_low(itime)  = sum(~isnan(dat(rating_split==1,1,itime)));
% end
% 
% ERF_probe.dof             = dof;
% ERF_probe_high.dof        = dof_high;
% ERF_probe_low.dof         = dof_low;
% ERF_probe_high.trialinfo  = rating_cnt_high;
% ERF_probe_low.trialinfo   = rating_cnt_low;
% 
% lastnan = find(isnan(ERF_probe_low.avg(1,:)),1,'last');
% if ~isempty(lastnan)
%     ERF_probe_low.avg   = ERF_probe_low.avg(:,lastnan+1:end);
%     ERF_probe_low.dof   = ERF_probe_low.dof(lastnan+1:end);
%     ERF_probe_low.time  = ERF_probe_low.time(lastnan+1:end);
% end
% lastnan = find(isnan(ERF_probe_high.avg(1,:)),1,'last');
% if ~isempty(lastnan)
%     ERF_probe_high.avg  = ERF_probe_high.avg(:,lastnan+1:end);
%     ERF_probe_high.dof  = ERF_probe_high.dof(lastnan+1:end);
%     ERF_probe_high.time = ERF_probe_high.time(lastnan+1:end);
% end
% 
% % plot nr of observations per timepoint
% %     figure; hold; plot(ERF_probe.dof); plot(ERF_probe_high.dof); plot(ERF_probe_low.dof); axis tight
% 
% clear avg data_epoch_MEG
% 
% %     cfg = [];
% %     cfg.latency = [-2 2];
% %     cfg.channel = 'all'; % if not, combined grad is lost
% %     ERF = ft_selectdata(cfg,ERF);
% %
% %     cfg = [];
% %     cfg.lpfilter = 'yes';
% %     cfg.lpfreq = 30;
% %     ERF = ft_preprocessing(cfg,ERF);
% %
% % cfg = [];
% % cfg.channel                 = 'MEG';
% % cfg.latency                 = [-1.5 10];
% % ERF_probe                   = ft_selectdata(cfg,ERF_probe);
% % ERF_probe_high              = ft_selectdata(cfg,ERF_probe_high);
% % ERF_probe_low               = ft_selectdata(cfg,ERF_probe_low);
% % ERF_cue                     = ft_selectdata(cfg,ERF_probe);
% % ERF_cue_high                = ft_selectdata(cfg,ERF_cue_high);
% % ERF_cue_low                 = ft_selectdata(cfg,ERF_cue_low);
% 
% % Calculate TFR
% cfg                         = [];
% cfg.pad                     = 'nextpow2';
% cfg.channel                 = 'all';
% cfg.method                  = 'mtmconvol';
% cfg.foi                     = 1:1:30;
% cfg.taper                   = 'hanning';
% cfg.t_ftimwin               = ones(size(cfg.foi));
% 
% cfg.toi                     = -1:0.05:30;
% TFR_cue                     = ft_freqanalysis(cfg, ERF_cue);
% TFR_cue_high                = ft_freqanalysis(cfg, ERF_cue_high);
% TFR_cue_low                 = ft_freqanalysis(cfg, ERF_cue_low);
% 
% cfg.toi                     = -30:0.05:1;
% TFR_probe                   = ft_freqanalysis(cfg, ERF_probe);
% TFR_probe_high              = ft_freqanalysis(cfg, ERF_probe_high);
% TFR_probe_low               = ft_freqanalysis(cfg, ERF_probe_low);
% 
% % Combine planar
% TFR_probe                   = ft_combineplanar([],TFR_probe);
% TFR_probe_high              = ft_combineplanar([],TFR_probe_high);
% TFR_probe_low               = ft_combineplanar([],TFR_probe_low);
% TFR_cue                     = ft_combineplanar([],TFR_cue);
% TFR_cue_high                = ft_combineplanar([],TFR_cue_high);
% TFR_cue_low                 = ft_combineplanar([],TFR_cue_low);
% 
% % Diff
% TFR_probe_diff              = TFR_probe_high;
% TFR_probe_diff.powspctrm    = TFR_probe_low.powspctrm   - TFR_probe_high.powspctrm;
% TFR_cue_diff                = TFR_cue_high;
% TFR_cue_diff.powspctrm      = TFR_cue_low.powspctrm     - TFR_cue_high.powspctrm;
% 
% fprintf('Saving data, this will take a while... \n');
% save(fname_ERF,'ERF*','TFR*','-v7.3');
% 
% end
% % %% plotting
% %
% % % TFR average over SSEF channels
% % cfg         = [];
% % cfg.layout  = 'neuromag306cmb';
% % cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143'};
% % % cfg.channel = 'MEG*3';
% % cfg.ylim    = [5 30];
% % cfg.xlim = [-10 1];
% % % % cfg.baselinetype = 'relative';
% % cfg.zlim = 'absmax';
% % % cfg.baseline = [-10 1];
% %
% % fig = figure; ft_singleplotTFR(cfg,TFR_probe_diff);
% %
% % cfg.xlim = [-1 10];
% % fig = figure; ft_singleplotTFR(cfg,TFR_cue_diff);
% %
% % title('');
% % saveas(fig,['d:\analysis\WANDER\images\ERF\s' num2str(isubject) '_SSEF_' timing '.jpg']);
% %
