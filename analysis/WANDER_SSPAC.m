% t = 1:2000;
% ti = 2000-t;
% p = mod(ti,160);
%
% figure; plot(p)

function [MI,MI_high,MI_low,PAC,PAC_high,PAC_low] = WANDER_SSPAC(isubject,force,timing,latency,rootpath)

MEGwindow = 0.3; % extending artefact rejection

if ~(strcmp(timing,'cue') || strcmp(timing,'probe'))
    fprintf('Use "cue" or "probe" as third argument.\n');
    return
end

if strcmp(latency,'all')
    if rootpath == 1
        fname_SSPAC = ['i:\analysis\WANDER\data\SSPAC\s' num2str(isubject) '_SSPAC_' timing '.mat'];
    else
        fname_SSPAC = ['/home/swhitmarsh/WANDER/data/SSPAC/s' num2str(isubject) '_SSPAC_' timing '.mat'];
    end
else
    if rootpath == 1
        fname_SSPAC = ['i:\analysis\WANDER\data\SSPAC\s' num2str(isubject) '_SSPAC_' num2str(latency(1)) '_' num2str(latency(2)) '.mat'];
    else
        fname_SSPAC = ['/home/swhitmarsh/WANDER/data/SSPAC/s' num2str(isubject) '_SSPAC_' num2str(latency(1)) '_' num2str(latency(2)) '.mat'];
    end  
end

if exist(fname_SSPAC,'file') && force ~= 1
    fprintf('Steady-state PAC \n');
    load(fname_SSPAC);
else
    fprintf('Steady-state PAC not found, creating it now! \n');
    data_epoch_MEG = WANDER_ICA(isubject,0,rootpath,0);
    
    % only correct rejections
    for ipart = 1 : 4
        cfg = [];
        cfg.trials = data_epoch_MEG{ipart}.trialinfo(:,3) == 4;
        data_epoch_MEG{ipart} = ft_selectdata(cfg,data_epoch_MEG{ipart});
    end
    
    % re-cut and re-align data to probe
    if strcmp(timing,'probe')
        for ipart = 1 : 4
            for itrial = 1 : size(data_epoch_MEG{ipart}.trial,2)
                fprintf('Cutting off end of trial %d block %d \n',itrial,ipart);
                data_epoch_MEG{ipart}.trial{itrial}         = data_epoch_MEG{ipart}.trial{itrial}(:,1:end-1000);
                data_epoch_MEG{ipart}.time{itrial}          = data_epoch_MEG{ipart}.time{itrial}(1:end-1000);
                data_epoch_MEG{ipart}.sampleinfo(itrial,2)  = data_epoch_MEG{ipart}.sampleinfo(itrial,2)-1000;
                
                fprintf('Cutting off beginning of trial %d block %d \n',itrial,ipart);
                data_epoch_MEG{ipart}.trial{itrial}         = data_epoch_MEG{ipart}.trial{itrial}(:,1501:end);
                data_epoch_MEG{ipart}.time{itrial}          = data_epoch_MEG{ipart}.time{itrial}(1501:end);
                data_epoch_MEG{ipart}.sampleinfo(itrial,1)  = data_epoch_MEG{ipart}.sampleinfo(itrial,1)+1500;
                
                fprintf('Shifting timecourse of trial %d block %d \n',itrial,ipart);
                data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial} - data_epoch_MEG{ipart}.time{itrial}(end);
            end
        end
    else
        for ipart = 1 : 4
            for itrial = 1 : size(data_epoch_MEG{ipart}.trial,2)
                fprintf('Cutting off end of trial %d block %d \n',itrial,ipart);
                data_epoch_MEG{ipart}.trial{itrial}         = data_epoch_MEG{ipart}.trial{itrial}(:,1:end-1000);
                data_epoch_MEG{ipart}.time{itrial}          = data_epoch_MEG{ipart}.time{itrial}(1:end-1000);
                data_epoch_MEG{ipart}.sampleinfo(itrial,2)  = data_epoch_MEG{ipart}.sampleinfo(itrial,2)-1000;
                
                fprintf('Cutting off beginning of trial %d block %d \n',itrial,ipart);
                data_epoch_MEG{ipart}.trial{itrial}         = data_epoch_MEG{ipart}.trial{itrial}(:,1501:end);
                data_epoch_MEG{ipart}.time{itrial}          = data_epoch_MEG{ipart}.time{itrial}(1501:end);
                data_epoch_MEG{ipart}.sampleinfo(itrial,1)  = data_epoch_MEG{ipart}.sampleinfo(itrial,1)+1500;
            end
        end
    end
    
    % combine blocks of artifact definitions
    temp          = WANDER_blink_detection(isubject,0,rootpath,0);
    artdef_EOG    = [temp{1}.zvalue.artifact; ...
        temp{2}.zvalue.artifact + data_epoch_MEG{1}.sampleinfo(2); ...
        temp{3}.zvalue.artifact + data_epoch_MEG{1}.sampleinfo(2) + data_epoch_MEG{2}.sampleinfo(2); ...
        temp{4}.zvalue.artifact + data_epoch_MEG{1}.sampleinfo(2) + data_epoch_MEG{2}.sampleinfo(2) + data_epoch_MEG{3}.sampleinfo(2)];
    
    temp          = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
    artdef_MEG    = [temp{1}.MEG.artifact; ...
        temp{2}.MEG.artifact    + data_epoch_MEG{1}.sampleinfo(2); ...
        temp{3}.MEG.artifact    + data_epoch_MEG{1}.sampleinfo(2) + data_epoch_MEG{2}.sampleinfo(2); ...
        temp{4}.MEG.artifact    + data_epoch_MEG{1}.sampleinfo(2) + data_epoch_MEG{2}.sampleinfo(2) + data_epoch_MEG{3}.sampleinfo(2)];
    
    artdef_muscle = [temp{1}.muscle.artifact; ...
        temp{2}.muscle.artifact + data_epoch_MEG{1}.sampleinfo(2); ...
        temp{3}.muscle.artifact + data_epoch_MEG{1}.sampleinfo(2) + data_epoch_MEG{2}.sampleinfo(2); ...
        temp{4}.muscle.artifact + data_epoch_MEG{1}.sampleinfo(2) + data_epoch_MEG{2}.sampleinfo(2) + data_epoch_MEG{3}.sampleinfo(2)];
    
    % extend artifact intervals according to predefined filtering window
    for iartefact = 1:size(artdef_MEG,1)
        artdef_ext_MEG(iartefact,1)    = artdef_MEG(iartefact,1)    - data_epoch_MEG{ipart}.fsample*MEGwindow/2;
        artdef_ext_MEG(iartefact,2)    = artdef_MEG(iartefact,2)    + data_epoch_MEG{ipart}.fsample*MEGwindow/2;
    end
    for iartefact = 1:size(artdef_muscle,1)
        artdef_ext_muscle(iartefact,1) = artdef_muscle(iartefact,1) - data_epoch_MEG{ipart}.fsample*MEGwindow/2;
        artdef_ext_muscle(iartefact,2) = artdef_muscle(iartefact,2) + data_epoch_MEG{ipart}.fsample*MEGwindow/2;
    end
    for iartefact = 1:size(artdef_EOG,1)
        artdef_ext_EOG(iartefact,1)    = artdef_EOG(iartefact,1)    - data_epoch_MEG{ipart}.fsample*MEGwindow/2;
        artdef_ext_EOG(iartefact,2)    = artdef_EOG(iartefact,2)    + data_epoch_MEG{ipart}.fsample*MEGwindow/2;
    end
    
    % create arrays of indexes for artefact samples - No EGG
    MEG_samples          = [];
    EOG_samples          = [];
    muscle_samples       = [];
    MEG_samples_ext      = [];
    EOG_samples_ext      = [];
    muscle_samples_ext   = [];
    try
        for iart = 1:size(artdef_ext_MEG,1)
            MEG_samples      = [MEG_samples     artdef_MEG(iart,1)     : artdef_MEG(iart,2)];
            MEG_samples_ext  = [MEG_samples_ext artdef_ext_MEG(iart,1) : artdef_ext_MEG(iart,2)];
        end
        fprintf('%d MEG artifacts found \n',size(artdef_ext_MEG,1));
    catch
        fprintf('No MEG artifacts found \n');
    end
    try
        for iart = 1:size(artdef_ext_EOG,1)
            EOG_samples      = [EOG_samples     artdef_EOG(iart,1)     : artdef_EOG(iart,2)];
            EOG_samples_ext  = [EOG_samples_ext artdef_ext_EOG(iart,1) : artdef_ext_EOG(iart,2)];
        end
        fprintf('%d blinks found \n',size(artdef_ext_EOG,1));
    catch
        fprintf('No blinks found \n');
    end
    try
        for iart = 1:size(artdef_ext_muscle,1)
            muscle_samples       = [muscle_samples     artdef_muscle(iart,1)     : artdef_muscle(iart,2)];
            muscle_samples_ext   = [muscle_samples_ext artdef_ext_muscle(iart,1) : artdef_ext_muscle(iart,2)];
        end
        fprintf('%d muscle artifacts found \n',size(artdef_ext_muscle,1));
    catch
        fprintf('No muscle artifacts found \n');
    end
    
    % concatinate data while maintain samplingo
    sampleinfo    = [data_epoch_MEG{1}.sampleinfo; ...
        data_epoch_MEG{2}.sampleinfo + data_epoch_MEG{1}.sampleinfo(2); ...
        data_epoch_MEG{3}.sampleinfo + data_epoch_MEG{1}.sampleinfo(2) + data_epoch_MEG{2}.sampleinfo(2); ...
        data_epoch_MEG{4}.sampleinfo + data_epoch_MEG{1}.sampleinfo(2) + data_epoch_MEG{2}.sampleinfo(2) + data_epoch_MEG{3}.sampleinfo(2)];
    
    data_epoch_MEG = ft_appenddata([],data_epoch_MEG{:});
    
    clear overlap*
    % map artefact samples onto boolean array of samples in trials
    for itrial = 1:size(data_epoch_MEG.time,2)
        fprintf('Calculating overlap with MEG artifacts, trial%d\n',itrial);
        overlap_MEG{itrial}          = ismember(sampleinfo(itrial,1):sampleinfo(itrial,2),MEG_samples);
        overlap_MEG_ext{itrial}      = ismember(sampleinfo(itrial,1):sampleinfo(itrial,2),MEG_samples_ext);
    end
    for itrial = 1:size(data_epoch_MEG.time,2)
        fprintf('Calculating overlap with EOG artifacts, trial%d\n',itrial);
        overlap_EOG{itrial}          = ismember(sampleinfo(itrial,1):sampleinfo(itrial,2),EOG_samples);
        overlap_EOG_ext{itrial}      = ismember(sampleinfo(itrial,1):sampleinfo(itrial,2),EOG_samples_ext);
    end
    for itrial = 1:size(data_epoch_MEG.time,2)
        fprintf('Calculating overlap with muscle artifacts, trial%d\n',itrial);
        overlap_muscle{itrial}       = ismember(sampleinfo(itrial,1):sampleinfo(itrial,2),muscle_samples);
        overlap_muscle_ext{itrial}   = ismember(sampleinfo(itrial,1):sampleinfo(itrial,2),muscle_samples_ext);
    end
    
    % mask whole trials that have a combined rejection of samples above threshold ratio
    for itrial = 1:size(data_epoch_MEG.time,2)
        if mean(overlap_muscle_ext{itrial} | overlap_MEG_ext{itrial}) > 0.5
            overlap_muscle_ext{itrial}  = true(size(overlap_muscle_ext{itrial}));
            overlap_MEG{itrial}         = true(size(overlap_MEG{itrial}));
            fprintf('Masking whole trial %d\n',itrial);
        end
    end
    
    matrix_EGG       = zeros(size(overlap_MEG,2),30000);
    matrix_MEG       = zeros(size(overlap_MEG,2),30000);
    matrix_EOG       = zeros(size(overlap_MEG,2),30000);
    matrix_jump      = zeros(size(overlap_MEG,2),30000);
    matrix_muscle    = zeros(size(overlap_MEG,2),30000);
    
    for itrial = 1:size(data_epoch_MEG.time,2)
        matrix_MEG(itrial,    ~overlap_MEG{itrial})        = 1;
        matrix_MEG(itrial,     overlap_MEG_ext{itrial})    = 3;
        matrix_MEG(itrial,     overlap_MEG{itrial})        = 2;
        matrix_EOG(itrial,    ~overlap_EOG{itrial})        = 1;
        matrix_EOG(itrial,     overlap_EOG_ext{itrial})    = 3;
        matrix_EOG(itrial,     overlap_EOG{itrial})        = 2;
        matrix_muscle(itrial, ~overlap_muscle{itrial})     = 1;
        matrix_muscle(itrial,  overlap_muscle_ext{itrial}) = 3;
        matrix_muscle(itrial,  overlap_muscle{itrial})     = 2;
    end
    
    figure;
    subplot(3,1,1); imagesc(matrix_MEG); colorbar;
    subplot(3,1,2); imagesc(matrix_EOG);colorbar;
    subplot(3,1,3); imagesc(matrix_muscle);colorbar;
% 
%     % maximum length trial
%     maxlength = 0;
%     minlength = 30000;
%     
%     for itrial = 1 : size(data_epoch_MEG.trial,2)
%         if size(data_epoch_MEG.trial{itrial},2) > maxlength
%             maxlength = size(data_epoch_MEG.trial{itrial},2);
%         end
%         if size(data_epoch_MEG.trial{itrial},2) < minlength
%             minlength = size(data_epoch_MEG.trial{itrial},2);
%         end
%     end
    
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
        if strcmp(latency,'all')
            if strcmp(timing,'probe')
                cfg.toi         = -30:0.013:0;
            else
                cfg.toi         = 0:0.013:30;
            end
        else
            cfg.toi         = latency(1):0.013:latency(2);
        end
        cfg.t_ftimwin   = (1./cfg.foi)*4;
        cfg.trials      = itrial;
        trialTFR        = ft_freqanalysis(cfg, data_epoch_MEG);
        
        % combine planar
        trialTFR        = ft_combineplanar([],trialTFR);
        
        % remove power estimates during artefact periods
        clear artefact_mask
        for itime = 1 : length(trialTFR.time)
                
                % find the closest sample in the timecourse corresponding to the current TFR bin
                [timediff,indx] = min(abs(data_epoch_MEG.time{itrial} - trialTFR.time(itime)));
                  
                % as a sanity check, save the distance in time between timecourse sample and time of TFR window
                dist(itrial,itime)      = timediff;
                
                % as a sanity check, save the time index
                time_indx(itrial,itime) = indx;
                
                % create artifact mask
                if size(overlap_MEG_ext{itrial},2) == size(data_epoch_MEG.time{itrial},2)
                    artefact_mask(itime) = overlap_MEG_ext{itrial}(indx) | overlap_muscle_ext{itrial}(indx);
                else
                    fprintf('LOST SYNCHRONIZATION BETWEEN ARTEFACTS AND TIME SERIES!!!\n');
                end
        end
        
        % remove artefacts
        artefact_mask = reshape(artefact_mask, 1, 1, size(artefact_mask,2));
        artefact_mask = repmat(artefact_mask, size(trialTFR.powspctrm,1), size(trialTFR.powspctrm,2), 1);
        trialTFR.powspctrm(artefact_mask) = NaN;
        clear artefact_mask
        
        % prepare phase bin
        phase_nrbins        = 18; %32;                                              % The number of phase bins = best to have a prime number
        phase_bins          = -pi:2*pi/phase_nrbins:pi;                        % The extreme phases of each bin
        x = trialTFR.time ./ (1/16);
        phase = (x - floor(x))*2*pi-pi;
        [~,~,phase_indx]    = histcounts(phase,phase_bins);                                                         % bin phase, check usage e.g. with: [count, edges, index] = histcounts([1 2 3 1 2 3 6 7 8 4 2 1],[1 2 3]);
        phase_indx(phase_indx == 0) = NaN;
        for ifreq = 1 : size(trialTFR.powspctrm,2)                                                                  % loop through all frequencies
            fprintf('working on freq %d, trial %d of %d \n',ifreq,itrial,size(data_epoch_MEG.trial,2));
            for ichan_MEG = 1 : size(trialTFR.powspctrm,1)                                                          % loop through all sensors
                for ibin = unique(phase_indx(~isnan(phase_indx)))
                    power_bincount(itrial,ichan_MEG,ifreq,ibin)  = sum(~isnan(trialTFR.powspctrm(ichan_MEG,ifreq,phase_indx == ibin)),3);
                    power_binmean(itrial,ichan_MEG,ifreq,ibin)   = nanmean(trialTFR.powspctrm(ichan_MEG,ifreq,phase_indx == ibin),3);
                    power_binmedian(itrial,ichan_MEG,ifreq,ibin) = nanmedian(trialTFR.powspctrm(ichan_MEG,ifreq,phase_indx == ibin),3);
                end
            end
        end
%         figure; hist(phase_indx);
%         figure; plot(phase_indx);
%         figure; plot(squeeze(power_bincount(2,1,1,:)))
% plot(squeeze(power_bincount(179,1,1,:)))
    end % itrial
    
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
    
    for irating = 1 : 8
        MI_low.trialinfo(irating)  = size(find(data_epoch_MEG.trialinfo(:,2) == irating & rating_split == 1),1);
        MI_high.trialinfo(irating) = size(find(data_epoch_MEG.trialinfo(:,2) == irating & rating_split == 2),1);
    end
    
    %normalize median split with average over all trials
%     PAC_binmean_norm      = nan(size(PAC_binmean));
%     PAC_binmean_high_norm = nan(size(PAC_binmean_high));
%     PAC_binmean_low_norm  = nan(size(PAC_binmean_low));
%     for ichan = 1:size(PAC_binmean,1)
%         for ifreq = 1 : size(PAC_binmean,2)
%             PAC_binmean_norm(ichan,ifreq,:)      = PAC_binmean(ichan,ifreq,:);     % ./ nanmean(PAC_binmean(ichan,ifreq,:));
%             PAC_binmean_high_norm(ichan,ifreq,:) = PAC_binmean_high(ichan,ifreq,:);% ./ nanmean(PAC_binmean(ichan,ifreq,:));
%             PAC_binmean_low_norm(ichan,ifreq,:)  = PAC_binmean_low(ichan,ifreq,:); % ./ nanmean(PAC_binmean(ichan,ifreq,:));
%         end
%     end
    
%     MI_diff = MI;
%     MI_diff.avg = (MI_low.avg - MI_high.avg) ./ (MI_low.avg + MI_high.avg) ;
    
    % put it in a FieldTrip frequency data format
    PAC                = [];
    PAC.label           = trialTFR.label;
    PAC.dimord          = 'chan_freq_time';
    PAC.freq            = 30:100;
    PAC.time            = [1 : phase_nrbins];
    PAC.powspctrm       = PAC_binmean;
    
    PAC_high            = PAC;
    PAC_low             = PAC;
    PAC_high.powspctrm  = PAC_binmean_high;
    PAC_low.powspctrm   = PAC_binmean_low;
    
%     cfg                 = [];
%     cfg.layout          = 'neuromag306cmb.lay';
%     cfg.channel         = 'MEG*3';
% %     cfg.xlim            = [1 : phase_nrbins];
% %     cfg.ylim            = [30 100];
%     %     cfg.zlim            = [0.99 1.01];
%     figure; ft_singleplotTFR(cfg, PAC_norm)
%     figure; ft_singleplotTFR(cfg, PAC_low_norm)
%     figure; ft_singleplotTFR(cfg, PAC_diff)
%     
%     
%     PAC_norm = PAC_low;
%     PAC_norm.powspctrm = PAC.powspctrm ./ mean(PAC.powspctrm,3);
%     PAC_low_norm = PAC_low;
%     PAC_low_norm.powspctrm = PAC_low.powspctrm ./ ((mean(PAC_low.powspctrm,3) + mean(PAC_high.powspctrm,3))./2) ;
%     PAC_high_norm = PAC_high;
%     PAC_high_norm.powspctrm = PAC_high.powspctrm ./ ((mean(PAC_low.powspctrm,3) + mean(PAC_high.powspctrm,3))./2);
%     
%     PAC_diff = PAC_low;
%     PAC_diff.powspctrm = PAC_low_norm.powspctrm - PAC_high_norm.powspctrm;
% 
% %     figure; ft_singleplotTFR(cfg, PAC_diff)
%     
%     cfg                 = [];
%     cfg.layout          = 'neuromag306cmb.lay';
%     cfg.channel         = 'MEG*3';
%     cfg.xlim            = [30 100];
%     cfg.ylim            = 'maxabs';
%     figure; ft_singleplotER(cfg, MI)
%     
%     MI_diff = MI_high;
%     MI_diff.avg = (MI_low.avg - MI_high.avg) ./ (MI_low.avg + MI_high.avg);
%         
%     figure; ft_singleplotER(cfg, MI_diff)
%     figure; ft_singleplotER(cfg, MI_low, MI_high)
%     
%     cfg.layout          = 'neuromag306mag.lay';
%     cfg.channel         = 'MEG*1';
%     cfg.xlim            = [30 100];
%     cfg.ylim            = 'maxabs';
%     figure; ft_singleplotER(cfg, MI)
%     figure; ft_singleplotER(cfg, MI_diff)
%     figure; ft_singleplotER(cfg, MI_low, MI_high)
    
    save(fname_SSPAC,'MI','MI_high','MI_low','PAC','PAC_high','PAC_low');
end