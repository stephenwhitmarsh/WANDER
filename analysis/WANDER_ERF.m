function ERF = WANDER_ERF(isubject,force,timing,rootpath)

if ~(strcmp(timing,'cue') || strcmp(timing,'probe'))
    fprintf('Use "cue" or "probe" as third argument.\n');
    return
end

if rootpath == 1
    fname_ERF = ['i:\analysis\WANDER\data\ERF\s' num2str(isubject) '_' timing '.mat'];
else
    fname_ERF = ['/shared/projects/project_wander/WANDER/data/ERF/s' num2str(isubject) '_' timing '.mat'];
end

if exist(fname_ERF,'file') && force ~= 1
    fprintf('Returning ERF\n');
    load(fname_ERF);
else
    fprintf('ERF not found, creating it now! \n');
    data_epoch_MEG = WANDER_ICA(isubject,0,rootpath,0);
    artdef = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
    
    % start with cuelocked
    for ipart = 1 : 4
        cfg = [];
        cfg.reject                      = 'nan'; % make bugreport, does not work with cfg.artfctdef.reject
        cfg.artfctdef.minaccepttim      = 0;
        cfg.artfctdef                   = artdef{ipart};
        data_epoch_MEG{ipart}           = ft_rejectartifact(cfg,data_epoch_MEG{ipart});
        for itrial = 1 : size(data_epoch_MEG{ipart}.trial,2)
            fprintf('Cutting off end of trial %d block %d \n',itrial,ipart);
            data_epoch_MEG{ipart}.trial{itrial} = data_epoch_MEG{ipart}.trial{itrial}(:,1:end-1000);
            data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial}(1:end-1000);
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
    for itrial = 1 : size(data_epoch_MEG.trial,2)
        if size(data_epoch_MEG.trial{itrial},2) > maxlength
            maxlength = size(data_epoch_MEG.trial{itrial},2);
        end
    end
    
    % put data in matrix
    dat = NaN(size(data_epoch_MEG.trial,2),size(data_epoch_MEG.label,1),maxlength);
    for itrial = 1 : size(data_epoch_MEG.trial,2)
        fprintf('Adding trial %d to data matrix \n',itrial);
        dat(itrial,:,1:size(data_epoch_MEG.trial{itrial},2)) = data_epoch_MEG.trial{itrial};
    end
    
    % make index according to median split based on ratings
    F = ceil(2 * tiedrank(data_epoch_MEG.trialinfo(:,2)) / length(data_epoch_MEG.trialinfo(:,2)));
    
    % before equalizing high/low bins
    rating_split        = ones(size(F));
    rating_split(F==2)  = 2;
    low_cnt             = size(find(rating_split == 1),1);
    high_cnt            = size(find(rating_split == 2),1);
    
    % manually average over all trials, excluding nans
    cfg = [];
    cfg.vartrllength    = 2;
    ERF_cue             = ft_timelockanalysis(cfg,data_epoch_MEG);
    ERF_cue.avg         = squeeze(nanmean(dat,1));
    ERF_cue             = rmfield(ERF_cue,'var');
    
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
        
        temp_low(irand,:,:)  = squeeze(nanmean(dat(rating_split==1,:,:),1));
        temp_high(irand,:,:) = squeeze(nanmean(dat(rating_split==2,:,:),1));
    end
    
    % average over randomizations
    ERF_cue_high        = ERF_cue;
    ERF_cue_high.avg    = squeeze(nanmean(temp_high,1));
    ERF_cue_low         = ERF_cue;
    ERF_cue_low.avg     = squeeze(nanmean(temp_low,1));
    clear temp_low temp_high
    
    for irating = 1 : 7
        rating_cnt_low(irating)       = size(find(data_epoch_MEG.trialinfo(:,2) == irating & rating_split == 1),1);
        rating_cnt_high(irating)      = size(find(data_epoch_MEG.trialinfo(:,2) == irating & rating_split == 2),1);
        rating_cnt_discarted(irating) = size(find(data_epoch_MEG.trialinfo(discarted,2) == irating),1);
    end
    
    % nr of observations per timepoint
    for itime = 1 : maxlength
        dof(itime)      = sum(~isnan(dat(:,1,itime)));
        dof_high(itime) = sum(~isnan(dat(rating_split==2,1,itime)));
        dof_low(itime)  = sum(~isnan(dat(rating_split==1,1,itime)));
    end
    
    ERF_cue.dof             = dof;
    ERF_cue_high.dof        = dof_high;
    ERF_cue_low.dof         = dof_low;
    ERF_cue_high.trialinfo  = rating_cnt_high;
    ERF_cue_low.trialinfo   = rating_cnt_low;
    
    firstnan = find(isnan(ERF_cue_low.avg(1,:)),1,'first');
    if ~isempty(firstnan)
        ERF_cue_low.avg     = ERF_cue_low.avg(:,1:firstnan-1);
        ERF_cue_low.dof     = ERF_cue_low.dof(1:firstnan-1);
        ERF_cue_low.time    = ERF_cue_low.time(1:firstnan-1);
    end
    firstnan = find(isnan(ERF_cue_high.avg(1,:)),1,'first');
    if ~isempty(firstnan)
        ERF_cue_high.avg    = ERF_cue_high.avg(:,1:firstnan-1);
        ERF_cue_high.dof    = ERF_cue_high.dof(1:firstnan-1);
        ERF_cue_high.time   = ERF_cue_high.time(1:firstnan-1);
    end
    
    % plot nr of observations per timepoint
%     figure; hold; plot(ERF_cue.dof); plot(ERF_cue_high.dof); plot(ERF_cue_low.dof); axis tight
    
    clear avg;
    
    % now lock to end of trial
    data_epoch_MEG          = WANDER_ICA(isubject,0,rootpath,0);
    artdef                  = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
    
    % repeat for probelocked
    for ipart = 1 : 4
        cfg = [];
        cfg.reject                      = 'nan'; % make bugreport, does not work with cfg.artfctdef.reject
        cfg.artfctdef.minaccepttim      = 0;
        cfg.artfctdef                   = artdef{ipart};
        data_epoch_MEG{ipart}           = ft_rejectartifact(cfg,data_epoch_MEG{ipart});
        for itrial = 1 : size(data_epoch_MEG{ipart}.trial,2)
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
    for itrial = 1 : size(data_epoch_MEG.trial,2)
        if size(data_epoch_MEG.trial{itrial},2) > maxlength
            maxlength = size(data_epoch_MEG.trial{itrial},2);
        end
    end
    
    % put data in matrix
    dat = NaN(size(data_epoch_MEG.trial,2),size(data_epoch_MEG.label,1),maxlength);
    for itrial = 1 : size(data_epoch_MEG.trial,2)
        fprintf('Adding trial %d to data matrix \n',itrial);
        triallength = size(data_epoch_MEG.trial{itrial},2);
        dat(itrial,:,maxlength-triallength+1:maxlength) = data_epoch_MEG.trial{itrial};
    end
    
    % make index according to median split based on ratings
    F = ceil(2 * tiedrank(data_epoch_MEG.trialinfo(:,2)) / length(data_epoch_MEG.trialinfo(:,2)));
    
    % before equalizing high/low bins
    rating_split        = ones(size(F));
    rating_split(F==2)  = 2;
    low_cnt             = size(find(rating_split == 1),1);
    high_cnt            = size(find(rating_split == 2),1);
    
    % manually average over all trials, excluding nans
    cfg = [];
    cfg.vartrllength    = 2;
    ERF_probe           = ft_timelockanalysis(cfg,data_epoch_MEG);
    ERF_probe.avg       = squeeze(nanmean(dat,1));
    ERF_probe           = rmfield(ERF_probe,'var');
    
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
        
        temp_low(irand,:,:)  = squeeze(nanmean(dat(rating_split==1,:,:),1));
        temp_high(irand,:,:) = squeeze(nanmean(dat(rating_split==2,:,:),1));
    end
    
    % average over randomizations
    ERF_probe_high        = ERF_probe;
    ERF_probe_high.avg    = squeeze(nanmean(temp_high,1));
    ERF_probe_low         = ERF_probe;
    ERF_probe_low.avg     = squeeze(nanmean(temp_low,1));
    clear temp_low temp_high
    
    % save nr of ratings per high/low bin
    for irating = 1 : 7
        rating_cnt_low(irating)       = size(find(data_epoch_MEG.trialinfo(:,2) == irating & rating_split == 1),1);
        rating_cnt_high(irating)      = size(find(data_epoch_MEG.trialinfo(:,2) == irating & rating_split == 2),1);
        rating_cnt_discarted(irating) = size(find(data_epoch_MEG.trialinfo(discarted,2) == irating),1);
    end
    
    % nr of observations per timepoint
    clear dof*
    for itime = 1 : maxlength
        dof(itime)      = sum(~isnan(dat(:,1,itime)));
        dof_high(itime) = sum(~isnan(dat(rating_split==2,1,itime)));
        dof_low(itime)  = sum(~isnan(dat(rating_split==1,1,itime)));
    end
    
    ERF_probe.dof             = dof;
    ERF_probe_high.dof        = dof_high;
    ERF_probe_low.dof         = dof_low;
    ERF_probe_high.trialinfo  = rating_cnt_high;
    ERF_probe_low.trialinfo   = rating_cnt_low;
    
    lastnan = find(isnan(ERF_probe_low.avg(1,:)),1,'last');
    if ~isempty(lastnan)
        ERF_probe_low.avg   = ERF_probe_low.avg(:,lastnan+1:end);
        ERF_probe_low.dof   = ERF_probe_low.dof(lastnan+1:end);
        ERF_probe_low.time  = ERF_probe_low.time(lastnan+1:end);
    end
    lastnan = find(isnan(ERF_probe_high.avg(1,:)),1,'last');
    if ~isempty(lastnan)
        ERF_probe_high.avg  = ERF_probe_high.avg(:,lastnan+1:end);
        ERF_probe_high.dof  = ERF_probe_high.dof(lastnan+1:end);
        ERF_probe_high.time = ERF_probe_high.time(lastnan+1:end);
    end
    
    % plot nr of observations per timepoint
%     figure; hold; plot(ERF_probe.dof); plot(ERF_probe_high.dof); plot(ERF_probe_low.dof); axis tight
    
    clear avg data_epoch_MEG
    
    %     cfg = [];
    %     cfg.latency = [-2 2];
    %     cfg.channel = 'all'; % if not, combined grad is lost
    %     ERF = ft_selectdata(cfg,ERF);
    %
    %     cfg = [];
    %     cfg.lpfilter = 'yes';
    %     cfg.lpfreq = 30;
    %     ERF = ft_preprocessing(cfg,ERF);
    %
    % cfg = [];
    % cfg.channel                 = 'MEG';
    % cfg.latency                 = [-1.5 10];
    % ERF_probe                   = ft_selectdata(cfg,ERF_probe);
    % ERF_probe_high              = ft_selectdata(cfg,ERF_probe_high);
    % ERF_probe_low               = ft_selectdata(cfg,ERF_probe_low);
    % ERF_cue                     = ft_selectdata(cfg,ERF_probe);
    % ERF_cue_high                = ft_selectdata(cfg,ERF_cue_high);
    % ERF_cue_low                 = ft_selectdata(cfg,ERF_cue_low);
    
    % Calculate TFR
    cfg                         = [];
    cfg.pad                     = 'nextpow2';
    cfg.channel                 = 'all';
    cfg.method                  = 'mtmconvol';
    cfg.foi                     = 1:1:30;
    cfg.taper                   = 'hanning';
    cfg.t_ftimwin               = ones(size(cfg.foi));
    
    cfg.toi                     = -1:0.05:30;
    TFR_cue                     = ft_freqanalysis(cfg, ERF_cue);
    TFR_cue_high                = ft_freqanalysis(cfg, ERF_cue_high);
    TFR_cue_low                 = ft_freqanalysis(cfg, ERF_cue_low);
    
    cfg.toi                     = -30:0.05:1;
    TFR_probe                   = ft_freqanalysis(cfg, ERF_probe);
    TFR_probe_high              = ft_freqanalysis(cfg, ERF_probe_high);
    TFR_probe_low               = ft_freqanalysis(cfg, ERF_probe_low);
    
    % Combine planar
    TFR_probe                   = ft_combineplanar([],TFR_probe);
    TFR_probe_high              = ft_combineplanar([],TFR_probe_high);
    TFR_probe_low               = ft_combineplanar([],TFR_probe_low);
    TFR_cue                     = ft_combineplanar([],TFR_cue);
    TFR_cue_high                = ft_combineplanar([],TFR_cue_high);
    TFR_cue_low                 = ft_combineplanar([],TFR_cue_low);
    
    % Diff
    TFR_probe_diff              = TFR_probe_high;
    TFR_probe_diff.powspctrm    = TFR_probe_low.powspctrm   - TFR_probe_high.powspctrm;
    TFR_cue_diff                = TFR_cue_high;
    TFR_cue_diff.powspctrm      = TFR_cue_low.powspctrm     - TFR_cue_high.powspctrm;
    
    fprintf('Saving data, this will take a while... \n');
    save(fname_ERF,'ERF*','TFR*','-v7.3');
    
end
    % %% plotting
    %
    % % TFR average over SSEF channels
    % cfg         = [];
    % cfg.layout  = 'neuromag306cmb';
    % cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143'};
    % % cfg.channel = 'MEG*3';
    % cfg.ylim    = [5 30];
    % cfg.xlim = [-10 1];
    % % % cfg.baselinetype = 'relative';
    % cfg.zlim = 'absmax';
    % % cfg.baseline = [-10 1];
    %
    % fig = figure; ft_singleplotTFR(cfg,TFR_probe_diff);
    %
    % cfg.xlim = [-1 10];
    % fig = figure; ft_singleplotTFR(cfg,TFR_cue_diff);
    %
    % title('');
    % saveas(fig,['d:\analysis\WANDER\images\ERF\s' num2str(isubject) '_SSEF_' timing '.jpg']);
    %
