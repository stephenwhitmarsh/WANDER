function WANDER_source_individual_alpha_GRAD(isubject,force,rootpath)

disp(['Working on Subjectnr. ' num2str(isubject)]);

if rootpath == 1
    fname_TFR_phaselocked = ['w:\WANDER\data\beamformer\s' num2str(isubject) '_MI_source_grad_individual.mat'];
else
    fname_TFR_phaselocked = ['/shared/projects/project_wander/WANDER/data/beamformer/s' num2str(isubject) '_MI_source_grad_individual.mat'];
end

if exist(fname_TFR_phaselocked,'file') && force ~= 1
    fprintf('Returning TFR-EGG phaselocking \n');
    load(fname_TFR_phaselocked);
else
    fprintf('TFR-EGG phaselocking not found, creating it now! \n');
    
    % load data
    data_MEG            = WANDER_ICA(isubject,0,rootpath,0);

    for iblock = 1 : 4
        
        % select all channels, only correct rejections
        cfg                 = [];
        cfg.channel         = {'MEG*2','MEG*3'};
        cfg.trials          = find(data_MEG{iblock}.trialinfo(:,3) == 4);
        data_MEG{iblock}    = ft_selectdata(cfg,data_MEG{iblock});
        data_MEG{iblock}    = rmfield(data_MEG{iblock},'cfg');
      
        % remove artefacts, as well as trials with > 20% rejection
        artdef      = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
        
        % remove artefactual segments
        cfg                             = [];
        cfg.reject                      = 'partial'; 
        cfg.artfctdef                   = artdef{iblock};
        cfg.artfctdef.minaccepttim      = 1;
        data_MEG_clean{iblock}          = ft_rejectartifact(cfg,data_MEG{iblock});
 
%         % find trials that have < 20% data removed que to artefacts
%         trials_to_keep = [];
%         for itrial = unique(data_MEG_clean{iblock}.trialinfo(:,6))'
%             if any(data_orig.trialinfo(:,8) == itrial)
%                 nr_segments_orig = sum(data_orig.trialinfo(:,8) == itrial);
%                 nr_segments_art  = sum(data_MEG_timebinned{iblock}.trialinfo(:,8) == itrial);
%                 if nr_segments_orig ~= nr_segments_art
%                     fprintf('%0.2f percent removed from trialnr %d', 100-nr_segments_art/nr_segments_orig*100, itrial);
%                     if 100-nr_segments_art/nr_segments_orig*100 > 20
%                         fprintf(': removing whole trial\n');
%                     else
%                         fprintf('\n');
%                         trials_to_keep = [trials_to_keep itrial];
%                     end
%                 else
%                     trials_to_keep = [trials_to_keep itrial];
%                 end
%             end
%         end
%         clear data_orig
        
%         % only keep trials that have < 20% data removed que to artefacts
%         cfg = [];
%         cfg.trials = find(ismember(data_MEG_timebinned{iblock}.trialinfo(:,8),trials_to_keep));
%         data_MEG_timebinned{iblock} = ft_selectdata(cfg,data_MEG_timebinned{iblock});
%         data_EGG_timebinned{iblock} = ft_selectdata(cfg,data_EGG_timebinned{iblock});

%         %% EGG phase binning
%         
%         % add phase of EGG to trialinfo (phase at middle of epoch)
%         for itrial = 1 : size(data_MEG_timebinned{iblock}.trial,2)
%             data_MEG_timebinned{iblock}.trialinfo(itrial,10) = data_EGG_timebinned{iblock}.trial{itrial}(1000);
%         end
%         
%         % convert and enter EGG phase bins
%         phase_nrbins    = 18;                                              % The number of phase bins
%         phase_bins      = -pi:2*pi/phase_nrbins:pi;                        % The extreme phases of each bin
%         phase_axis      = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;   % The midpoint of each phase bin
%         
%         % bin phase, check usage with e.g.: [count, edges, index] = histcounts([1 2 3 1 2 3 6 7 8 4 2 1],[1 2 3]);
%         [~,~,data_MEG_timebinned{iblock}.trialinfo(:,11)] = histcounts(data_MEG_timebinned{iblock}.trialinfo(:,10),phase_bins);                                         
        
    end
    
    clear data_MEG
    data_MEG_append = ft_appenddata([],data_MEG_clean{:});
    clear data_MEG_clean 
%     
%     % give unique trialnumber
%     itrial = 1;
%     data_MEG_timebinned_append.trialinfo(1,12) = itrial;  
%     for i = 2 : size(data_MEG_timebinned_append.trialinfo,1)
%         if data_MEG_timebinned_append.trialinfo(i,6) ~= data_MEG_timebinned_append.trialinfo(i-1,6)
%             itrial = itrial + 1;
%         end
%         data_MEG_timebinned_append.trialinfo(i,12) = itrial;
%     end
    
    % load headmodel
    [mri_segmented, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(isubject,rootpath,0);
    
    % load sensor information
    [dataset_task, ~] = WANDER_subjectinfo;
    hdr = ft_read_header(dataset_task{isubject,1});
    
    % load common filter
    [common_filter, leadfield] = WANDER_common_filter_DICS_individual_onlyalpha(isubject,rootpath,0);

    cfg = [];
    % get individual alpha
    if isubject ~= 13
        [FFT_MEG_mag, ~] = WANDER_FFT_MEG_alpha(isubject,0,'cue',rootpath,1);
        cfg.foi = FFT_MEG_mag.max_freq;
    else
        cfg.foi = 11;
    end
    
    cfg.output          = 'powandcsd';
    cfg.method          = 'mtmfft';
    cfg.channel         = {'MEG*2', 'MEG*3'};
    cfg.keeptrials      = 'yes'; 
    cfg.tapsmofrq       = 0.5;
    
    cfg.pad             = 'nextpow2';
    FFT                 = ft_freqanalysis(cfg, data_MEG_append);
    
    
    
    
    
    
    
    
    
    
    source              = cell(1,size(FFT10.trialinfo,1));

    for itrial = 1 : size(FFT10.trialinfo,1)
        
        fprintf('Trialnr. %d\n',itrial);
        
        cfg              = [];
        cfg.trials       = itrial;
        FFT_sel10        = ft_selectdata(cfg, FFT10);
%         FFT_sel11        = ft_selectdata(cfg, FFT11);
      
        cfg              = [];
        cfg.channel      = {'MEG*2', 'MEG*3'};
        cfg.grad         = hdr.grad;
        cfg.headmodel    = headmodel;
        cfg.keeptrials   = 'yes';
        cfg.rawtrial     = 'no';
        cfg.method       = 'dics';
        cfg.grid         = leadfield;
        cfg.grid.filter  = common_filter;
        
        source10{itrial} = ft_sourceanalysis(cfg, FFT_sel10);
%         source11{itrial} = ft_sourceanalysis(cfg, FFT_sel11);
    end
    
    %% reorganize into pos_rpt dimord AND LOG-TRANSFORM

    source10_powrpt = rmfield(source10{1},{'avg','method','cfg'});
    source10_powrpt.dimord = 'pos_rpt';
    source10_powrpt.freq = source10{1}.freq;
%     source11_powrpt = rmfield(source11{1},{'avg','method','cfg'});
%     source11_powrpt.dimord = 'pos_rpt';
%     source11_powrpt.freq = source11{1}.freq;
    
    for itrial = 1 : size(source10,2)
        source10_powrpt.pow(:,itrial) = source10{itrial}.avg.pow;
        source10_powrpt.trialinfo(itrial,:) = FFT10.trialinfo(itrial,:);
%         source11_powrpt.pow(:,itrial) = log(source11{itrial}.avg.pow);
%         source11_powrpt.trialinfo(itrial,:) = FFT11.trialinfo(itrial,:);        
    end
    
    clear FFT* source10 source11 data_MEG_timebinned_append
    
    % create trial-based trialinfo
    for itrial = 1 : max(source10_powrpt.trialinfo(:,12))
        trials_trialinfo(itrial,:)  = source10_powrpt.trialinfo(find(source10_powrpt.trialinfo(:,12) == itrial,1,'first'),:);
        nrofobs(itrial)             = size(find(source10_powrpt.trialinfo(:,12) == itrial),1);
    end
    
    % median split based on observations (1-second segments)
    F = ceil(2 * tiedrank(source10_powrpt.trialinfo(:,2)) / length(source10_powrpt.trialinfo(:,2)));
    rating_split_obs        = ones(size(F));
    rating_split_obs(F==2)  = 2;
    
    % median split of trials based on median split from observations
    max_low                 = max(source10_powrpt.trialinfo(rating_split_obs == 1,2));
    rating_split_trial      = ones(1,size(trials_trialinfo,1));
    rating_split_trial(trials_trialinfo(:,2) > max_low) = 2;

    % average power over median split (without balancing) 
    source10_avg           = source10_powrpt;
    source10_avg.dimord    = 'pos';
    source10_avg.pow       = nanmean(source10_powrpt.pow,2);
    source10_avg_low       = source10_powrpt;
    source10_avg_low.pow   = nanmean(source10_powrpt.pow(:,rating_split_obs == 1),2);
    source10_avg_high      = source10_powrpt;
    source10_avg_high.pow  = nanmean(source10_powrpt.pow(:,rating_split_obs == 2),2);
%     
%     source11_avg           = source11_powrpt;
%     source11_avg.dimord    = 'pos';    
%     source11_avg.pow       = nanmean(source11_powrpt.pow,2);
%     source11_avg_low       = source11_powrpt;
%     source11_avg_low.pow   = nanmean(source11_powrpt.pow(:,rating_split_obs == 1),2);
%     source11_avg_high      = source11_powrpt;
%     source11_avg_high.pow  = nanmean(source11_powrpt.pow(:,rating_split_obs == 2),2);
    
    %% PAC per trial (without balancing) 
    PAC_trial_10 = nan(size(trials_trialinfo,1),phase_nrbins,size(source10_avg.pow,1));
%     PAC_trial_11 = nan(size(trials_trialinfo,1),phase_nrbins,size(source10_avg.pow,1));
    PAC_bins_10  = nan(size(trials_trialinfo,1),phase_nrbins);
%     PAC_bins_11  = nan(size(trials_trialinfo,1),phase_nrbins);
    
    for itrial = 1 : max(trials_trialinfo(:,12))
        fprintf('Trialnr. %d of %d\n',itrial,max(trials_trialinfo(:,12)));
        for ibin = 1 : phase_nrbins
            binindx = source10_avg.trialinfo(:,12) == itrial & source10_avg.trialinfo(:,11) == ibin;
            PAC_trial_10(itrial,ibin,:) = nanmean(source10_powrpt.pow(:,binindx),2);
            PAC_bins_10(itrial,ibin)    = size(find(binindx),1);
%             PAC_trial_11(itrial,ibin,:) = nanmean(source11_powrpt.pow(:,binindx),2);
%             PAC_bins_11(itrial,ibin)    = size(find(binindx),1);            
        end
    end

    PAC_orig10                  = source10_avg;
    PAC_orig10.pow              = permute(nanmean(PAC_trial_10,1),[3,2,1]);
    PAC_orig10.time             = 1:18; %bins
    PAC_orig10.dimord           = 'pos_time';
    PAC_orig10.nrobs            = PAC_bins_10;
    
    PAC_orig10_low              = source10_avg;
    PAC_orig10_low.pow          = permute(nanmean(PAC_trial_10(rating_split_trial == 1,:,:),1),[3,2,1]);
    PAC_orig10_low.time         = 1:18;
    PAC_orig10_low.dimord       = 'pos_time';
    PAC_orig10_low.nrobs        = PAC_bins_10(rating_split_trial == 1,:);
    
    PAC_orig10_high             = source10_avg;
    PAC_orig10_high.pow         = permute(nanmean(PAC_trial_10(rating_split_trial == 2,:,:),1),[3,2,1]);
    PAC_orig10_high.time        = 1:18;
    PAC_orig10_high.dimord      = 'pos_time';
    PAC_orig10_high.nrobs       = PAC_bins_10(rating_split_trial == 2,:);
    
%     PAC_orig11                  = source11_avg;
%     PAC_orig11.pow              = permute(nanmean(PAC_trial_11,1),[3,2,1]);
%     PAC_orig11.time             = 1:18; %bins
%     PAC_orig11.dimord           = 'pos_time';
%     PAC_orig11.nrobs            = PAC_bins_11;
%     
%     PAC_orig11_low              = source11_avg;
%     PAC_orig11_low.pow          = permute(nanmean(PAC_trial_11(rating_split_trial == 1,:,:),1),[3,2,1]);
%     PAC_orig11_low.time         = 1:18;
%     PAC_orig11_low.dimord       = 'pos_time';
%     PAC_orig11_low.nrobs        = PAC_bins_11(rating_split_trial == 1,:);
%     
%     PAC_orig11_high             = source11_avg;
%     PAC_orig11_high.pow         = permute(nanmean(PAC_trial_11(rating_split_trial == 2,:,:),1),[3,2,1]);
%     PAC_orig11_high.time        = 1:18;
%     PAC_orig11.dimord           = 'pos_time';
%     PAC_orig11_high.nrobs       = PAC_bins_11(rating_split_trial == 2,:);   
   
    clear PAC_trial* PAC_bins* rand*
    
    nrand = 24;
    
    % preallocate memory
    rand_MI10       = nan(size(PAC_orig10.pow,1),nrand);
    rand_MI10_high  = nan(size(PAC_orig10.pow,1),nrand);
    rand_MI10_low   = nan(size(PAC_orig10.pow,1),nrand);
%     rand_MI11       = nan(size(PAC_orig10.pow,1),nrand);
%     rand_MI11_high  = nan(size(PAC_orig10.pow,1),nrand);
%     rand_MI11_low   = nan(size(PAC_orig10.pow,1),nrand); 
    
    for irand = 1 : nrand
        
        fprintf('Calculating MI for randomization %d \n',irand);
        
        % refresh median split (based on nr. of obs.) because it will be edited each irand
        rating_split_trial = ones(1,size(trials_trialinfo,1));
        rating_split_trial(trials_trialinfo(:,2) > max_low) = 2;
        
        % compare nr. of observations
        nrofobs_sum_low    = nansum(nrofobs(rating_split_trial == 1));
        nrofobs_sum_high   = nansum(nrofobs(rating_split_trial == 2));
        nr_obs_to_remove = abs(sum(nrofobs_sum_high) - sum(nrofobs_sum_low));

        if nrofobs_sum_high > nrofobs_sum_low
            % if there are more high observations than low observations remove trials from lowest of high ratings
            diff_bar = min(source10_powrpt.trialinfo(rating_split_obs == 2,2)); % rating bin between high/low split
        else
            % if there are more low observations than high observations remove trials from highest of low ratings
            diff_bar = max(source10_powrpt.trialinfo(rating_split_obs == 1,2)); % rating bin between high/low split
        end
        
        % calculate if nr of observation in middle rating are enough to balance
        if nansum(nrofobs(trials_trialinfo(:,2) == diff_bar)) - nr_obs_to_remove >= 0
            fprintf('Cool... Enough bins in difference rating to balance number of observations \n');
            
            % randomize trials in bin that will be discarted
            trials_that_can_be_discarded = find(trials_trialinfo(:,2) == diff_bar);
            trials_that_can_be_discarded = trials_that_can_be_discarded(randperm(size(trials_that_can_be_discarded,1)));
            
            % cumulative distribution of total nr. of observations
            obs_cumsum = cumsum(nrofobs(trials_that_can_be_discarded));
            
            % find max nr. of trials that can be fully removed
            trials_to_discard = find(obs_cumsum <= nr_obs_to_remove);
            
            if ~isempty(trials_to_discard)
                % remove trials
                rating_split_trial(trials_that_can_be_discarded(trials_to_discard)) = 0;
            else
                fprintf('Do not need to remove a trial \n');
            end
            
            % how many observations left to remove?
            nr_obs_to_remove = nr_obs_to_remove - sum(nrofobs(trials_that_can_be_discarded(trials_to_discard)));
            obs_to_remain = ones(length(source10_powrpt.trialinfo),1);
          
            if nr_obs_to_remove > 0
                % determine trial from which observations will be removed
                if ~isempty(trials_to_discard)
                    % remove remaining observations from next trial in cumsum
                    crack_trial = trials_that_can_be_discarded(trials_to_discard(end) + 1);
                else
                    % remove remaining observations from next trial in cumsum
                    crack_trial = trials_that_can_be_discarded(1);
                end
                
                % randomly rotate observation indexes
                randobs = circshift(1:nrofobs(crack_trial),randi(nrofobs(crack_trial)));
                
                % remove observations from trial
                starttrial_indx = find(source10_powrpt.trialinfo(:,12) == crack_trial,1,'first');
                obs_to_remain(starttrial_indx + randobs(1:nr_obs_to_remove) - 1) = 0;
            else
                fprintf('******** No need to remove extra observations ********** \n');
            end
        else
            fprintf('******** Not enough bins in difference rating to balance number of observations ********** \n');
            break
        end
        
        % PAC per trial
        
        PAC_trial_10 = nan(size(source10_avg.pow,1),size(trials_trialinfo,1),phase_nrbins);
%         PAC_trial_11 = nan(size(source10_avg.pow,1),size(trials_trialinfo,1),phase_nrbins);        

        for itrial = 1 : max(source10_powrpt.trialinfo(:,12))
            for ibin = 1 : 18
                binindx = source10_powrpt.trialinfo(:,12) == itrial & source10_powrpt.trialinfo(:,11) == ibin & obs_to_remain;
                PAC_trial_10(:,itrial,ibin) = nanmean(source10_powrpt.pow(:,binindx),2);
%                 PAC_trial_11(:,itrial,ibin) = nanmean(source11_powrpt.pow(:,binindx),2);                
            end
        end
        
        PAC10         = squeeze(nanmean(PAC_trial_10,2));
        PAC10_high    = squeeze(nanmean(PAC_trial_10(:,rating_split_trial == 2,:),2));
        PAC10_low     = squeeze(nanmean(PAC_trial_10(:,rating_split_trial == 1,:),2));
%         PAC11         = squeeze(nanmean(PAC_trial_11,2));
%         PAC11_high    = squeeze(nanmean(PAC_trial_11(:,rating_split_trial == 2,:),2));
%         PAC11_low     = squeeze(nanmean(PAC_trial_11(:,rating_split_trial == 1,:),2));
%         
        clear PAC_trial_10 PAC_trial_11
        
        % Calculate MI per frequency and channel
        fprintf('Calculating MI, randomization %d \n',irand);   
        for ipos = 1 : size(PAC10,1)
            rand_MI10(ipos,irand)      = (log(phase_nrbins) - (-sum((PAC10(ipos,:)      / sum(PAC10(ipos,:)))      .* log((PAC10(ipos,:)      /sum(PAC10(ipos,:)))))))      /log(phase_nrbins);
            rand_MI10_high(ipos,irand) = (log(phase_nrbins) - (-sum((PAC10_high(ipos,:) / sum(PAC10_high(ipos,:))) .* log((PAC10_high(ipos,:) /sum(PAC10_high(ipos,:))))))) /log(phase_nrbins);
            rand_MI10_low(ipos,irand)  = (log(phase_nrbins) - (-sum((PAC10_low(ipos,:)  / sum(PAC10_low(ipos,:)))  .* log((PAC10_low(ipos,:)  /sum(PAC10_low(ipos,:)))))))  /log(phase_nrbins);
%             rand_MI11(ipos,irand)      = (log(phase_nrbins) - (-sum((PAC11(ipos,:)      / sum(PAC11(ipos,:)))      .* log((PAC11(ipos,:)      /sum(PAC11(ipos,:)))))))      /log(phase_nrbins);
%             rand_MI11_high(ipos,irand) = (log(phase_nrbins) - (-sum((PAC11_high(ipos,:) / sum(PAC11_high(ipos,:))) .* log((PAC11_high(ipos,:) /sum(PAC11_high(ipos,:))))))) /log(phase_nrbins);
%             rand_MI11_low(ipos,irand)  = (log(phase_nrbins) - (-sum((PAC11_low(ipos,:)  / sum(PAC11_low(ipos,:)))  .* log((PAC11_low(ipos,:)  /sum(PAC11_low(ipos,:)))))))  /log(phase_nrbins);
        end
        
        clear PAC10* PAC11*
    end % irand
    
    % prepare MI data structures
    MI10          = keepfields(source10_powrpt,{'dim','inside','pos'});
    MI10.dimord   = 'pos';
    MI10_high     = MI10;
    MI10_low      = MI10;
%     MI11          = MI10;    
%     MI11_high     = MI10;
%     MI11_low      = MI10;
%     
    % average MI over randomizations
    MI10.MI       = squeeze(nanmean(rand_MI10,2));
    MI10_high.MI  = squeeze(nanmean(rand_MI10_high,2));
    MI10_low.MI   = squeeze(nanmean(rand_MI10_low,2));
%     MI11.MI       = squeeze(nanmean(rand_MI11,2));
%     MI11_high.MI  = squeeze(nanmean(rand_MI11_high,2));
%     MI11_low.MI   = squeeze(nanmean(rand_MI11_low,2));
    
    % save data
    fprintf('Saving data in %s. This can take a while.\n',fname_TFR_phaselocked);
  
    save(fname_TFR_phaselocked, 'MI10','MI10_low','MI10_high', ...
                                'source10_avg','source10_avg_low','source10_avg_high', ... 
                                'PAC_orig10','PAC_orig10_low','PAC_orig10_high','-v7.3');
end