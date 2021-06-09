function [TFR_median_split] = WANDER_TFR_median_split_GA(force,timing,rootpath)

if ~(strcmp(timing,'cue') || strcmp(timing,'probe'))
    fprintf('Use "cue" or "probe" as third argument.\n');
    return
end

% timing = 'probe';

if rootpath == 1
    fname_TFR_median_split = ['h:/analysis/WANDER/data/TFR/TFR_median_split_' timing '.mat'];
else
    fname_TFR_median_split = ['/home/swhitmarsh/WANDER/data/TFR/TFR_median_split_' timing '.mat'];
end

if exist(fname_TFR_median_split,'file') && force ~= 1
    fprintf('Returning TFR median split \n');
    load(fname_TFR_median_split);
else
    fprintf('TFR correlation not found, creating it now! \n');
    
%     for isubject = [1:20 22:26]
    for isubject = [1:26]
        
        fprintf('Subject %d \n',isubject);

        TFR = WANDER_TFR(isubject,0,timing,rootpath,0);
        
        % only correct rejections
        % only selected timesegment, to be implemented (just add to argin)
        for ipart = 1 : 4
            TFR{ipart}  = ft_combineplanar([],TFR{ipart});
            cfg         = [];
            cfg.trials  = TFR{ipart}.trialinfo(:,3) == 4;
%             cfg.latency = latency;
            TFR{ipart} = ft_selectdata(cfg,TFR{ipart});
        end
        
        % combine blocks
        TFR_comb{isubject}             = TFR{1};
        TFR_comb{isubject}.powspctrm   = [TFR{1}.powspctrm; TFR{2}.powspctrm; TFR{3}.powspctrm; TFR{4}.powspctrm];
        TFR_comb{isubject}.trialinfo   = [TFR{1}.trialinfo; TFR{2}.trialinfo; TFR{3}.trialinfo; TFR{4}.trialinfo];        
        TFR_comb{isubject}             = rmfield(TFR_comb{isubject},'cumtapcnt');
        TFR_comb{isubject}             = rmfield(TFR_comb{isubject},'grad');
        clear TFR
        
        % load artifact definition and reshape to match TFR
        fname_TFR_phaselocking = ['h:\analysis\WANDER\data\TFR\s' num2str(isubject) '_TFR_phaselocked_' timing '.mat'];
        load(fname_TFR_phaselocking,'artefact_mask');
                
        % remove measures during artefact periods
        TFR_comb{isubject}.powspctrm(repmat(artefact_mask,1,size(TFR_comb{isubject}.powspctrm,2),size(TFR_comb{isubject}.powspctrm,3),1)) = NaN;
        
        % make index according to median split based on ratings
        F                            = ceil(2 * tiedrank(TFR_comb{isubject}.trialinfo(:,2)) / length(TFR_comb{isubject}.trialinfo(:,2)));
        rating_split{isubject}       = ones(size(F));
        rating_split{isubject}(F==2) = 2;
        
        % before rejection
        low_cnt  = size(find(rating_split{isubject} == 1),1);
        high_cnt = size(find(rating_split{isubject} == 2),1);
        
        if high_cnt > low_cnt
            diff_bar                = min(TFR_comb{isubject}.trialinfo(rating_split{isubject} == 2,2));
            nr_to_remove            = high_cnt-low_cnt;
            discarted               = randsample(find(TFR_comb{isubject}.trialinfo(:,2) == diff_bar & rating_split{isubject} == 2),nr_to_remove);
            rating_split{isubject}(discarted) = 0;
        elseif high_cnt < low_cnt
            diff_bar                = max(TFR_comb{isubject}.trialinfo(rating_split{isubject} == 1,2));
            nr_to_remove            = low_cnt-high_cnt;
            discarted               = randsample(find(TFR_comb{isubject}.trialinfo(:,2) == diff_bar & rating_split{isubject} == 1),nr_to_remove);
            rating_split{isubject}(discarted) = 0;
        else
            discarted = [];
        end
        
        for irating = 1 : 7
            rating_cnt_low(irating)       = size(find(TFR_comb{isubject}.trialinfo(:,2) == irating & rating_split{isubject} == 1),1);           
            rating_cnt_high(irating)      = size(find(TFR_comb{isubject}.trialinfo(:,2) == irating & rating_split{isubject} == 2),1);
            rating_cnt_discarted(irating) = size(find(TFR_comb{isubject}.trialinfo(discarted,2) == irating),1);           
        end
        
        % after rejection
        disc_cnt = size(find(rating_split{isubject} == 0),1);
        low_cnt  = size(find(rating_split{isubject} == 1),1);               
        high_cnt = size(find(rating_split{isubject} == 2),1);
        
        
        fig = figure;
        bar([rating_cnt_low;rating_cnt_high;rating_cnt_discarted]','stacked');
        title(['Subjectnr ' num2str(isubject)]);
        legend({['low: ' num2str(low_cnt)],['high: ' num2str(high_cnt)],['discarted: ' num2str(disc_cnt)]});
        print(fig,'-dpng',['d:\analysis\WANDER\images\TFR\s' num2str(isubject) '_' timing '_TFR_median_split.png']);
        
        TFR_comb{isubject}.trialinfo(:,7)  = rating_split{isubject};
        
        TFR_comb_low{isubject}             = TFR_comb{isubject};
        TFR_comb_low{isubject}.powspctrm   = TFR_comb{isubject}.powspctrm(rating_split{isubject}==1,:,:,:);
        TFR_comb_low{isubject}.trialinfo   = TFR_comb{isubject}.trialinfo(rating_split{isubject}==1,:);
        
        TFR_comb_high{isubject}            = TFR_comb{isubject};
        TFR_comb_high{isubject}.powspctrm  = TFR_comb{isubject}.powspctrm(rating_split{isubject}==2,:,:,:);
        TFR_comb_high{isubject}.trialinfo  = TFR_comb{isubject}.trialinfo(rating_split{isubject}==2,:);
        
        % average over trials
        cfg = [];
        cfg.avgoverrpt = 'yes';
        TFR_comb_high{isubject}.powspctrm  = squeeze(nanmean(TFR_comb_high{isubject}.powspctrm,1));
        TFR_comb_high{isubject}.dimord     = 'chan_freq_time';
        TFR_comb_low{isubject}.powspctrm   = squeeze(nanmean(TFR_comb_low{isubject}.powspctrm,1));
        TFR_comb_low{isubject}.dimord      = 'chan_freq_time';              
        TFR_comb_diff{isubject}            = TFR_comb_high{isubject};
        TFR_comb_diff{isubject}.powspctrm  = TFR_comb_high{isubject}.powspctrm ./ TFR_comb_low{isubject}.powspctrm;
        clear TFR_comb
    end
    
    cfg = [];
    TFR_comb_high_GA = ft_freqgrandaverage(cfg,TFR_comb_high{[1:20 22:26]});
    TFR_comb_low_GA  = ft_freqgrandaverage(cfg,TFR_comb_low{[1:20 22:26]});
    TFR_comb_diff_GA = ft_freqgrandaverage(cfg,TFR_comb_diff{[1:20 22:26]});
    
    save(fname_TFR_median_split,'TFR_comb_high_GA','TFR_comb_low_GA','TFR_comb_diff_GA','-v7.3');
end


% statistics
cfg = [];
cfg.channel     = 'all';
cfg.latency     = [-5 -1];
cfg.avgovertime = 'no';
cfg.frequency   = [8 14];
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'no';
cfg.parameter   = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.method           = 'analytic';

cfg.statistic           = 'ft_statfun_depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, TFR{1});

Nsub = 25;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_freqstatistics(cfg,TFR_comb_high{[1:20 22:end]},TFR_comb_low{[1:20 22:end]});


cfg             = [];
cfg.layout      = 'neuromag306cmb';
figure; ft_clusterplot(cfg,stat)

MI_diff = MI_high_GA{1};
MI_diff.avg = MI_high_GA{1}.avg - MI_low_GA{1}.avg;






% TFR multiplot average
fig             = figure;
cfg             = [];
cfg.layout      = 'neuromag306cmb';
% cfg.layout      = 'neuromag306mag';
cfg.zlim        = [0.85 1.15];
cfg.xlim        = [-10 0];
figure; ft_singleplotTFR(cfg,TFR_comb_diff_GA);

cfg.parameter = 'stat';
figure; ft_singleplotTFR(cfg,stat);

% 9, 10, 18 is atypical, 11,12, 25 somewhat
figure; ft_singleplotTFR(cfg,TFR_comb_diff{2});


for i = [1:20 22:26]
subplot(5,5,i);
cfg             = [];
cfg.layout      = 'neuromag306cmb';
cfg.zlim        = 'maxabs';

ft_singleplotTFR(cfg,TFR_comb_diff_GA);
end

% ax = axis;
% perc = squeeze(TFR_corr.nr) ./ max(TFR_corr.nr) * ax(4);
% plot(TFR_corr.time,perc,'r','linewidth',2);
% print(fig,'-dpdf',['d:\analysis\WANDER\images\TFR\s' num2str(isubject) '_' timing '_TFR_corr_cue.pdf']);
% print(fig,'-dpng',['d:\analysis\WANDER\images\TFR\s' num2str(isubject) '_' timing '_TFR_corr_cue.png']);


