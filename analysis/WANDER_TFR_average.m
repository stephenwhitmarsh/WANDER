function WANDER_TFR_average(isubject, force, timing, rootpath)

timing = 'cue';
addpath d:/fieldtrip/
addpath('D:/analysis/WANDER/scripts/');
ft_defaults

if isempty(rootpath)
    rootpath = 'd:\analysis\WANDER\data\';
end

fname_TFR_avg = [rootpath filesep 'TFR_avg\s' num2str(isubject) '_TFR_avg_' timing '.mat'];

if exist(fname_TFR_avg,2) && force ~= 1
    load(fname_TFR_avg);
else
    TFR = WANDER_TFR(isubject,0,timing);
    
    for ipart = 1 : 4
        TFR_comb{ipart} = ft_combineplanar([],TFR{ipart});
    end
    
    % concatinate blocks
    TFR_conc = TFR_comb{1};
    TFR_conc.powspctrm = [TFR_comb{1}.powspctrm; TFR_comb{2}.powspctrm; TFR_comb{3}.powspctrm; TFR_comb{4}.powspctrm];
    TFR_conc.trialinfo = [TFR_comb{1}.trialinfo; TFR_comb{2}.trialinfo; TFR_comb{3}.trialinfo; TFR_comb{4}.trialinfo];
    TFR_conc.cumtapcnt = [TFR_comb{1}.cumtapcnt; TFR_comb{2}.cumtapcnt; TFR_comb{3}.cumtapcnt; TFR_comb{4}.cumtapcnt];
    
    % average over correct trials
    cfg                = [];
    cfg.avgoverrpt     = 'yes';
    cfg.trials         = find(TFR_conc.trialinfo(:,3) == 4);
    TFR_avg            = ft_selectdata(cfg,TFR_conc);
    
    % effecting nanmean over whole available time period
    TFR_avg.powspctrm = squeeze(nanmean(TFR_conc.powspctrm,1));
    
    % number of observations per timepoint
    for itime = 1:size(TFR_conc.time,2)
        TFR_avg.nr(itime) = sum(~isnan(TFR_conc.powspctrm(:,1,1,itime)));
    end
    
    % save data
    save(fname_TFR_avg,'TFR_avg','-v7.3');
end
