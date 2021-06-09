
i = 1;
for isubject = [1:20 22:26]
    TFR_corr{i} = WANDER_TFR_corr(isubject,0,'probe');
    TFR_corr{i} = rmfield(TFR_corr{i},'powspctrm');
    TFR_corr{i}.dimord = 'chan_freq_time';
    
    % average over trials
%     TFR_corr{i}.rho       = squeeze(nanmean(TFR_corr{i}.rho,1));
%     TFR_corr{i}.pval      = squeeze(nanmean(TFR_corr{i}.pval,1));
%     TFR_corr{i}.dimord     = 'chan_freq_time';
    i = i + 1;
%     
end

cfg = [];
cfg.parameter = 'rho';
TFR_corr_GA = ft_freqgrandaverage(cfg,TFR_corr{:})


    % TFR multiplot average
    fig         = figure;
    cfg         = [];
    cfg.layout  = 'neuromag306cmb';
%     cfg.channel = 'MEG*1';
    cfg.parameter = 'rho'
    cfg.xlim = [-5 0];
    cfg.ylim = [16 16];
    cfg.zlim = 'maxabs';
    ft_topoplotTFR(cfg,TFR_corr{1});
    
        TFR_corr{i}.dimord = 'chan_freq_time';

        ft_topoplotTFR(cfg,TFR_corr_avg);

    
    
    fig         = figure;
    cfg         = [];
    cfg.layout  = 'neuromag306cmb';
    cfg.channel = 'MEG*3';
    cfg.parameter = 'rho';
%     cfg.zlim    = [0.3 1.7];
%     cfg.baseline = [-1 0];
%     cfg.xlim = [-1 30];
%     cfg.baselinetype = 'relative';
    ft_singleplotTFR(cfg,TFR_corr_avg);