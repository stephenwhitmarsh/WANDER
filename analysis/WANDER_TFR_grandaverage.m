function WANDER_TFR_grandaverage(timing)
i = 1;

% load single subject TFR averages
for isubject = 1:26
    fname_TFR_avg = ['h:\analysis\WANDER\data\TFR_avg\s' num2str(isubject) '_TFR_avg_' timing '.mat'];
    fprintf('Trying to load subject %d\n',isubject);
    try
        temp = load(fname_TFR_avg);
        TFR_avg{i} = temp.TFR_avg;
        i = i + 1;
        fprintf('OK\n');
        clear temp
    catch
        fprintf('FAIL\n');
    end
end

% baseline corrections
for isubject = 1:size(TFR_avg,2)
    cfg = [];
    cfg.baseline = [-1 -0.5];
    cfg.baselinetype = 'relative';
    TFR_avg_bl{isubject} = ft_freqbaseline(cfg,TFR_avg{isubject}) ;
end

% grand average
GA_bl = ft_freqgrandaverage(cfg,TFR_avg_bl{:});

% TFR multiplot average
fig         = figure;
cfg         = [];
cfg.layout  = 'neuromag306cmb';
cfg.channel = 'MEGGRAD';
cfg.zlim    = [0.3 1.7];
%     cfg.baseline = [-1 0];
cfg.xlim = [-1 30];
cfg.baselinetype = 'relative';
ft_singleplotTFR(cfg,GA_bl);

