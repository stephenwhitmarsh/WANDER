function [FFT_MEG] = WANDER_FFT_MEG(isubject,force,timing,rootpath,restingstate)

if restingstate == 1
    if rootpath == 1
        fname_FFT_MEG = ['d:\analysis\WANDER\data\FFT\s' num2str(isubject) '_FFT_MEG_rs.mat'];
    else
        fname_FFT_MEG = ['/shared/projects/project_wander/WANDER/data/FFT\s' num2str(isubject) '_FFT_MEG_rs.mat'];
    end
else
    if rootpath == 1
        fname_FFT_MEG = ['d:\analysis\WANDER\data\FFT\s' num2str(isubject) '_FFT_MEG.mat'];
    else
        fname_FFT_MEG = ['/shared/projects/project_wander/WANDER/data/FFT\s' num2str(isubject) '_FFT_MEG.mat'];
    end
end

if exist(fname_FFT_MEG,'file') && force ~= 1
    fprintf('Returning FFT\n');
    load(fname_FFT_MEG);
else
    
    % process EGG with Welch method to determine max EGG freq and sensor
    fprintf('FFT not found. Making it now!\n');
    
    if restingstate == 1
        data_epoch_MEG = WANDER_ICA(isubject,0,rootpath,restingstate);

        cfg = [];
        cfg.length                  = 200;
        cfg.overlap                 = 0.75;
        data_cut                    = ft_redefinetrial(cfg,data_epoch_MEG);
        clear data_epoch_MEG
        
        % FFT
        cfg             = [];
        cfg.pad         = 'nextpow2';
        cfg.channel     = 'all';
        cfg.method      = 'mtmfft';
        cfg.keeptrials  = 'no';
        cfg.channel     = 'all';
        cfg.foilim      = [1 30];
        cfg.taper       = 'hanning';
        FFT_MEG         = ft_freqanalysis(cfg, data_cut);
        
    else
        if strcmp(timing,'cue')
            data_epoch_MEG = WANDER_ICA_round2(isubject,0);
            data_epoch_MEG = data_epoch_MEG.data_ICA;
        else
            data_epoch_MEG = WANDER_redefine_MEG_to_probe(isubject,0);
        end
        
        data_epoch_MEG_comb = ft_appenddata([],data_epoch_MEG{:});
        clear data_epoch_MEG
        
        % select time period
        cfg = [];
        if strcmp(timing,'cue')
            cfg.latency      = [0 30];
        end
        if strcmp(timing,'probe')
            cfg.latency      = [-30 0];
        end
        data_epoch_MEG_comb = ft_selectdata(cfg,data_epoch_MEG_comb);
        
        % FFT
        cfg             = [];
        cfg.pad         = 32; % nextpow2(30)
        cfg.channel     = 'all';
        cfg.method      = 'mtmfft';
        cfg.keeptrials  = 'no';
        cfg.channel     = {'MEG'};
        cfg.foilim      = [1 30];
        cfg.taper       = 'dpss';
        cfg.tapsmofrq   = 1;
        cfg.trials      = data_epoch_MEG_comb.trialinfo(:,3) == 4; % only correct rejections
        FFT_MEG         = ft_freqanalysis(cfg, data_epoch_MEG_comb);
        clear data_epoch_MEG_comb
        
    end
    
    
    % average over channels
    cfg = [];
    cfg.avgoverchan = 'yes';
    cfg.parameter   = {'powspctrm','powspctrm_log'};
    cfg.channel     = {'MEG0431', 'MEG0441', 'MEG0731', 'MEG0741', 'MEG1141', 'MEG1631', 'MEG1731', 'MEG1741', 'MEG1811', 'MEG1821', 'MEG1831', 'MEG1841', 'MEG1911', 'MEG1921', 'MEG1931', 'MEG1941', 'MEG2011', 'MEG2021', 'MEG2031', 'MEG2041', 'MEG2111', 'MEG2121', 'MEG2211', 'MEG2231', 'MEG2241', 'MEG2311', 'MEG2331', 'MEG2341'};
    % more posterior
    cfg.channel     = {'MEG1631', 'MEG1741', 'MEG1843', 'MEG1912', 'MEG1913', 'MEG1911', 'MEG1923', 'MEG1922', 'MEG1921', 'MEG1932', 'MEG1933', 'MEG1931', 'MEG1941', 'MEG2013', 'MEG2012', 'MEG2011', 'MEG2023', 'MEG2022', 'MEG2021', 'MEG2032', 'MEG2033', 'MEG2031', 'MEG2042', 'MEG2043', 'MEG2041', 'MEG2113', 'MEG2112', 'MEG2111', 'MEG2122', 'MEG2123', 'MEG2121', 'MEG2133', 'MEG2132', 'MEG2131', 'MEG2143', 'MEG2142', 'MEG2141', 'MEG2233', 'MEG2312', 'MEG2313', 'MEG2311', 'MEG2323', 'MEG2322', 'MEG2321', 'MEG2332', 'MEG2333', 'MEG2331', 'MEG2343', 'MEG2342', 'MEG2341', 'MEG2442', 'MEG2441', 'MEG2512', 'MEG2513', 'MEG2543', 'MEG2542', 'MEG2541'};
    FFT_MEG_avg     = ft_selectdata(cfg,FFT_MEG);
    FFT_MEG_avg.powspctrm_log = log(FFT_MEG_avg.powspctrm);
    
    % get peak frequency EGG bounds
    low_freq_indx   = find(FFT_MEG_avg.freq >= 7, 1, 'first');
    high_freq_indx  = find(FFT_MEG_avg.freq <= 14, 1, 'last');
    
    % find maximum power within bounds for each electrode
    [max_pow,max_freq_indx] = max(FFT_MEG_avg.powspctrm_log(:,low_freq_indx:high_freq_indx)');
    
    % recover frequency index
    max_freq_indx           = max_freq_indx + low_freq_indx - 1;
    
    % get corresponding frequency
    max_freq                = FFT_MEG.freq(max_freq_indx);
    
    % return to datastructure
    FFT_MEG.max_freq     = max_freq;
    
    % save data
    save(fname_FFT_MEG,'FFT_MEG');
    
    % plot FFT
    fig = figure;
    plot(FFT_MEG_avg.freq,FFT_MEG_avg.powspctrm_log);
    ax = axis;
    line([FFT_MEG.max_freq, FFT_MEG.max_freq],[ax(3) ax(4)]);
    saveas(fig,['/shared/projects/project_wander/WANDER\images\FFT\s' num2str(isubject) '_FFT_MEG.jpg']);
    
end
    % TFR multiplot average
%     fig         = figure;
%     cfg         = [];
%     cfg.layout  = 'neuromag306all';
%     cfg.channel = 'MEG*1';
%     FFT_MEG.powspctrm_log = log(FFT_MEG.powspctrm);
%     cfg.parameter = 'powspctrm_log';
%     cfg.freq = FFT_MEG.max_freq;
%     ft_singleplotER(cfg,FFT_MEG);

