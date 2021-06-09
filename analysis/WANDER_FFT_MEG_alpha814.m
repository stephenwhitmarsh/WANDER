function [FFT_MEG] = WANDER_FFT_MEG_alpha(isubject,force,timing,rootpath,restingstate)

if restingstate == 1
    if rootpath == 1
        fname_FFT_MEG = ['w:\WANDER\data\FFT\s' num2str(isubject) '_FFT_MEG_rs.mat'];
    else
        fname_FFT_MEG = ['/shared/projects/project_wander/WANDER/data/FFT/s' num2str(isubject) '_FFT_MEG_rs.mat'];
    end
else
    if rootpath == 1
        fname_FFT_MEG = ['w:\WANDER\data\FFT\s' num2str(isubject) '_FFT_MEG.mat'];
    else
        fname_FFT_MEG = ['/shared/projects/project_wander/WANDER/data/FFT/s' num2str(isubject) '_FFT_MEG.mat'];
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
        
        % load artefact definition
        artdef_MEG = WANDER_artefact_detection_MEG(isubject,0,rootpath,restingstate);
        
        % artifact rejection - replace with NaN
        cfg = [];
        cfg.artfctdef           = artdef_MEG;
        cfg.artfctdef.reject    = 'nan';
        data_epoch_MEG  = ft_rejectartifact(cfg,data_epoch_MEG);
        
        % Welch method for FFT
        cfg = [];
        cfg.length                  = 8;
        cfg.overlap                 = 0.75;
        data_cut                    = ft_redefinetrial(cfg,data_epoch_MEG);
        %         clear data_epoch_MEG
        
        % select magnetometers
        cfg                     = [];
        cfg.channel             = 'MEG*1';
        data_cut_mag            = ft_selectdata(cfg,data_cut);
        
        % FFT
        cfg                     = [];
        cfg.pad                 = 'nextpow2';
        cfg.method              = 'mtmfft';
        cfg.keeptrials          = 'yes';
        cfg.channel             = 'all';
        cfg.foilim              = [1, 30];
        cfg.taper               = 'hanning';
        FFT_MEG_mag             = ft_freqanalysis(cfg, data_cut_mag);
        FFT_MEG_mag.dimord      = 'chan_freq';
        FFT_MEG_mag.powspctrm   = squeeze(nanmean(FFT_MEG_mag.powspctrm,1));
        
        % select magnetometers
        cfg                     = [];
        cfg.channel             = {'MEG*2','MEG*3'};
        data_cut_grad           = ft_selectdata(cfg,data_cut);
        
        % FFT
        cfg                     = [];
        cfg.pad                 = 'nextpow2';
        cfg.method              = 'mtmfft';
        cfg.keeptrials          = 'yes';
        cfg.channel             = 'all';
        cfg.foilim              = [1, 30];
        cfg.taper               = 'hanning';
        FFT_MEG_grad            = ft_freqanalysis(cfg, data_cut_grad);
        FFT_MEG_grad.dimord     = 'chan_freq';
        FFT_MEG_grad.powspctrm  = squeeze(nanmean(FFT_MEG_grad.powspctrm,1));
        
    else
        if strcmp(timing,'cue')
            data_epoch_MEG            = WANDER_ICA(isubject,0,rootpath,0);
            
            %             data_epoch_MEG = data_epoch_MEG.data_ICA;
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
        cfg.taper       = 'hanning';
%         cfg.tapsmofrq   = 1;
        cfg.trials      = data_epoch_MEG_comb.trialinfo(:,3) == 4; % only correct rejections
        FFT_MEG         = ft_freqanalysis(cfg, data_epoch_MEG_comb);
        clear data_epoch_MEG_comb
        
    end
    
    %
    % if isubject == 13
    %     low_freq = 9;
    %     high_freq = 11;
    % elseif isubject == 23
    %      low_freq = 10;
    %     high_freq = 11;
    % elseif isubject == 7
    %      low_freq = 13;
    %     high_freq = 14;
    % else

    low_freq = 7;
    high_freq = 14;
    if isubject == 16 || isubject == 22
        low_freq = 8;
    end    
    % average over FREQ
    cfg = [];
    cfg.avgoverfreq = 'yes';
    cfg.avgoverchan = 'no';
    cfg.frequency   = [low_freq high_freq];
    % cfg.channel     = {'MEG1631', 'MEG1741', 'MEG1843', 'MEG1912', 'MEG1913', 'MEG1911', 'MEG1923', 'MEG1922', 'MEG1921', 'MEG1932', 'MEG1933', 'MEG1931', 'MEG1941', 'MEG2013', 'MEG2012', 'MEG2011', 'MEG2023', 'MEG2022', 'MEG2021', 'MEG2032', 'MEG2033', 'MEG2031', 'MEG2042', 'MEG2043', 'MEG2041', 'MEG2113', 'MEG2112', 'MEG2111', 'MEG2122', 'MEG2123', 'MEG2121', 'MEG2133', 'MEG2132', 'MEG2131', 'MEG2143', 'MEG2142', 'MEG2141', 'MEG2233', 'MEG2312', 'MEG2313', 'MEG2311', 'MEG2323', 'MEG2322', 'MEG2321', 'MEG2332', 'MEG2333', 'MEG2331', 'MEG2343', 'MEG2342', 'MEG2341', 'MEG2442', 'MEG2441', 'MEG2512', 'MEG2513', 'MEG2543', 'MEG2542', 'MEG2541'};
    cfg.channel = {'MEG0131', 'MEG0141', 'MEG0211', 'MEG0221', 'MEG0231', 'MEG0241', 'MEG0411', 'MEG0421', 'MEG0431', 'MEG0441', 'MEG0631', 'MEG0641', 'MEG0711', 'MEG0721', 'MEG0731', 'MEG0741', 'MEG1041', 'MEG1111', 'MEG1121', 'MEG1131', 'MEG1141', 'MEG1311', 'MEG1321', 'MEG1331', 'MEG1341', 'MEG1431', 'MEG1441', 'MEG1511', 'MEG1541', 'MEG1611', 'MEG1621', 'MEG1811', 'MEG1821', 'MEG1831', 'MEG2211', 'MEG2221', 'MEG2241', 'MEG2411', 'MEG2421', 'MEG2611', 'MEG2621'};
    
    FFT_MEG_freqavg = ft_selectdata(cfg,FFT_MEG);
    
    % get peak channels
    [sorted_pow, sorted_idx] = sort(FFT_MEG_freqavg.powspctrm,'descend');
    max_chani  = sorted_idx(1:5);
    max_chan = [];
    for i = 1 : 5
        max_chan = [max_chan; FFT_MEG_freqavg.label{max_chani(i)}];
    end
    max_chan = cellstr(max_chan);
    
    % average over max channels
    cfg = [];
    cfg.avgoverchan = 'yes';
    cfg.channel = max_chan;
    FFT_MEG_chanavg = ft_selectdata(cfg,FFT_MEG);
    
    
    % get peak frequency edges
    low_freq_indx   = find(FFT_MEG.freq >= low_freq, 1, 'first');
    high_freq_indx  = find(FFT_MEG.freq <= high_freq, 1, 'last');
    
    % find maximum power within bounds for max electrode
%     [~,max_freq_indx] = max(FFT_MEG_chanavg.powspctrm(low_freq_indx:high_freq_indx));
    
    % replace the above with this, on 3-9-2018
    [~,max_freq_indx] = findpeaks(FFT_MEG_chanavg.powspctrm(low_freq_indx:high_freq_indx),'SortStr','descend','NPeaks',1);    
    
    % recover frequency index
    max_freq_indx           = max_freq_indx + low_freq_indx - 1;
    
    % get corresponding frequency
    max_freq                = FFT_MEG.freq(max_freq_indx);
    
    % return to datastructure
    FFT_MEG.max_freq     = max_freq;
    FFT_MEG.max_chan     = max_chan;
    
    % save data
    save(fname_FFT_MEG,'FFT_MEG');
    
    % TFR multiplot average
    if rootpath == 1
        fig         = figure;
        
        subplot(2,1,2);
        cfg         = [];
        cfg.layout  = 'neuromag306mag';
        cfg.parameter = 'powspctrm';
        cfg.channel = FFT_MEG.max_chan;
        ft_singleplotER(cfg,FFT_MEG);
        ax = axis;
        line([FFT_MEG.max_freq, FFT_MEG.max_freq],[ax(3) ax(4)],'color','k');
        
        subplot(2,1,1);
        cfg         = [];
        cfg.layout  = 'neuromag306mag';
        cfg.parameter = 'powspctrm';
        cfg.xlim = [FFT_MEG.max_freq FFT_MEG.max_freq];
        cfg.marker = 'off';
        cfg.highlight = 'on';
        cfg.highlightchannel = FFT_MEG.max_chan;
        cfg.highlightsymbol = 'pentagram';
        cfg.highlightsize = 14;
        cfg.comment = 'no';
        ft_topoplotER(cfg,FFT_MEG);
        
        saveas(fig,['/shared/projects/project_wander/WANDER\images\FFT\s' num2str(isubject) '_max_alpha_FFT.jpg']);
    end
end