function [FFT_EGG] = WANDER_FFT_EGG(isubject,force,rootpath,restingstate)

[dataset_task, dataset_rs] = WANDER_subjectinfo;

if restingstate == 1
    if rootpath == 1
        fname_FFT_EGG = ['W:\WANDER\data\FFT\s' num2str(isubject) '_FFT_EGG_rs.mat'];
    else
        fname_FFT_EGG = ['/shared/projects/project_wander/WANDER/data/FFT/s' num2str(isubject) '_FFT_EGG_rs.mat'];
    end
else
    if rootpath == 1
        fname_FFT_EGG = ['W:\WANDER\data\FFT\s' num2str(isubject) '_FFT_EGG.mat'];
    else
        fname_FFT_EGG = ['/shared/projects/project_wander/WANDER/data/FFT/s' num2str(isubject) '_FFT_EGG.mat'];
    end
end

if exist(fname_FFT_EGG,'file') && force ~= 1
    fprintf('Returning FFT\n');
    load(fname_FFT_EGG);
else
    
    % process EGG with Welch method to determine max EGG freq and sensor
    fprintf('FFT not found. Making it now!\n');
    
    if restingstate == 1
        fprintf('Loading restingstate data\n');
        cfg                         = [];
        cfg.continuous              = 'yes';
        cfg.dataset                 = dataset_rs{isubject};
        cfg.channel                 = {'BIO004','BIO005','BIO006','BIO007','BIO008','BIO009','BIO010','BIO011'};
        datapart                    = ft_preprocessing(cfg);
        
        cfg = [];
        cfg.length                  = 200;
        cfg.overlap                 = 0.75;
        data_cut                    = ft_redefinetrial(cfg,datapart);
    else
        for ipart = 1 : 4
            fprintf('Loading datapart (task) %s\n',num2str(ipart));
            cfg                     = [];
            cfg.continuous          = 'yes';
            cfg.dataset             = dataset_task{isubject,ipart};
            cfg.channel             = {'BIO004','BIO005','BIO006','BIO007','BIO008','BIO009','BIO010','BIO011'};
            datapart{ipart}         = ft_preprocessing(cfg);
            
            cfg = [];
            cfg.length              = 200;
            cfg.overlap             = 0.75;
            datapart_cut{ipart}     = ft_redefinetrial(cfg,datapart{ipart});
        end
        
        data_cut = ft_appenddata([],datapart_cut{:});
        clear datapart_cut datapart
    end
    
    cfg                 = [];
    cfg.output          = 'pow';
    cfg.channel         = 'all';
    cfg.method          = 'mtmfft';
    cfg.taper           = 'hann';
    cfg.keeptrials      = 'no';
    cfg.foilim          = [0 0.25];
    cfg.pad             = 1000;
    FFT_EGG             = ft_freqanalysis(cfg, data_cut);
    clear data_cut
    
    % get peak frequency EGG bounds
    low_freq_indx             = find(FFT_EGG.freq > 0.04, 1, 'first');
    high_freq_indx            = find(FFT_EGG.freq < 0.06, 1, 'last');
    
    % find maximum power within bounds for each electrode
    [max_pow,max_freq_indx] = max(FFT_EGG.powspctrm(:,low_freq_indx:high_freq_indx)');
    
    % recover frequency index
    max_freq_indx           = max_freq_indx + low_freq_indx - 1;
    
    % get corresponding frequency
    max_freq                = FFT_EGG.freq(max_freq_indx);
    
    % get electrode with maximum power
    [~,max_chan_indx]         = max(max_pow);
    
    % exceptions
    if isubject == 5
        max_chan_indx = 3;
    end
    if isubject == 26
        max_chan_indx = 6;
    end
    
    % return to datastructure
    FFT_EGG.max_chan     = FFT_EGG.label{max_chan_indx};
    FFT_EGG.max_freq     = max_freq(max_chan_indx);
    
    % save data
    save(fname_FFT_EGG,'FFT_EGG');
end

% plot FFT
% fig = figure;
% plot(FFT_EGG.freq,FFT_EGG.powspctrm);
%saveas(fig,['d:\analysis\WANDER\images\EGG\s' num2str(isubject) '_FFT_EGG.jpg']);

