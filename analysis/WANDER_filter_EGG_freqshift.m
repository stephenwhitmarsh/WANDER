function [data_EGG] = WANDER_filter_EGG_freqshift(isubject,force,rootpath,restingstate)

[dataset_task, dataset_rs] = WANDER_subjectinfo;

if restingstate == 1
    if rootpath == 1
        fname_EGGfilt = ['d:\analysis\WANDER\data\raw\s' num2str(isubject) '_EGG_rs_freqshift.mat'];
    else
        fname_EGGfilt = ['/home/swhitmarsh/WANDER/data/raw/s' num2str(isubject) '_EGG_rs_freqshift.mat'];
    end
else
    if rootpath == 1
        fname_EGGfilt = ['d:\analysis\WANDER\data\raw\s' num2str(isubject) '_EGG_freqshift.mat'];
    else
        fname_EGGfilt = ['/home/swhitmarsh/WANDER/data/raw/s' num2str(isubject) '_EGG_freqshift.mat'];
    end
end

if exist(fname_EGGfilt,'file') && force~=1
    fprintf('Returning filtered EGG data\n');
    load(fname_EGGfilt);
else
    fprintf('Filtered EGG data not found. Making it now!\n');
    FFT_EGG = WANDER_FFT_EGG(isubject,0,rootpath,restingstate);
    % plus/minus 0.015 Hz, in steps of 0.005 Hz)
    EGG_freq_offset = [-0.005 : 0.005 : 0.03];
    
    if restingstate == 1
        fprintf('Loading datafile \n');
        
        % load EGG channels
        cfg                 = [];
        cfg.continuous      = 'yes';
        cfg.dataset         = dataset_rs{isubject};
        cfg.channel         = FFT_EGG.max_chan;
        EGG                 = ft_preprocessing(cfg);
        
        for ifreq = 1:length(EGG_freq_offset)
            % set up filter
            srate               = 1000;
            center_frequency    = FFT_EGG.max_freq + EGG_freq_offset(ifreq); % depends on data
            bandwidth           = 0.02;
            transition_width    = 0.15;
            nyquist             = srate/2;
            ffreq(1)            = 0;
            ffreq(2)            = (1-transition_width)*(center_frequency-bandwidth);
            ffreq(3)            = (center_frequency-bandwidth);
            ffreq(4)            = (center_frequency+bandwidth);
            ffreq(5)            = (1+transition_width)*(center_frequency+bandwidth);
            ffreq(6)            = nyquist;
            ffreq               = ffreq/nyquist;
            fOrder              = 3; % in cycles changed from 7
            filterOrder         = fOrder*fix(srate/(center_frequency - bandwidth)); %in samples
            idealresponse       = [ 0 0 1 1 0 0 ];
            filterweights       = fir2(filterOrder,ffreq,idealresponse);
            disp(size(EGG.trial{1})/size(filterweights));
            
            % filter
            EGG_filt = EGG;
            disp('Filtering EGG - this will take some time');
            EGG_filt.trial{1}   = filtfilt(filterweights,1,EGG.trial{1});
            EGG_filt.label      = {'filtered'};
            
            % analytic EGG
            EGG_phase = EGG;
            EGG_phase.trial{1}   = angle(hilbert(EGG_filt.trial{1}'))';
            EGG_phase.label      = {'phase'};
            
            % append filtered EGG and EGG phase to rest of data
            data_EGG{ifreq} = ft_appenddata([],EGG,EGG_phase,EGG_filt);
            clear EGG_filt EGG_phase
        end
    else
        for ipart = 1 : 4
            fprintf('Loading datafile %d\n',ipart);
            
            % load EGG channels
            cfg                 = [];
            cfg.continuous      = 'yes';
            cfg.dataset         = dataset_task{isubject,ipart};
            cfg.channel         = FFT_EGG.max_chan;
            EGG                 = ft_preprocessing(cfg);
            for ifreq = 1:length(EGG_freq_offset)
                
                % set up filter
                srate               = 1000;
                center_frequency    = FFT_EGG.max_freq + EGG_freq_offset(ifreq); % depends on data
                bandwidth           = 0.02;
                transition_width    = 0.15;
                nyquist             = srate/2;
                ffreq(1)            = 0;
                ffreq(2)            = (1-transition_width)*(center_frequency-bandwidth);
                ffreq(3)            = (center_frequency-bandwidth);
                ffreq(4)            = (center_frequency+bandwidth);
                ffreq(5)            = (1+transition_width)*(center_frequency+bandwidth);
                ffreq(6)            = nyquist;
                ffreq               = ffreq/nyquist;
                fOrder              = 3; % in cycles
                filterOrder         = fOrder*fix(srate/(center_frequency - bandwidth)); %in samples
                idealresponse       = [ 0 0 1 1 0 0 ];
                filterweights       = fir2(filterOrder,ffreq,idealresponse);
                
                % filter
                EGG_filt = EGG;
                disp('Filtering EGG - this will take some time');
                EGG_filt.trial{1}   = filtfilt(filterweights,1,EGG.trial{1});
                EGG_filt.label      = {'filtered'};
                
                % analytic EGG
                EGG_phase = EGG;
                EGG_phase.trial{1}   = angle(hilbert(EGG_filt.trial{1}'))';
                EGG_phase.label      = {'phase'};
                
                % append filtered EGG and EGG phase to rest of data
                data_EGG{ifreq}{iparam} = ft_appenddata([],EGG,EGG_phase,EGG_filt);
                clear EGG_filt EGG_phase
            end
        end
    end
    save(fname_EGGfilt,'data_EGG','-v7.3');
end

% FFT_EGG = WANDER_FFT_EGG(isubject,0);
%
% for ipart = 1 : nrofsets(isubject)
%     fig = figure;
%     x = [1 100000];
%     subplot(3,1,1); plot(data_EGG{ipart}.time{1}(x(1):x(2)),data_EGG{ipart}.trial{1}(strcmp([FFT_EGG.max_chan], data_EGG{ipart}.label),x(1):x(2)));             title('raw');
%     subplot(3,1,2); plot(data_EGG{ipart}.time{1}(x(1):x(2)),data_EGG{ipart}.trial{1}(strcmp([FFT_EGG.max_chan '_filt'], data_EGG{ipart}.label),x(1):x(2)));     title('filtered');
%     subplot(3,1,3); plot(data_EGG{ipart}.time{1}(x(1):x(2)),data_EGG{ipart}.trial{1}(strcmp([FFT_EGG.max_chan '_phase'], data_EGG{ipart}.label),x(1):x(2)));    title('phase');
%     saveas(fig,['d:\analysis\WANDER\images\EGG\s' num2str(isubject) '_EGG_phase_' num2str(ipart) '.jpg']);
% end
