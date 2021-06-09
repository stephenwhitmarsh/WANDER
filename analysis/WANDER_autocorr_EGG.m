rootpath = 1;
restingstate = 1;

clear EGGac* data*

for isubject = 1 : 26
    
%     resting state
    FFT_EGG_rs{isubject}            = WANDER_FFT_EGG(isubject,0,rootpath,1);
    data_EGG_rs{isubject}           = WANDER_filter_EGG_param(isubject,0,rootpath,1);
    
    for iparam = 1 : 5 
        cfg = [];
        data_EGG_rs_demeaned{isubject} = data_EGG_rs{isubject};
        data_EGG_rs_demeaned{isubject}{iparam}.trial{1}(3,:)  = data_EGG_rs{isubject}{iparam}.trial{1}(3,:) - mean(data_EGG_rs{isubject}{iparam}.trial{1}(3,:));
        
        [EGGac_rs{iparam}(isubject,:) lags] = xcorr(data_EGG_rs_demeaned{isubject}{iparam}.trial{1}(3,:),10*60*1000,'coeff');
    end
    %
%     %task
%     FFT_EGG_task{isubject}            = WANDER_FFT_EGG(isubject,0,rootpath,0);
%     data_EGG_task{isubject}           = WANDER_filter_EGG(isubject,0,rootpath,0);
%     for iblock = 1 : 4
%         
%         cfg = [];
%         cfg.channel                             = {[FFT_EGG_task{isubject}.max_chan '_filt'],'filtered'};
%         data_EGG_task{isubject}{iblock}           = ft_selectdata(cfg,data_EGG_task{isubject}{iblock});
%         data_EGG_task_demeaned{isubject}{iblock}  = data_EGG_task{isubject}{iblock}.trial{1} - mean(data_EGG_task{isubject}{iblock}.trial{1});
%         
%         [EGGac_task(isubject,iblock,:) lags] = xcorr(data_EGG_task_demeaned{isubject}{iblock},10*60*1000,'coeff');
%     end
end

for iparam = 1 : 5
    EGGac_rs_avg{iparam} = squeeze(mean(EGGac_rs{iparam},1));
end

figure; hold;
for iparam = 1 : 5
    plot(lags,EGGac_rs_avg{iparam}); axis tight;
end
legend({'0.005','0.01','0.015','0.02','0.025'});

set(gcf, 'PaperSize', [8.5 11]*3)
set(gcf, 'paperposition', [0 0 8.5 11]*3);
print -dpdf 'd:\analysis\WANDER\images\autocorrelation_EGG_filter_halfwidth_averaged.pdf'

fig = figure;
% set(fig, 'PaperSize', [8.5 11]*3)    % Same, but for PDF output
i = 1;
for isubject = MI_index_diff(1:end-1)
% for isubject = 1 : 26
    subplot(5,5,i);
    hold;
    for iparam = [1, 5]
        plot(lags,EGGac_rs{iparam}(isubject,:)); axis tight;
    end
    i = i + 1;
    title(['S' num2str(isubject) ', MI diff: ' num2str(MI_value_diff(i))]);
    
end

fig = figure;
subplot(3,1,1); hold;
for iparam = [1, 5]
    plot(lags,EGGac_rs{iparam}(21,:),'linewidth',1.3); axis tight;
end
subplot(3,1,2); hold;
for iparam = [1, 5]
    plot(lags,EGGac_rs{iparam}(16,:),'linewidth',1.3); axis tight;
end
subplot(3,1,3); hold;
for iparam = [1, 5]
    plot(lags,EGGac_rs{iparam}(11,:),'linewidth',1.3); axis tight;
end

subplot(6,5,27);
hold;
for iparam = 1 : 5
    plot(lags,EGGac_rs{iparam}(isubject,:)); axis tight;
end
legend({'0.005','0.01','0.015','0.02','0.025'});

set(gcf, 'PaperSize', [8.5 11]*3)
set(gcf, 'paperposition', [0 0 8.5 11]*3);
print -dpdf 'd:\analysis\WANDER\images\autocorrelation_EGG_filter_halfwidth.pdf'

%000000
[dataset_task, dataset_rs] = WANDER_subjectinfo;


for isubject = [21,11]
    for iparam = 1 : 5
        cfg = [];
        cfg.length                  = 200;
        cfg.overlap                 = 0.75;
        data_cut                    = ft_redefinetrial(cfg,data_EGG_rs{isubject}{iparam});
        
        cfg                 = [];
        cfg.output          = 'pow';
        cfg.channel         = 'all';
        cfg.method          = 'mtmfft';
        cfg.taper           = 'hann';
        cfg.keeptrials      = 'no';
        cfg.foilim          = [0 0.25];
        cfg.pad             = 1024;
        FFT_EGG_allchan{isubject}{iparam} = ft_freqanalysis(cfg, data_cut);
        clear data_cut
    end
end


bandwidth_param = [0.005 : 0.005 : 0.025];

for isubject = [21,11]
    fig = figure;
    subplot(2,1,1);
    hold;
    for iparam = [1, 5]
        plot(FFT_EGG_allchan{isubject}{iparam}.freq, FFT_EGG_allchan{isubject}{iparam}.powspctrm(3,:),'linewidth',2);
        plot(FFT_EGG_allchan{isubject}{iparam}.freq, FFT_EGG_allchan{isubject}{iparam}.powspctrm(1,:),'k:','linewidth',2);
        axis tight        
        xlim([0.0 0.2]);
    end
    
    subplot(2,1,2);
    hold;
    for iparam = [1, 5]
        bandwidth           = bandwidth_param(iparam);
        center_frequency    = FFT_EGG_rs{isubject}.max_freq;
        transition_width    = 0.15;
        fOrder              = 3;
        srate = 1000;
        nyquist = srate/2;
        ffreq(1)            = 0;
        ffreq(2)            = (1-transition_width)*(center_frequency-bandwidth);
        ffreq(3)            = (center_frequency-bandwidth);
        ffreq(4)            = (center_frequency+bandwidth);
        ffreq(5)            = (1+transition_width)*(center_frequency+bandwidth);
        ffreq(6)            = nyquist;
        ffreq               = ffreq/nyquist;
        filterOrder         = fOrder*fix(srate/(center_frequency - bandwidth)); %in samples
        idealresponse       = [ 0 0 1 1 0 0 ];
        filterweights       = fir2(filterOrder,ffreq,idealresponse);
        y = abs(fft(filterweights));
        
%         subplot(8,2,i*2-1);
        plot(linspace(0,nyquist,length(y)/2),y(1:size(y,2)/2),'linewidth',2); 
        axis tight;
        xlim([0.0 0.2]);
        
    end
    legend({'0.005','0.025'});

%     legend({'0.005','0.01','0.015','0.02','0.025'});
%     set(gcf, 'PaperSize', [8.5 11]*1)
%     set(gcf, 'paperposition', [0 0 8.5 11]*1);
%     saveas(fig,['d:\analysis\WANDER\images\filtered_data_EGG' num2str(isubject) '_zoomed'],'pdf');
end





















figure;
plot(EGGac_task_avg')

figure;
plot(lags, mean(EGGac_rs)); axis tight; hold;
plot(lags, mean(EGGac_task_avg));

figure;
plot(lags, mean(abs(EGGac_rs))); axis tight; hold;
plot(lags, mean(abs(EGGac_task_avg)));


figure;
for isubject = 1 : 26
subplot(6,5,isubject);
plot(lags,EGGac_rs(isubject,:)); axis tight; hold;
plot(lags,EGGac_task_avg(isubject,:)); axis tight
end

figure;
for isubject = 1 : 26
subplot(6,5,isubject);
plot(lags,EGGac_task_avg(isubject,:)); axis tight
end
        
        
        