slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

force = 0;
rootpath = 1;
restingstate = 0;

addpath('D:/analysis/WANDER/scripts/');
% addpath('D:\fieldtrip\fieldtrip.git\trunk');

for isubject = slist
    
    disp(['subject ', num2str(isubject)]);
    fname_SSEF_phase = ['w:\WANDER\data\SSEF\SSEF_phase_' num2str(isubject) '.mat'];
    
    % load directly from file
    temp = load(fname_SSEF_phase,'data_MEG_stim_avg');
    data_MEG_stim_avg{isubject} = temp.data_MEG_stim_avg;
    temp = load(fname_SSEF_phase,'data_MEG_stim_avg_high');
    data_MEG_stim_avg_high{isubject} = temp.data_MEG_stim_avg_high;
    temp = load(fname_SSEF_phase,'data_MEG_stim_avg_low');
    data_MEG_stim_avg_low{isubject} = temp.data_MEG_stim_avg_low;
    temp = load(fname_SSEF_phase,'data_MEG_stim_MI');
    data_MEG_stim_MI{isubject} = temp.data_MEG_stim_MI;
    temp = load(fname_SSEF_phase,'data_MEG_stim_MI_high');
    data_MEG_stim_MI_high{isubject} = temp.data_MEG_stim_MI_high;
    temp = load(fname_SSEF_phase,'data_MEG_stim_MI_low');
    data_MEG_stim_MI_low{isubject} = temp.data_MEG_stim_MI_low;
    
    data_MEG_stim_MI_low{isubject}.fsample  = 1000;
    data_MEG_stim_MI_high{isubject}.fsample = 1000;
    data_MEG_stim_MI{isubject}.fsample      = 1000;
    
    % for some reason it reverts to single trial rather than average data
    % structure - this is a fix
    %     data_MEG_stim_MI_low{isubject}          = ft_timelockanalysis([],data_MEG_stim_MI_low{isubject});
    %     data_MEG_stim_MI_high{isubject}         = ft_timelockanalysis([],data_MEG_stim_MI_high{isubject});
    %     data_MEG_stim_MI{isubject}              = ft_timelockanalysis([],data_MEG_stim_MI{isubject});
    
    % normalize data as percentage of total
    data_MEG_MI_diff{isubject}              = data_MEG_stim_MI{isubject};
    data_MEG_MI_diff{isubject}.avg          = (data_MEG_stim_MI_high{isubject}.avg - data_MEG_stim_MI_low{isubject}.avg) ./ (data_MEG_stim_MI_high{isubject}.avg + data_MEG_stim_MI_low{isubject}.avg);
    %     data_MEG_MI_diff{isubject}.avg          = (data_MEG_stim_MI_high{isubject}.avg - data_MEG_stim_MI_low{isubject}.avg);
    
    data_MEG_stim_MI_low_norm{isubject} = data_MEG_stim_MI_low{isubject};
    data_MEG_stim_MI_high_norm{isubject} = data_MEG_stim_MI_high{isubject};
    data_MEG_stim_MI_low_norm{isubject}.avg      = data_MEG_stim_MI_low{isubject}.avg  ./ (data_MEG_stim_MI_high{isubject}.avg + data_MEG_stim_MI_low{isubject}.avg);
    data_MEG_stim_MI_high_norm{isubject}.avg     = data_MEG_stim_MI_high{isubject}.avg ./ (data_MEG_stim_MI_high{isubject}.avg + data_MEG_stim_MI_low{isubject}.avg);
    %
    data_MEG_stim_diff{isubject}            = data_MEG_stim_avg{isubject};
    data_MEG_stim_diff{isubject}.avg        = (data_MEG_stim_avg_low{isubject}.avg - data_MEG_stim_avg_high{isubject}.avg);
    
end

% calculate rms and time of maximum response
for isubject = slist
    
    % square ERF
    data_MEG_stim_rms{isubject}          = data_MEG_stim_avg{isubject};
    data_MEG_stim_rms{isubject}.avg      = data_MEG_stim_rms{isubject}.avg .^2;
    data_MEG_stim_rms_high{isubject}     = data_MEG_stim_avg_high{isubject};
    data_MEG_stim_rms_high{isubject}.avg = data_MEG_stim_rms_high{isubject}.avg.^2;
    data_MEG_stim_rms_low{isubject}      = data_MEG_stim_avg_low{isubject};
    data_MEG_stim_rms_low{isubject}.avg  = data_MEG_stim_rms_low{isubject}.avg.^2;
    
    % root mean ERF
    cfg = [];
    cfg.avgovertime = 'no';
    data_MEG_stim_rms{isubject}          = ft_selectdata(cfg,data_MEG_stim_rms{isubject});
    data_MEG_stim_rms{isubject}.avg      = sqrt(data_MEG_stim_rms{isubject}.avg);
    data_MEG_stim_rms_high{isubject}     = ft_selectdata(cfg,data_MEG_stim_rms_high{isubject});
    data_MEG_stim_rms_high{isubject}.avg = sqrt(data_MEG_stim_rms_high{isubject}.avg);
    data_MEG_stim_rms_low{isubject}      = ft_selectdata(cfg,data_MEG_stim_rms_low{isubject});
    data_MEG_stim_rms_low{isubject}.avg  = sqrt(data_MEG_stim_rms_low{isubject}.avg);
    
    % diff
    data_MEG_stim_rms_diff{isubject}     = data_MEG_stim_rms_high{isubject};
    %     data_MEG_stim_rms_diff{isubject}.avg = (data_MEG_stim_rms_low{isubject}.avg - data_MEG_stim_rms_high{isubject}.avg) ./ (data_MEG_stim_rms_low{isubject}.avg + data_MEG_stim_rms_high{isubject}.avg) ;
    data_MEG_stim_rms_diff{isubject}.avg = (data_MEG_stim_rms_low{isubject}.avg - data_MEG_stim_rms_high{isubject}.avg);
    
end

for isubject = slist
    
    % time of maximum response
    
    cfg = [];
    cfg.channel = 'MEG*1';
    %     cfg.avgovertime = 'yes';
    cfg.avgoverchan = 'yes';
    temp = ft_selectdata(cfg,data_MEG_stim_rms{isubject});
    %     [~, max_sens_indx{isubject}] = max(temp.avg);
    %     max_sens_indx{isubject} = temp.label(max_sens_indx{isubject});
    % %
    %     cfg = [];
    %     cfg.channel = max_sens_indx{isubject};
    %     temp = ft_selectdata(cfg,data_MEG_stim_rms{isubject});
    [~, max_time_indx] = max(temp.avg);
    max_time(isubject) = temp.time(max_time_indx);
end

figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    
    cfg = [];
    cfg.zlim    = 'maxabs';
    cfg.layout  = 'neuromag306mag';
    cfg.channel = {'MEG*1'};
    cfg.xlim = [max_time(isubject) max_time(isubject)];
    %     cfg.xlim = [max_time(isubject) max_time(isubject)];
    
    ft_topoplotER(cfg,data_MEG_stim_avg{slist(i)});
    title(num2str(isubject));
    %         ft_topoplotER(cfg,data_MEG_stim_rms{slist(i)});
    %     ft_topoplotER(cfg,data_MEG_stim_MI{slist(i)});
    %     ft_topoplotER(cfg,data_MEG_MI_diff{slist(i)});
    %     ft_topoplotER(cfg,data_MEG_MI_diff{slist(i)});
    i = i + 1;
end

% dipole fitting
for isubject = slist
    
    % load headmodel
    [mri_segmented{isubject}, headmodel{isubject}, subject_grid{isubject}, template_grid, mri_realigned{isubject}] = WANDER_grid(isubject,rootpath,0);
    
    % load sensor information
    [dataset_task, ~] = WANDER_subjectinfo;
    hdr{isubject} = ft_read_header(dataset_task{isubject,1});
    
    cfg = [];
    cfg.avgovertime = 'yes';
    cfg.latency = [max_time(isubject) max_time(isubject)];
    data_MEG_stim_maxtime{isubject} = ft_selectdata(cfg,data_MEG_stim_avg{isubject});
    
    % dipole fit
    cfg = [];
    cfg.numdipoles = 1;
    cfg.grid        = subject_grid{isubject};
    %   cfg.grid.pos    = template_grid.pos;
    cfg.vol         = headmodel{isubject};
    cfg.grad        = hdr{isubject}.grad;
    cfg.channel     = 'MEG*1';
    cfg.nonlinear   = 'no';
    source_mag{isubject} = ft_dipolefitting(cfg, data_MEG_stim_maxtime{isubject});
end

template_mri            = ft_read_mri('D:\analysis\WANDER\scripts\colin27_t1_tal_lin.nii');
for isubject = slist
    
    pos_subject{isubject} = source_mag{isubject}.dip.pos;
    ori_subject{isubject} = source_mag{isubject}.dip.mom;
    
    index = find(subject_grid{isubject}.pos(:,1) == source_mag{isubject}.dip.pos(:,1) & subject_grid{isubject}.pos(:,2) == source_mag{isubject}.dip.pos(:,2) & subject_grid{isubject}.pos(:,3) == source_mag{isubject}.dip.pos(:,3));
    pos_template{isubject} = template_grid.pos(index,:);
    ori_template{isubject} = source_mag{isubject}.dip.mom;
    
end


i = 1;
fig = figure;
for isubject = slist
    subplot(5,5,i); hold;
    ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos_template{isubject}, 'plotmarker' ,pos_template{isubject},'markersize',30,'markercolor',[1 1 1])
    
    view(0,90);
    i = i + 1;
    %     title(num2str(isubject));
end

% write positions
% https://www.nitrc.org/docman/view.php/504/1190/BrainNet%20Viewer%20Manual%201.41.pdf
% The node file is defined as an ASCII text file with the suffix ‘node’. In the node file, there
% are 6 columns: columns 1-3 represent node coordinates, column 4 represents node
% colors, column 5 represents node sizes, and the last column represents node labels.
% Please note, a symbol ‘-‘(no ‘’) in column 6 means no labels. The user may put the
% modular information of the nodes into column 4, like ‘1, 2, 3…’ or other information to
% be shown by color. Column 5 could be set as nodal degree, centrality, T-value, etc. to
% emphasize nodal differences by size. You can generate your nodal file according to the
% requirements.

i = 1;
for isubject = slist
    d(i,:) = [pos_template{isubject},1,0.1,isubject];
    i = i + 1;
end

% tabdelimited
dlmwrite('D:\analysis\WANDER\dipoles.node',d,'\t')

% extract timecourses
for isubject = slist
    
    fprintf('Loading subject %d \n',isubject);
    temp = load(['W:\WANDER\data\ERF\s' num2str(isubject) '_cue.mat'],'ERF_cue','ERF_cue_high','ERF_cue_low');
    ERF_cue{isubject}       = temp.ERF_cue;
    ERF_cue_high{isubject}  = temp.ERF_cue_high;
    ERF_cue_low{isubject}   = temp.ERF_cue_low;
    clear temp
    
    fprintf('Cutting subject %d to size \n',isubject);
    cfg = [];
%     cfg.latency = [-0.1 1.5]; % wider for TFR
    
    % only magnetometers, since dipole is based on those
    cfg.channel = 'MEG*1';
    ERF_cue{isubject}       = ft_selectdata(cfg,ERF_cue{isubject});
    ERF_cue_low{isubject}   = ft_selectdata(cfg,ERF_cue_low{isubject});
    ERF_cue_high{isubject}  = ft_selectdata(cfg,ERF_cue_high{isubject});
    
    
    cfg = [];
    cfg.latency = 'all';
    cfg.numdipoles = 1;
    cfg.symmetry = [];
    cfg.nonlinear = 'no';  % use a fixed position
    cfg.gridsearch = 'no';
    cfg.dip.pos = pos_subject{isubject};
    cfg.vol = headmodel{isubject};
    cfg.channel = {'MEG*1'};
    cfg.senstype = 'meg';
    cfg.grad = hdr{isubject}.grad;
    source_all{isubject}    = ft_dipolefitting(cfg, ERF_cue{isubject}); % estimate the amplitude and orientation
    source_high{isubject}   = ft_dipolefitting(cfg, ERF_cue_high{isubject}); % estimate the amplitude and orientation
    source_low{isubject}    = ft_dipolefitting(cfg, ERF_cue_low{isubject}); % estimate the amplitude and orientation
    
    % project to first component
    source_all{isubject}.dip.svd    = svdfft(source_all{isubject}.dip.mom,1);
    source_high{isubject}.dip.svd   = svdfft(source_high{isubject}.dip.mom,1);
    source_low{isubject}.dip.svd    = svdfft(source_low{isubject}.dip.mom,1);
    
end

% put in data matrix
i = 1;
tc      = nan(length(slist),40000);
tc_high = nan(length(slist),40000);
tc_low  = nan(length(slist),40000);

for isubject = slist
    tc(i,1:length(source_all{isubject}.dip.svd)) = source_all{isubject}.dip.svd;
    tc_high(i,1:length(source_high{isubject}.dip.svd)) = source_high{isubject}.dip.svd;
    tc_low(i,1:length(source_low{isubject}.dip.svd)) = source_low{isubject}.dip.svd;  
    i = i + 1;  
end

% orient according to correlation with mean
latency = [1500:2000];
temp_avg = mean(tc,1);
for iteration = 1 : 100
    disp(['orienting iteration',num2str(iteration)]);
    for isubject = 1 : size(tc,1)
        if corr(tc(isubject,latency)',temp_avg(latency)','rows','complete') < 0
            flip(isubject) = 1;
            tc(isubject,:) = -tc(isubject,:);
            tc_high(isubject,:) = -tc_high(isubject,:);
            tc_low(isubject,:) = -tc_low(isubject,:); 
        end
    end
    temp_avg = mean(tc,1);
end

figure; imagesc(tc(:,latency))
xticklabels = latency(1:25:end)/1000-1.5;
xticks = linspace(1, length(latency), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

latency = 1500:3000;
tc_avg = nanmean(tc(:,latency));
tc_avg_high = nanmean(tc_high(:,latency));
tc_avg_low = nanmean(tc_low(:,latency));
time = latency/1000-1.5-0.033;

figure; hold;
plot(time,tc_avg);
axis tight;

figure; hold;
plot(time,tc_avg_high,'b');
plot(time,tc_avg_low,'r');
axis tight


%% FFT

latency = 1633:3000;

% tc_avg      = nanmean(tc(:,latency));
% tc_avg_high = nanmean(tc_high(:,latency));
% tc_avg_low  = nanmean(tc_low(:,latency));

Fs = 1000;
L = length(latency);
f = Fs*(0:(L/2))/L;

% FFT_avg               = fft(tc_avg);
% FFT_avg_high          = fft(tc_avg_high);
% FFT_avg_low           = fft(tc_avg_low);

FFT_avg               = fft(tc(:,latency)');
FFT_avg_high          = fft(tc_high(:,latency)');
FFT_avg_low           = fft(tc_low(:,latency)');

FFT_avg               = abs(FFT_avg/L);
FFT_avg               = FFT_avg(1:L/2+1,:);
FFT_avg(2:end-1,:)    = 2*FFT_avg(2:end-1,:);

FFT_avg_high          = abs(FFT_avg_high/L);
FFT_avg_high          = FFT_avg_high(1:L/2+1,:);
FFT_avg_high(2:end-1,:) = 2*FFT_avg_high(2:end-1,:);

FFT_avg_low           = abs(FFT_avg_low/L);
FFT_avg_low           = FFT_avg_low(1:L/2+1,:);
FFT_avg_low(2:end-1,:)  = 2*FFT_avg_low(2:end-1,:);

[~, ilow] = find(f>1,1,'first');
[~, ihigh] = find(f>30,1,'first');
    
figure;
plot(f(ilow:ihigh),FFT_avg(ilow:ihigh,:)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure;
plot(f(ilow:ihigh),mean(FFT_avg(ilow:ihigh,:),2)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


figure; hold;
plot(f(ilow:ihigh),mean(FFT_avg_high(ilow:ihigh,:),2),'b');
plot(f(ilow:ihigh),mean(FFT_avg_low(ilow:ihigh,:),2),'r');
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% relative difference
FFT_avg_diff = (FFT_avg_low - FFT_avg_high) ./ (FFT_avg_low + FFT_avg_high);

figure; hold;
plot(f(ilow:ihigh),FFT_avg_diff(ilow:ihigh,:));
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure; hold;
plot(f(ilow:ihigh),mean(FFT_avg_diff(ilow:ihigh,:),2));
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')



%% FFT welch method

% put in FT format
for isubject = slist
    data{isubject}              = rmfield(ERF_cue{isubject},{'cfg','avg','dof','label'});
    data{isubject}.avg          = source_all{isubject}.dip.svd;
    data{isubject}.label        = {'dipole'};
    
    data_high{isubject}         = rmfield(ERF_cue_high{isubject},{'cfg','avg','dof','label'});
    data_high{isubject}.avg     = source_high{isubject}.dip.svd;
    data_high{isubject}.label   = {'dipole'};   

    data_low{isubject}          = rmfield(ERF_cue_low{isubject},{'cfg','avg','dof','label'});
    data_low{isubject}.avg      = source_low{isubject}.dip.svd;
    data_low{isubject}.label    = {'dipole'};  
    
    cfg                 = [];
    cfg.output          = 'pow';
    cfg.method          = 'mtmconvol';
    cfg.keeptrials      = 'no';
    cfg.taper           = 'hanning';
    cfg.foi             = [1:30];
    cfg.t_ftimwin       = ones(size(cfg.foi))*1;
    cfg.toi             = [-1.5 : 0.01 : 10];
    TFR{isubject}       = ft_freqanalysis(cfg, data{isubject});
    TFR_high{isubject}  = ft_freqanalysis(cfg, data_high{isubject});
    TFR_low{isubject}   = ft_freqanalysis(cfg, data_low{isubject});   
    TFR_diff{isubject}  = TFR_high{isubject};
    TFR_diff{isubject}.powspctrm = (TFR_low{isubject}.powspctrm - TFR_high{isubject}.powspctrm) ./ (TFR_low{isubject}.powspctrm + TFR_high{isubject}.powspctrm);

end

TFR_avg = ft_freqgrandaverage([],TFR{slist});
TFR_diff_avg = ft_freqgrandaverage([],TFR_diff{slist});

figure; 
cfg = [];
cfg.colormap = parula(1000);
ft_singleplotTFR(cfg,TFR_avg);

figure; 
cfg = [];
ft_singleplotTFR(cfg,TFR_diff_avg);



latency = 1633:3000;

% tc_avg      = nanmean(tc(:,latency));
% tc_avg_high = nanmean(tc_high(:,latency));
% tc_avg_low  = nanmean(tc_low(:,latency));

Fs = 1000;
L = length(latency);
f = Fs*(0:(L/2))/L;

% FFT_avg               = fft(tc_avg);
% FFT_avg_high          = fft(tc_avg_high);
% FFT_avg_low           = fft(tc_avg_low);

FFT_avg               = fft(tc(:,latency)');
FFT_avg_high          = fft(tc_high(:,latency)');
FFT_avg_low           = fft(tc_low(:,latency)');

FFT_avg               = abs(FFT_avg/L);
FFT_avg               = FFT_avg(1:L/2+1,:);
FFT_avg(2:end-1,:)    = 2*FFT_avg(2:end-1,:);

FFT_avg_high          = abs(FFT_avg_high/L);
FFT_avg_high          = FFT_avg_high(1:L/2+1,:);
FFT_avg_high(2:end-1,:) = 2*FFT_avg_high(2:end-1,:);

FFT_avg_low           = abs(FFT_avg_low/L);
FFT_avg_low           = FFT_avg_low(1:L/2+1,:);
FFT_avg_low(2:end-1,:)  = 2*FFT_avg_low(2:end-1,:);

[~, ilow] = find(f>1,1,'first');
[~, ihigh] = find(f>30,1,'first');
    
figure;
plot(f(ilow:ihigh),FFT_avg(ilow:ihigh,:)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure;
plot(f(ilow:ihigh),mean(FFT_avg(ilow:ihigh,:),2)) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


figure; hold;
plot(f(ilow:ihigh),mean(FFT_avg_high(ilow:ihigh,:),2),'b');
plot(f(ilow:ihigh),mean(FFT_avg_low(ilow:ihigh,:),2),'r');
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% relative difference
FFT_avg_diff = (FFT_avg_low - FFT_avg_high) ./ (FFT_avg_low + FFT_avg_high);

figure; hold;
plot(f(ilow:ihigh),FFT_avg_diff(ilow:ihigh,:));
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure; hold;
plot(f(ilow:ihigh),mean(FFT_avg_diff(ilow:ihigh,:),2));
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')






time = latency/1000-1.5-0.033;

figure; hold;
plot(time,tc_avg);
axis tight;

figure; hold;
plot(time,tc_avg_high,'b');
plot(time,tc_avg_low,'r');
axis tight







s = svdfft(source_all{isubject}.dip.mom,1);
figure; hold;
plot(s(1:1000));
plot(source_all{isubject}.dip.mom(:,1:1000)')

data_MEG_stim_GA        = ft_timelockgrandaverage([],data_MEG_stim_avg{slist});
data_MEG_stim_high_GA   = ft_timelockgrandaverage([],data_MEG_stim_avg_high{slist});
data_MEG_stim_low_GA    = ft_timelockgrandaverage([],data_MEG_stim_avg_low{slist});
data_MEG_stim_diff_GA   = ft_timelockgrandaverage([],data_MEG_stim_diff{slist});
data_MEG_MI_GA          = ft_timelockgrandaverage([],data_MEG_stim_MI{slist});
data_MEG_MI_high_GA     = ft_timelockgrandaverage([],data_MEG_stim_MI_high{slist});
data_MEG_MI_low_GA      = ft_timelockgrandaverage([],data_MEG_stim_MI_low{slist});
data_MEG_MI_diff_GA     = ft_timelockgrandaverage([],data_MEG_MI_diff{slist});
data_MEG_rms_GA         = ft_timelockgrandaverage([],data_MEG_stim_rms{slist});
data_MEG_rms_high_GA    = ft_timelockgrandaverage([],data_MEG_stim_rms_high{slist});
data_MEG_rms_low_GA     = ft_timelockgrandaverage([],data_MEG_stim_rms_low{slist});
data_MEG_rms_diff_GA    = ft_timelockgrandaverage([],data_MEG_stim_rms_diff{slist});

% plot timecourses
cfg         = [];
cfg.zlim    = 'maxabs';
cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEG*1';
figure; ft_multiplotER(cfg,data_MEG_stim_avg{25});

% shows no clear difference, perhaps because these are averages of absolute
figure; ft_multiplotER(cfg,data_MEG_stim_GA);
figure; ft_multiplotER(cfg,data_MEG_stim_high_GA,data_MEG_stim_low_GA);
figure; ft_multiplotER(cfg,data_MEG_stim_diff_GA);

cfg.layout  = 'neuromag306cmb';
cfg.channel = 'MEG*3';

% shows no clear difference, perhaps because these are averages of absolute
figure; ft_multiplotER(cfg,data_MEG_stim_high_GA,data_MEG_stim_low_GA);

% plot topography
cfg = [];
cfg.zlim    = 'maxabs';
cfg.layout  = 'neuromag306mag';
cfg.channel = {'MEG*1'};

% topography
figure; ft_topoplotER(cfg,data_MEG_stim_GA);

cfg.xlim = [0 : 0.01 : 0.06];
figure; ft_topoplotER(cfg,data_MEG_stim_diff_GA);

figure; ft_topoplotER(cfg,data_MEG_MI_GA);
figure; ft_topoplotER(cfg,data_MEG_MI_diff_GA);
figure; ft_topoplotER(cfg,data_MEG_MI_high_GA);
figure; ft_topoplotER(cfg,data_MEG_MI_low_GA);

figure; ft_topoplotER(cfg,data_MEG_rms_GA);
figure; ft_topoplotER(cfg,data_MEG_rms_diff_GA);
figure; ft_topoplotER(cfg,data_MEG_rms_high_GA);
figure; ft_topoplotER(cfg,data_MEG_rms_low_GA);




cfg         = [];
cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEGMAG';
cfg.zlim    = 'absmax';
figure; ft_multiplotER(cfg,data_MEG_stim_bin{isubject}{:});

cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEGMAG';
cfg.zlim    = 'absmax';
cfg.colorbar = 'yes';
cfg.parameter = 'rms';
figure;
ft_multiplotER(cfg,data_MEG_stim_binhist{isubject});
ft_multiplotER(cfg,data_MEG_stim_binhist_low{isubject},data_MEG_stim_binhist_high{isubject});


cfg         = [];
cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEG*1';
cfg.zlim    = [-0.02 0.02];
figure; ft_topoplotER(cfg,data_MEG_stim_MI_high{isubject},data_MEG_stim_MI_low{isubject});
cfg.zlim    = 'absmax';
figure; ft_topoplotER(cfg,data_MEG_stim_MI_diff{isubject});


% STATISTICS

load('neuromag306mag_neighb_last');

cfg = [];
cfg.channel             = 'MEG*1';
cfg.neighbours          = neighbours;
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 1000;
cfg.alpha               = 0.05;
cfg.clusteralpha        = 0.05;
cfg.tail                = 0;
cfg.design(1,:)         = [1:length(slist) 1:length(slist)];
cfg.design(2,:)         = [ones(size(slist)) ones(size(slist))*2];
cfg.uvar                = 1;
cfg.ivar                = 2;
cfg.parameter           = 'avg';
stat_MI                 = ft_timelockstatistics(cfg, data_MEG_stim_MI_high{slist}, data_MEG_stim_MI_low{slist});
stat_MI                 = ft_timelockstatistics(cfg, data_MEG_stim_MI_high_norm{slist}, data_MEG_stim_MI_low_norm{slist});


stat_rms                = ft_timelockstatistics(cfg, data_MEG_stim_rms_low{slist}, data_MEG_stim_rms_high{slist});


cfg         = [];
cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEG*1';
cfg.highlight = 'on';
cfg.parameter = 'stat';
cfg.highlightchannel = {stat_MI.label{stat_MI.mask}};
% cfg.zlim    = [-0.02 0.02];
figure; ft_topoplotER(cfg,stat_MI);

cfg.highlightchannel = {stat_rms.label{stat_rms.mask}};

cfg.zlim    = 'maxabs';

figure; ft_topoplotER(cfg,stat_rms);

stat_MI = ft_combineplanar([],stat_MI);
cfg         = [];
cfg.layout  = 'neuromag306grad';
cfg.channel = 'MEG*1';
cfg.highlight = 'on';
cfg.parameter = 'stat';
cfg.highlightchannel = {stat_MI.label{stat_MI.mask}};
% cfg.zlim    = [-0.02 0.02];
figure; ft_topoplotER(cfg,stat_MI);


