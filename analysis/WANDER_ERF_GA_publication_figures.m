function [TFR_median_split] = WANDER_ERF_GA(force,timing,rootpath)

if ~(strcmp(timing,'cue') || strcmp(timing,'probe'))
    fprintf('Use "cue" or "probe" as third argument.\n');
    return
end
addpath('D:\fieldtrip-master_21092017\fieldtrip-master\');
addpath('D:/analysis/WANDER/scripts/');
ft_defaults

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD
timing = 'cue';

for isubject = slist
    fprintf('Loading subject %d \n',isubject);    
    temp = load(['W:\WANDER\data\ERF\s' num2str(isubject) '_' timing '.mat'],'ERF_cue','ERF_cue_high','ERF_cue_low');
    ERF_cue{isubject}       = temp.ERF_cue;
    ERF_cue_high{isubject}  = temp.ERF_cue_high;
    ERF_cue_low{isubject}   = temp.ERF_cue_low;
    
    fprintf('Cutting subject %d to size \n',isubject);    
    cfg = [];
    cfg.latency = [-0.1 1.5]; % wider for TFR
    
    cfg.channel = 'MEG*';
    ERF_cue{isubject}       = ft_selectdata(cfg,ERF_cue{isubject});   
    ERF_cue_low{isubject}   = ft_selectdata(cfg,ERF_cue_low{isubject});
    ERF_cue_high{isubject}  = ft_selectdata(cfg,ERF_cue_high{isubject});
end

% for filtered ERF, so not matter we already combine
clear temp
for isubject = slist
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 30; 
    ERF_cue_filt{isubject}       = ft_preprocessing(cfg,ERF_cue{isubject});
    ERF_cue_low_filt{isubject}   = ft_preprocessing(cfg,ERF_cue_low{isubject});
    ERF_cue_high_filt{isubject}  = ft_preprocessing(cfg,ERF_cue_high{isubject});  
    
    cfg = [];
    ERF_cue_filt{isubject}       = ft_combineplanar(cfg,ERF_cue_filt{isubject});
    ERF_cue_low_filt{isubject}   = ft_combineplanar(cfg,ERF_cue_low_filt{isubject});
    ERF_cue_high_filt{isubject}  = ft_combineplanar(cfg,ERF_cue_high_filt{isubject});       
end

% grand averages for plotting - have to reload for TFR
i = 1; clear avg_cue*
for isubject = slist
    cfg                         = [];
    cfg.latency                 = [-0.5 1.5];
    temp                        = ft_selectdata(cfg,ERF_cue_filt{isubject});
    avg_cue_filt(i,:,:)         = temp.avg;
    temp                        = ft_selectdata(cfg,ERF_cue_high_filt{isubject});
    avg_cue_high_filt(i,:,:)    = temp.avg;
    temp                        = ft_selectdata(cfg,ERF_cue_low_filt{isubject});
    avg_cue_low_filt(i,:,:)     = temp.avg;
    i = i + 1;
end

cfg = [];
cfg.latency                     = [-0.5 1.5];
ERF_cue_filt_GA                 = ft_timelockgrandaverage(cfg,ERF_cue_filt{slist});
ERF_cue_high_filt_GA            = ft_timelockgrandaverage(cfg,ERF_cue_high_filt{slist});
ERF_cue_low_filt_GA             = ft_timelockgrandaverage(cfg,ERF_cue_low_filt{slist});
ERF_cue_filt_GA.std             = squeeze(std(avg_cue_filt,1));
ERF_cue_high_filt_GA.std        = squeeze(std(avg_cue_high_filt,1));
ERF_cue_low_filt_GA.std         = squeeze(std(avg_cue_low_filt,1));

%% get highest pos. and neg. deflections at FIRST COMPONENT

peak_latency            = [0 0.100];

cfg                     = [];
cfg.latency             = peak_latency;
cfg.avgovertime         = 'yes';
cfg.avgovertime         = 'no';

cfg.channel             = 'MEG*3';
temp                    = ft_selectdata(cfg,ERF_cue_filt_GA);

% max timepoint
[Y_time,I_time]         = max(mean(abs(temp.avg)));
time_max_cmb            = temp.time(I_time);

% max sensor at max timepoint
[Y_sens,I_sens]         = max(temp.avg(:,I_time));
chan_pos_cmb            = ft_channelselection(I_sens, temp.label);

cfg.channel             = 'MEG*1';
temp                    = ft_selectdata(cfg,ERF_cue_filt_GA);

% max/min timepoint
[Y_time,I_time]         = max(mean(abs(temp.avg)));
time_max_mag            = temp.time(I_time);

% max sensor at max timepoint
[Y_sens,I_sens]         = max(temp.avg(:,I_time));
chan_pos_mag            = ft_channelselection(I_sens, temp.label);

% min sensor at max timepoint
[Y_sens,I_sens]         = min(temp.avg(:,I_time));
chan_neg_mag            = ft_channelselection(I_sens, temp.label);

% select data only of max min sensors
cfg = [];
cfg.parameter                       = {'avg','std'};
cfg.avgoverchan                     = 'yes';
cfg.channel                         = chan_neg_mag;
ERF_cue_low_filt_sel_neg_GA_mag     = ft_selectdata(cfg,ERF_cue_low_filt_GA);
ERF_cue_high_filt_sel_neg_GA_mag    = ft_selectdata(cfg,ERF_cue_high_filt_GA);
cfg.channel                         = chan_pos_mag;
ERF_cue_low_filt_sel_pos_GA_mag     = ft_selectdata(cfg,ERF_cue_low_filt_GA);
ERF_cue_high_filt_sel_pos_GA_mag    = ft_selectdata(cfg,ERF_cue_high_filt_GA);
cfg.channel                         = chan_pos_cmb;
ERF_cue_low_filt_sel_pos_GA_cmb     = ft_selectdata(cfg,ERF_cue_low_filt_GA);
ERF_cue_high_filt_sel_pos_GA_cmb    = ft_selectdata(cfg,ERF_cue_high_filt_GA);


%% FOR ARTICLE: plot overview with ERF pos, neg and topo FIRST COMPONENT

% Fig 1a
fig = figure; 
cfg                     = [];
cfg.layout              = 'neuromag306mag';
cfg.xlim                = [time_max_mag time_max_mag];
cfg.zlim                = 'maxabs';
cfg.highlightchannel    = chan_pos_mag; % positive channel is largest (absolute)
cfg.highlight           = 'yes';
cfg.highlightsymbol     = '.';
cfg.highlightcolor      = [0 0 0];
cfg.highlightsize       = 20;
cfg.marker              = 'off';
cfg.comment             = 'no';
cfg.colorbar            = 'yes';
cfg.colormap            = colormap(jet(1024));
% cfg.gridscale           = 256;
ft_topoplotER(cfg,ERF_cue_filt_GA);
title('');
% 
% set(fig,'PaperSize',[11*2 8.5*2 ]);
% set(fig,'PaperPosition',[0 0 11*4 8.5*2]);
% set(fig,'PaperOrientation','landscape');
% set(fig,'Position',[50 50 1200 800]);

savefig('w:\WANDER\images\FIGa_onset_topo');

% saveas(fig,'w:\WANDER\images\FIGa_onset_topo.eps','eps2c');
% saveas(fig,'w:\WANDER\images\FIGa_onset_topo.png','png');

% Fig 1b
fig = figure; 
t = 1; 

% patch([ERF_cue_low_filt_sel_pos_GA_mag.time, ERF_cue_low_filt_sel_pos_GA_mag.time(end:-1:1)],[ERF_cue_low_filt_sel_pos_GA_mag.avg+ERF_cue_low_filt_sel_pos_GA_mag.std ./sqrt(22)*t, ERF_cue_low_filt_sel_pos_GA_mag.avg(end:-1:1)-ERF_cue_low_filt_sel_pos_GA_mag.std(end:-1:1) ./sqrt(22)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');

hold;
patch([ERF_cue_low_filt_sel_pos_GA_mag.time, ERF_cue_low_filt_sel_pos_GA_mag.time(end:-1:1)],[ERF_cue_low_filt_sel_pos_GA_mag.avg+ERF_cue_low_filt_sel_pos_GA_mag.std ./sqrt(22)*t, ERF_cue_low_filt_sel_pos_GA_mag.avg(end:-1:1)-ERF_cue_low_filt_sel_pos_GA_mag.std(end:-1:1) ./sqrt(21)*t],      [0   0  0.3 ],'LineStyle','None');
patch([ERF_cue_high_filt_sel_pos_GA_mag.time, ERF_cue_high_filt_sel_pos_GA_mag.time(end:-1:1)],[ERF_cue_high_filt_sel_pos_GA_mag.avg+ERF_cue_high_filt_sel_pos_GA_mag.std ./sqrt(22)*t, ERF_cue_high_filt_sel_pos_GA_mag.avg(end:-1:1)-ERF_cue_high_filt_sel_pos_GA_mag.std(end:-1:1) ./sqrt(21)*t],[0.3 0  0.3 ],'LineStyle','None');
ax = axis;

patch([time_max_mag-0.01 time_max_mag+0.01 time_max_mag+0.01 time_max_mag-0.01],[ax(3) ax(3) ax(4) ax(4)],[0.5 0.5 0.5],'LineStyle','None');

line(ERF_cue_low_filt_sel_pos_GA_mag.time,  ERF_cue_low_filt_sel_pos_GA_mag.avg, 'linewidth',1,'color',[0 0 0.3]);
line(ERF_cue_high_filt_sel_pos_GA_mag.time, ERF_cue_high_filt_sel_pos_GA_mag.avg,'linewidth',1,'color',[0 0 0.6]);
axis tight

% stim times
ax = axis;
stim_times = [0:1/16:1.5];
y1 = ones(size(stim_times))*ax(3);
y2 = ones(size(stim_times))*ax(4);
for i = 1 : length(stim_times)
    line([stim_times(i),stim_times(i)],[y1(i), y1(i) + abs(y1(i))*0.1],'color',[0 0 0]);
end
set(gca,'PlotBoxAspectRatio',[8/3,1,1]);
savefig(fig,'w:\WANDER\images\poster\FIGb_onset_timecourse');

% saveas(fig,'w:\WANDER\images\FIGb_onset_timecourse.eps','eps');

%% statistics FIRST COMPONENT

neigh_cmb = load('neuromag306cmb_neighb');
neigh_mag = load('neuromag306mag_neighb');

% statistics
Nsub                    = length(slist);
cfg                     = [];
cfg.latency             = [time_max_mag-0.01 time_max_mag+0.01];
cfg.avgovertime         = 'yes';
cfg.avgoverchan         = 'yes';
cfg.method              = 'analytic';
cfg.statistic           = 'ft_statfun_depsamplesT';
% cfg.correctm            = 'cluster';
% cfg.clusteralpha        = 0.05;
% cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.05;
cfg.numrandomization    = 2000;
cfg.neighbours          = neigh_mag.neighbours; %  ft_prepare_neighbours(cfg_neighb, grad);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

cfg.channel             = chan_pos_mag;
stat                    = ft_timelockstatistics(cfg,ERF_cue_low_filt{slist},ERF_cue_high_filt{slist});

cfg.channel             = chan_neg_mag;
stat                    = ft_timelockstatistics(cfg,ERF_cue_low_filt{slist},ERF_cue_high_filt{slist});

cfg.channel             = chan_pos_cmb;
stat                    = ft_timelockstatistics(cfg,ERF_cue_low_filt{slist},ERF_cue_high_filt{slist});

%% get highest pos. and neg. deflections at SECOND COMPONENT
peak_latency            = [0.100 0.200];

cfg                     = [];
cfg.latency             = peak_latency;
cfg.avgovertime         = 'yes';
cfg.avgovertime         = 'no';

cfg.channel             = 'MEG*3';
temp                    = ft_selectdata(cfg,ERF_cue_filt_GA);

% max timepoint
[Y_time,I_time]         = max(mean(abs(temp.avg)));
time_max_cmb            = temp.time(I_time);

% max sensor at max timepoint
[Y_sens,I_sens]         = max(temp.avg(:,I_time));
chan_pos_cmb            = ft_channelselection(I_sens, temp.label);

cfg.channel             = 'MEG*1';
temp                    = ft_selectdata(cfg,ERF_cue_filt_GA);

% max/min timepoint
[Y_time,I_time]         = max(mean(abs(temp.avg)));
time_max_mag            = temp.time(I_time);

% max sensor at max timepoint
[Y_sens,I_sens]         = max(temp.avg(:,I_time));
chan_pos_mag            = ft_channelselection(I_sens, temp.label);

% min sensor at max timepoint
[Y_sens,I_sens]         = min(temp.avg(:,I_time));
chan_neg_mag            = ft_channelselection(I_sens, temp.label);

% select data only of max min sensors
cfg = [];
cfg.parameter                       = {'avg','std'};
cfg.avgoverchan                     = 'yes';
cfg.channel                         = chan_neg_mag;
ERF_cue_low_filt_sel_neg_GA_mag     = ft_selectdata(cfg,ERF_cue_low_filt_GA);
ERF_cue_high_filt_sel_neg_GA_mag    = ft_selectdata(cfg,ERF_cue_high_filt_GA);
cfg.channel                         = chan_pos_mag;
ERF_cue_low_filt_sel_pos_GA_mag     = ft_selectdata(cfg,ERF_cue_low_filt_GA);
ERF_cue_high_filt_sel_pos_GA_mag    = ft_selectdata(cfg,ERF_cue_high_filt_GA);
cfg.channel                         = chan_pos_cmb;
ERF_cue_low_filt_sel_pos_GA_cmb     = ft_selectdata(cfg,ERF_cue_low_filt_GA);
ERF_cue_high_filt_sel_pos_GA_cmb    = ft_selectdata(cfg,ERF_cue_high_filt_GA);

%% statistics SECOND COMPONENT

neigh_cmb = load('neuromag306cmb_neighb');
neigh_mag = load('neuromag306mag_neighb');

% statistics
Nsub                    = length(slist);
cfg                     = [];
cfg.latency             = [time_max_mag time_max_mag];
cfg.avgovertime         = 'yes';
cfg.avgoverchan         = 'yes';
cfg.method              = 'analytic';
cfg.statistic           = 'ft_statfun_depsamplesT';
% cfg.correctm            = 'cluster';
% cfg.clusteralpha        = 0.05;
% cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.05;
cfg.numrandomization    = 2000;
cfg.neighbours          = neigh_mag.neighbours; %  ft_prepare_neighbours(cfg_neighb, grad);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
% 
% cfg.channel             = chan_pos_mag;
% stat                    = ft_timelockstatistics(cfg,ERF_cue_low_filt{slist},ERF_cue_high_filt{slist});

cfg.channel             = chan_neg_mag;
stat                    = ft_timelockstatistics(cfg,ERF_cue_low_filt{slist},ERF_cue_high_filt{slist});

% cfg.channel             = chan_pos_cmb;
% stat                    = ft_timelockstatistics(cfg,ERF_cue_low_filt{slist},ERF_cue_high_filt{slist});

%% Calculate TFR at SSSF frequency
for isubject = slist
    cfg                         = [];
    cfg.pad                     = 'nextpow2';
    cfg.channel                 = 'all';
    cfg.method                  = 'mtmconvol';
    cfg.foi                     = [1:40];
    cfg.taper                   = 'hanning';
    cfg.t_ftimwin               = ones(size(cfg.foi));
    cfg.toi                     = -1:0.05:10;
    
    TFR{isubject}           = ft_freqanalysis(cfg, ERF_cue{isubject});
    TFR_high{isubject}      = ft_freqanalysis(cfg, ERF_cue_high{isubject});
    TFR_low{isubject}       = ft_freqanalysis(cfg, ERF_cue_low{isubject});
end

% combine planar
for isubject = slist
    TFR{isubject}      = ft_combineplanar([],TFR{isubject});
    TFR_high{isubject} = ft_combineplanar([],TFR_high{isubject});
    TFR_low{isubject}  = ft_combineplanar([],TFR_low{isubject});
end

% difference per subject
for isubject = slist
    TFR_diff{isubject}                  =  TFR_low{isubject};
    TFR_diff{isubject}.powspctrm        =  (TFR_low{isubject}.powspctrm - TFR_high{isubject}.powspctrm)  ./ (TFR_low{isubject}.powspctrm + TFR_high{isubject}.powspctrm);
    TFR_low_norm{isubject}              =  TFR_low{isubject};
    TFR_low_norm{isubject}.powspctrm    = (TFR_low{isubject}.powspctrm) ./ (TFR_low{isubject}.powspctrm + TFR_high{isubject}.powspctrm);    
    TFR_high_norm{isubject}             =  TFR_low{isubject};
    TFR_high_norm{isubject}.powspctrm   = (TFR_high{isubject}.powspctrm) ./ (TFR_low{isubject}.powspctrm + TFR_high{isubject}.powspctrm);        
end

% average over subjects
TFR_GA                  = ft_freqgrandaverage([],TFR{slist});
TFR_high_GA             = ft_freqgrandaverage([],TFR_high{slist});
TFR_low_GA              = ft_freqgrandaverage([],TFR_low{slist});
TFR_diff_GA             = ft_freqgrandaverage([],TFR_diff{slist});

% get maximum sensor
peak_latency_TFR        = [1 9.5];

cfg                     = [];
cfg.latency             = peak_latency_TFR;
cfg.avgovertime         = 'yes';
cfg.channel             = 'MEG*1';
cfg.frequency           = [16 16];
cfg.avgoverfreq         = 'yes';

temp                    = ft_selectdata(cfg,TFR_GA);
[Y I]                   = sort(temp.powspctrm,'ascend');
chan_TFR_mag            = ft_channelselection(I(end), temp.label);

cfg.channel             = 'MEG*3';
temp                    = ft_selectdata(cfg,TFR_GA);
[Y I]                   = sort(temp.powspctrm,'ascend');
chan_TFR_cmb            = ft_channelselection(I(end), temp.label);

% grand averages for plotting 16Hz
i = 1; clear TFR_ind*
for isubject = slist
    cfg                     = [];
    cfg.frequency           = [16 16];
    cfg.avgoverfreq         = 'yes';
    cfg.channel             = chan_TFR_mag;
    cfg.latency             = [-1 9.5];

    temp                    = ft_selectdata(cfg,TFR{isubject});
    TFR_ind_mag(i,:)        = squeeze(temp.powspctrm);
    temp                    = ft_selectdata(cfg,TFR_high{isubject});
    TFR_ind_high_mag(i,:)   = squeeze(temp.powspctrm);
    temp                    = ft_selectdata(cfg,TFR_low{isubject});
    TFR_ind_low_mag(i,:)    = squeeze(temp.powspctrm);
    temp                    = ft_selectdata(cfg,TFR_diff{isubject});
    TFR_ind_diff_mag(i,:)   = squeeze(temp.powspctrm);
    
    cfg.channel             = chan_TFR_cmb;
    temp                    = ft_selectdata(cfg,TFR{isubject});
    TFR_ind_cmb(i,:)        = squeeze(temp.powspctrm);
    temp                    = ft_selectdata(cfg,TFR_high{isubject});
    TFR_ind_high_cmb(i,:)   = squeeze(temp.powspctrm);
    temp                    = ft_selectdata(cfg,TFR_low{isubject});
    TFR_ind_low_cmb(i,:)    = squeeze(temp.powspctrm);    
    temp                    = ft_selectdata(cfg,TFR_diff{isubject});
    TFR_ind_diff_cmb(i,:)   = squeeze(temp.powspctrm);    
    i = i + 1;
end

cfg = [];
cfg.foilim                      = [16 16];
cfg.channel                     = chan_TFR_cmb;
TFR_16hz_cmb_GA                 = ft_freqgrandaverage(cfg,TFR{slist});
TFR_16hz_high_cmb_GA            = ft_freqgrandaverage(cfg,TFR_high{slist});
TFR_16hz_low_cmb_GA             = ft_freqgrandaverage(cfg,TFR_low{slist});
TFR_16hz_diff_cmb_GA            = ft_freqgrandaverage(cfg,TFR_diff{slist});

cfg.channel                     = chan_TFR_mag;
TFR_16hz_mag_GA                 = ft_freqgrandaverage(cfg,TFR{slist});
TFR_16hz_high_mag_GA            = ft_freqgrandaverage(cfg,TFR_high{slist});
TFR_16hz_low_mag_GA             = ft_freqgrandaverage(cfg,TFR_low{slist});
TFR_16hz_diff_mag_GA            = ft_freqgrandaverage(cfg,TFR_diff{slist});

% the rest is just to remove the first timepoints which are nan
TFR_16hz_cmb_GA.time            = temp.time;
TFR_16hz_high_cmb_GA.time       = temp.time;
TFR_16hz_low_cmb_GA.time        = temp.time;
TFR_16hz_mag_GA.time            = temp.time;
TFR_16hz_high_mag_GA.time       = temp.time;
TFR_16hz_low_mag_GA.time        = temp.time;

TFR_16hz_mag_GA.powspctrm       = nanmean(TFR_ind_mag);
TFR_16hz_high_mag_GA.powspctrm  = nanmean(TFR_ind_high_mag);
TFR_16hz_low_mag_GA.powspctrm   = nanmean(TFR_ind_low_mag);
TFR_16hz_diff_mag_GA.powspctrm  = nanmean(TFR_ind_diff_mag);

TFR_16hz_cmb_GA.powspctrm       = nanmean(TFR_ind_cmb);
TFR_16hz_high_cmb_GA.powspctrm  = nanmean(TFR_ind_high_cmb);
TFR_16hz_low_cmb_GA.powspctrm   = nanmean(TFR_ind_low_cmb);
TFR_16hz_diff_cmb_GA.powspctrm  = nanmean(TFR_ind_diff_cmb);

TFR_16hz_cmb_GA.std             = nanstd(TFR_ind_cmb,1);
TFR_16hz_high_cmb_GA.std        = nanstd(TFR_ind_high_cmb,1);
TFR_16hz_low_cmb_GA.std         = nanstd(TFR_ind_low_cmb,1);
TFR_16hz_diff_cmb_GA.std        = nanstd(TFR_ind_diff_cmb,1);

TFR_16hz_mag_GA.std             = nanstd(TFR_ind_mag,1);
TFR_16hz_high_mag_GA.std        = nanstd(TFR_ind_high_mag,1);
TFR_16hz_low_mag_GA.std         = nanstd(TFR_ind_low_mag,1);
TFR_16hz_diff_mag_GA.std        = nanstd(TFR_ind_diff_mag,1);

%% statistics TFR

neigh_cmb = load('neuromag306cmb_neighb');
neigh_mag = load('neuromag306mag_neighb');

% statistics
Nsub                    = length(slist);
cfg                     = [];
cfg.latency             = [1 9.5];
cfg.avgovertime         = 'yes';
cfg.avgoverchan         = 'yes';
cfg.frequency           = [16 16];
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 1;
cfg.tail                = 1;
cfg.clustertail         = 1;
cfg.minnbchan           = 1;
cfg.alpha               = 0.05;
cfg.numrandomization    = 2000;
cfg.neighbours          = neigh_mag.neighbours; %  ft_prepare_neighbours(cfg_neighb, grad);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

cfg.channel             = chan_TFR_mag;
stat_TFR_mag            = ft_freqstatistics(cfg,TFR_low_norm{slist},TFR_high_norm{slist});

cfg.channel             = chan_TFR_cmb;
stat_TFR_cmb            = ft_freqstatistics(cfg,TFR_low_norm{slist},TFR_high_norm{slist});

%% Plot overview SSEF TFR results
fig = figure;

subplot(4,2,1);
cfg                     = [];
cfg.layout              = 'neuromag306mag';
cfg.xlim                = peak_latency_TFR;
cfg.zlim                = 'maxabs';
cfg.highlightchannel    = [chan_TFR_mag];
cfg.highlight           = 'yes';
cfg.highlightsymbol     = 'pentagram';
cfg.highlightcolor      = [0 0 0];
cfg.highlightsize       = 16;
cfg.marker              = 'off';
cfg.comment             = 'no';
cfg.colorbar            = 'yes';
cfg.ylim                = [16 16];
ft_topoplotTFR(cfg,TFR_GA);

subplot(4,2,2);
cfg                     = [];
cfg.layout              = 'neuromag306cmb';
cfg.xlim                = peak_latency_TFR;
cfg.zlim                = 'maxabs';
cfg.highlightchannel    = chan_TFR_cmb;
cfg.highlight           = 'yes';
cfg.highlightsymbol     = '.';
cfg.highlightcolor      = [0 0 0];
cfg.highlightsize       = 20;
cfg.marker              = 'off';
cfg.comment             = 'no';
cfg.colorbar            = 'yes';
cfg.ylim                = [16 16];
ft_topoplotTFR(cfg,TFR_GA);

subplot(4,2,3);
cfg                     = [];
cfg.layout              = 'neuromag306cmb';
cfg.xlim                = [-1 10];
cfg.zlim                = 'maxabs';
cfg.channel             = chan_TFR_mag;
cfg.baseline            = [-1 -0.5];

cfg.baselinetype        = 'relative';
cfg.interactive         = 'no';
cfg.ylim                = [1 30];
ft_singleplotTFR(cfg,TFR_GA);

subplot(4,2,4);
cfg                     = [];
cfg.layout              = 'neuromag306cmb';
cfg.xlim                = [-1 10];
cfg.zlim                = 'maxabs';
cfg.channel             = chan_TFR_cmb;
cfg.baseline            = [-1 -0.5];
cfg.baselinetype        = 'relative';
cfg.interactive         = 'no';
cfg.ylim                = [1 30];
ft_singleplotTFR(cfg,TFR_GA);

subplot(4,2,5);
stim_times = [0:1/16:10];
y1 = ones(size(stim_times))*0;
y2 = ones(size(stim_times))*1.3e-27;
for i = 1 : length(stim_times)
    line([stim_times(i),stim_times(i)],[y1', y2'.*0.05],'color',[0 0 0]);
end
ax = axis;
patch([peak_latency_TFR, peak_latency_TFR(2:-1:1)],[0,0 1.3e-27,1.3e-27],[0 0 0 ],'FaceAlpha',.2,'LineStyle','None');
patch([TFR_16hz_high_mag_GA.time, TFR_16hz_high_mag_GA.time(end:-1:1)], [TFR_16hz_high_mag_GA.powspctrm+TFR_16hz_high_mag_GA.std ./sqrt(22)*t,    TFR_16hz_high_mag_GA.powspctrm(end:-1:1)- TFR_16hz_high_mag_GA.std(end:-1:1) ./sqrt(22)*t],[1 0 0 ],'FaceAlpha',.5,'LineStyle','None');
patch([TFR_16hz_low_mag_GA.time,  TFR_16hz_low_mag_GA.time(end:-1:1)],  [TFR_16hz_low_mag_GA.powspctrm+TFR_16hz_low_mag_GA.std ./sqrt(22)*t,      TFR_16hz_low_mag_GA.powspctrm(end:-1:1)-  TFR_16hz_low_mag_GA.std(end:-1:1)  ./sqrt(22)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
axis tight

subplot(4,2,6);
stim_times = [0:1/16:10];
y1 = ones(size(stim_times))*0;
y2 = ones(size(stim_times))*1.7e-24;
for i = 1 : length(stim_times)
    line([stim_times(i),stim_times(i)],[y1', y2'.*0.05],'color',[0 0 0]);
end
ax = axis;
patch([peak_latency_TFR, peak_latency_TFR(2:-1:1)],[0,0 1.7e-24,1.7e-24],[0 0 0 ],'FaceAlpha',.2,'LineStyle','None');
patch([TFR_16hz_high_cmb_GA.time, TFR_16hz_high_cmb_GA.time(end:-1:1)], [TFR_16hz_high_cmb_GA.powspctrm+TFR_16hz_high_cmb_GA.std ./sqrt(22)*t,    TFR_16hz_high_cmb_GA.powspctrm(end:-1:1)- TFR_16hz_high_cmb_GA.std(end:-1:1) ./sqrt(22)*t],[1 0 0 ],'FaceAlpha',.5,'LineStyle','None');
patch([TFR_16hz_low_cmb_GA.time,  TFR_16hz_low_cmb_GA.time(end:-1:1)],  [TFR_16hz_low_cmb_GA.powspctrm+TFR_16hz_low_cmb_GA.std ./sqrt(22)*t,      TFR_16hz_low_cmb_GA.powspctrm(end:-1:1)-  TFR_16hz_low_cmb_GA.std(end:-1:1)  ./sqrt(22)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
axis tight

subplot(4,2,7); hold;
line(TFR_16hz_high_mag_GA.time,zeros(size(TFR_16hz_diff_mag_GA.powspctrm,2)),'color','k','linestyle','--');
patch([TFR_16hz_high_mag_GA.time, TFR_16hz_high_mag_GA.time(end:-1:1)], [TFR_16hz_diff_mag_GA.powspctrm+TFR_16hz_diff_mag_GA.std ./sqrt(22)*t, TFR_16hz_diff_mag_GA.powspctrm(end:-1:1)- TFR_16hz_diff_mag_GA.std(end:-1:1) ./sqrt(22)*t],[0 0 0 ],'FaceAlpha',.5,'LineStyle','None');
plot(TFR_16hz_high_mag_GA.time,TFR_16hz_diff_mag_GA.powspctrm,'k');
scatter(TFR_16hz_high_mag_GA.time(stat_TFR_mag.mask),TFR_16hz_diff_mag_GA.powspctrm(stat_TFR_mag.mask),8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
axis tight

subplot(4,2,8); hold;
line(TFR_16hz_high_mag_GA.time,zeros(size(TFR_16hz_diff_mag_GA.powspctrm,2)),'color','k','linestyle','--');
patch([TFR_16hz_high_mag_GA.time, TFR_16hz_high_mag_GA.time(end:-1:1)], [TFR_16hz_diff_cmb_GA.powspctrm+TFR_16hz_diff_cmb_GA.std ./sqrt(22)*t, TFR_16hz_diff_cmb_GA.powspctrm(end:-1:1)- TFR_16hz_diff_cmb_GA.std(end:-1:1) ./sqrt(22)*t],[0 0 0 ],'FaceAlpha',.5,'LineStyle','None');
plot(TFR_16hz_high_mag_GA.time,TFR_16hz_diff_cmb_GA.powspctrm,'k');
scatter(TFR_16hz_high_cmb_GA.time(stat_TFR_cmb.mask),TFR_16hz_diff_cmb_GA.powspctrm(stat_TFR_cmb.mask),8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
axis tight

set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0 11*4 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['d:\analysis\WANDER\images\SSTFR_overview_cue.pdf']);
print(fig,'-dpng',['d:\analysis\WANDER\images\SSTFR_overview_cue.png']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NOW FOR PROBELOCKED

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD
timing = 'cue';
for isubject = slist
    fprintf('Loading subject %d \n',isubject);    
    temp = load(['d:\analysis\WANDER\data\ERF\s' num2str(isubject) '_' timing '.mat'],'ERF_probe','ERF_probe_high','ERF_probe_low');
    ERF{isubject}       = temp.ERF_probe;
    ERF_high{isubject}  = temp.ERF_probe_high;
    ERF_low{isubject}   = temp.ERF_probe_low;
    
    fprintf('Cutting subject %d to size \n',isubject);    
    cfg = [];
    cfg.latency = [-1.5 4]; % wider for TFR - CUE
    cfg.latency = [-10 1]; % wider for TFR - PROBE

    cfg.channel = 'MEG*';
    ERF{isubject}       = ft_selectdata(cfg,ERF{isubject});   
    ERF_low{isubject}   = ft_selectdata(cfg,ERF_low{isubject});
    ERF_high{isubject}  = ft_selectdata(cfg,ERF_high{isubject});
end

% for filtered ERF, so not matter we already combine
clear temp
for isubject = slist
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 30; 
    ERF_filt{isubject}          = ft_preprocessing(cfg,ERF{isubject});
    ERF_low_filt{isubject}      = ft_preprocessing(cfg,ERF_low{isubject});
    ERF_high_filt{isubject}     = ft_preprocessing(cfg,ERF_high{isubject});  
    ERF_diff_filt{isubject}     = ERF_low_filt{isubject};  
    ERF_diff_filt{isubject}.avg = (ERF_low_filt{isubject}.avg - ERF_low_filt{isubject}.avg) ./ (ERF_low_filt{isubject}.avg + ERF_low_filt{isubject}.avg);  

    cfg = [];   
    ERF_filt{isubject}       = ft_combineplanar(cfg,ERF_filt{isubject});
    ERF_low_filt{isubject}   = ft_combineplanar(cfg,ERF_low_filt{isubject});
    ERF_high_filt{isubject}  = ft_combineplanar(cfg,ERF_high_filt{isubject}); 
    ERF_diff_filt{isubject}  = ft_combineplanar(cfg,ERF_diff_filt{isubject});    
end

% grand averages for plotting power + std
i = 1; clear avg_*
for isubject = slist
    cfg                     = [];
    cfg.latency             = [-0.5 3]; % for plotting CUELOCKED
    cfg.latency             = [-3 1]; % for plotting PROBELOCKED
    cfg.latency             = [-0.5 1]; % for plotting PROBELOCKED - specific for offset response    
    temp                    = ft_selectdata(cfg,ERF_filt{isubject});
    avg_filt(i,:,:)         = temp.avg;
    temp                    = ft_selectdata(cfg,ERF_high_filt{isubject});
    avg_high_filt(i,:,:)    = temp.avg;
    temp                    = ft_selectdata(cfg,ERF_low_filt{isubject});
    avg_low_filt(i,:,:)     = temp.avg;
    temp                    = ft_selectdata(cfg,ERF_diff_filt{isubject});
    avg_diff_filt(i,:,:)    = temp.avg;    
    i = i + 1;
end

cfg = [];
cfg.latency                 = [-0.5 3];% for plotting CUELOCKED
cfg.latency                 = [-3 1];% for plotting PROBELOCKED
cfg.latency                 = [-0.5 1]; % for plotting PROBELOCKED - specific for offset response    

ERF_filt_GA                 = ft_timelockgrandaverage(cfg,ERF_filt{slist});
ERF_high_filt_GA            = ft_timelockgrandaverage(cfg,ERF_high_filt{slist});
ERF_low_filt_GA             = ft_timelockgrandaverage(cfg,ERF_low_filt{slist});
ERF_diff_filt_GA            = ft_timelockgrandaverage(cfg,ERF_diff_filt{slist});
ERF_filt_GA.std             = squeeze(std(avg_filt,1));
ERF_high_filt_GA.std        = squeeze(std(avg_high_filt,1));
ERF_low_filt_GA.std         = squeeze(std(avg_low_filt,1));
ERF_diff_filt_GA.std        = squeeze(std(avg_diff_filt,1));

%% get highest pos. and neg. deflections at offest response FIRST COMPONENT

peak_latency            = [0.200 0.270]; % for PROBELOCKED

cfg                     = [];
cfg.latency             = peak_latency;
cfg.avgovertime         = 'yes';
cfg.avgovertime         = 'no';

cfg.channel             = 'MEG*3';
temp                    = ft_selectdata(cfg,ERF_filt_GA);

% max timepoint
[Y_time,I_time]         = max(mean(abs(temp.avg)));
time_max_cmb            = temp.time(I_time);

% max sensor at max timepoint
[Y_sens,I_sens]         = max(temp.avg(:,I_time));
chan_pos_cmb            = ft_channelselection(I_sens, temp.label);

cfg.channel             = 'MEG*1';
temp                    = ft_selectdata(cfg,ERF_filt_GA);

% max/min timepoint
[Y_time,I_time]         = max(mean(abs(temp.avg)));
time_max_mag            = temp.time(I_time);

% max sensor at max timepoint
[Y_sens,I_sens]         = max(temp.avg(:,I_time));
chan_pos_mag            = ft_channelselection(I_sens, temp.label);

% min sensor at max timepoint
[Y_sens,I_sens]         = min(temp.avg(:,I_time));
chan_neg_mag            = ft_channelselection(I_sens, temp.label);

% select data only of max min sensors
cfg = [];
cfg.parameter                   = {'avg','std'};
cfg.avgoverchan                 = 'yes';
cfg.channel                     = chan_neg_mag;
ERF_low_filt_sel_neg_GA_mag     = ft_selectdata(cfg,ERF_low_filt_GA);
ERF_high_filt_sel_neg_GA_mag    = ft_selectdata(cfg,ERF_high_filt_GA);
cfg.channel                     = chan_pos_mag;
ERF_low_filt_sel_pos_GA_mag     = ft_selectdata(cfg,ERF_low_filt_GA);
ERF_high_filt_sel_pos_GA_mag    = ft_selectdata(cfg,ERF_high_filt_GA);
cfg.channel                     = chan_pos_cmb;
ERF_low_filt_sel_pos_GA_cmb     = ft_selectdata(cfg,ERF_low_filt_GA);
ERF_high_filt_sel_pos_GA_cmb    = ft_selectdata(cfg,ERF_high_filt_GA);

%% statistics OFFSET response, FIRST COMPONENT 

neigh_cmb = load('neuromag306cmb_neighb');
neigh_mag = load('neuromag306mag_neighb');
pad = 0.01;
% statistics
Nsub                    = length(slist);
cfg                     = [];
cfg.latency             = [time_max_mag-pad time_max_mag+pad];
cfg.avgovertime         = 'yes';
cfg.avgoverchan         = 'yes';
cfg.method              = 'montecarlo';
cfg.method              = 'analytic';

cfg.statistic           = 'ft_statfun_depsamplesT';
% cfg.correctm            = 'cluster';
% cfg.clusteralpha        = 0.05;
% cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.05;
cfg.numrandomization    = 2000;
cfg.neighbours          = neigh_mag.neighbours; %  ft_prepare_neighbours(cfg_neighb, grad);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

cfg.channel             = chan_pos_mag;
stat_pos_mag_1          = ft_timelockstatistics(cfg,ERF_low_filt{slist},ERF_high_filt{slist});

cfg.channel             = chan_neg_mag;
stat_neg_mag_1          = ft_timelockstatistics(cfg,ERF_low_filt{slist},ERF_high_filt{slist});

cfg.channel             = chan_pos_cmb;
stat_cmb_1              = ft_timelockstatistics(cfg,ERF_low_filt{slist},ERF_high_filt{slist});

%% plot overview with ERF pos, neg and topo OFFSET response FIRST COMPONENT

% for article
fig = figure;
cfg                     = [];
cfg.layout              = 'neuromag306mag';
cfg.xlim                = [time_max_mag-pad time_max_mag+pad];
cfg.zlim                = 'maxabs';
% cfg.highlightchannel    = chan_neg_mag;
% cfg.highlight           = 'yes';
% cfg.highlightsymbol     = '.';
% cfg.highlightcolor      = [0 0 0];
% cfg.highlightsize       = 20;
cfg.marker              = 'off';
cfg.comment             = 'no';
cfg.colorbar            = 'yes';
cfg.colormap            = colormap(jet(1024));
cfg.gridscale           = 64;
ft_topoplotER(cfg,ERF_filt_GA);

% set(gca,'PlotBoxAspectRatio',[8/3,1,1]);
saveas(fig,'w:\WANDER\images\FIGd_offset_topo_1st.eps','eps');

fig = figure;
t = 1;
stim_times = -0.5:1/16:0;

axis([-0.5,1,-10e-14,20e-14]); ax = axis; hold;
patch([time_max_mag-pad time_max_mag+pad time_max_mag+pad time_max_mag-pad],[ax(3),ax(3) ax(4),ax(4)],[0 0 0 ],'FaceAlpha',.2,'LineStyle','None');
patch([ERF_low_filt_sel_pos_GA_mag.time, ERF_low_filt_sel_pos_GA_mag.time(end:-1:1)],[ERF_low_filt_sel_pos_GA_mag.avg+ERF_low_filt_sel_pos_GA_mag.std ./sqrt(22)*t, ERF_low_filt_sel_pos_GA_mag.avg(end:-1:1)-ERF_low_filt_sel_pos_GA_mag.std(end:-1:1) ./sqrt(21)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
patch([ERF_high_filt_sel_pos_GA_mag.time, ERF_high_filt_sel_pos_GA_mag.time(end:-1:1)],[ERF_high_filt_sel_pos_GA_mag.avg+ERF_high_filt_sel_pos_GA_mag.std ./sqrt(22)*t, ERF_high_filt_sel_pos_GA_mag.avg(end:-1:1)-ERF_high_filt_sel_pos_GA_mag.std(end:-1:1) ./sqrt(21)*t],[1 0 0 ],'FaceAlpha',.5,'LineStyle','None');
y1 = ones(size(stim_times))*ax(3);
y2 = ones(size(stim_times))*ax(4);
for i = 1 : length(stim_times)
    line([stim_times(i),stim_times(i)],[y1(i), y1(i) + abs(y1(i))*0.1],'color',[0 0 0]);
end
title(sprintf('Positive channel, t=%f, p=%f\n',stat_pos_mag_1.stat,stat_pos_mag_1.prob));

set(gca,'PlotBoxAspectRatio',[8/3,1,1]);
savefig(fig,'w:\WANDER\images\poster\offset_timecourse_1st');




% 
% saveas(fig,'w:\WANDER\images\FIG1c_offset_time.eps','eps2c');
% saveas(fig,'w:\WANDER\images\FIG1c_offset_time.png','png');

subplot(3,2,5);
axis([-0.5,1,-20e-14,15e-14]); ax = axis; hold;
patch([time_max_mag-pad time_max_mag+pad time_max_mag+pad time_max_mag-pad],[ax(3),ax(3) ax(4),ax(4)],[0 0 0 ],'FaceAlpha',.2,'LineStyle','None');
patch([ERF_low_filt_sel_neg_GA_mag.time, ERF_low_filt_sel_neg_GA_mag.time(end:-1:1)],[ERF_low_filt_sel_neg_GA_mag.avg+ERF_low_filt_sel_neg_GA_mag.std ./sqrt(22)*t, ERF_low_filt_sel_neg_GA_mag.avg(end:-1:1)-ERF_low_filt_sel_neg_GA_mag.std(end:-1:1) ./sqrt(21)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
patch([ERF_high_filt_sel_neg_GA_mag.time, ERF_high_filt_sel_neg_GA_mag.time(end:-1:1)],[ERF_high_filt_sel_neg_GA_mag.avg+ERF_high_filt_sel_neg_GA_mag.std ./sqrt(22)*t, ERF_high_filt_sel_neg_GA_mag.avg(end:-1:1)-ERF_high_filt_sel_neg_GA_mag.std(end:-1:1) ./sqrt(21)*t],[1 0 0 ],'FaceAlpha',.5,'LineStyle','None');
y1 = ones(size(stim_times))*ax(3);
y2 = ones(size(stim_times))*ax(4);
for i = 1 : length(stim_times)
    line([stim_times(i),stim_times(i)],[y1(i), y1(i) + abs(y1(i))*0.1],'color',[0 0 0]);
end
title(sprintf('Negative channel, t=%f, p=%f\n',stat_neg_mag_1.stat,stat_neg_mag_1.prob));

subplot(3,2,4);
axis([-0.5,1,0,8e-12]); ax = axis; hold;
patch([time_max_cmb-pad time_max_cmb+pad time_max_cmb+pad time_max_cmb-pad],[ax(3),ax(3) ax(4),ax(4)],[0 0 0 ],'FaceAlpha',.2,'LineStyle','None');
patch([ERF_low_filt_sel_pos_GA_cmb.time,  ERF_low_filt_sel_pos_GA_cmb.time(end:-1:1)], [ERF_low_filt_sel_pos_GA_cmb.avg  + ERF_low_filt_sel_pos_GA_cmb.std  ./sqrt(22)*t, ERF_low_filt_sel_pos_GA_cmb.avg(end:-1:1)-ERF_low_filt_sel_pos_GA_cmb.std(end:-1:1) ./sqrt(21)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
patch([ERF_high_filt_sel_pos_GA_cmb.time, ERF_high_filt_sel_pos_GA_cmb.time(end:-1:1)],[ERF_high_filt_sel_pos_GA_cmb.avg + ERF_high_filt_sel_pos_GA_cmb.std ./sqrt(22)*t, ERF_high_filt_sel_pos_GA_cmb.avg(end:-1:1)-ERF_high_filt_sel_pos_GA_cmb.std(end:-1:1) ./sqrt(21)*t],[1 0 0 ],'FaceAlpha',.5,'LineStyle','None');
y1 = ones(size(stim_times))*ax(3);
y2 = ones(size(stim_times))*ax(4);
for i = 1 : length(stim_times)
    line([stim_times(i),stim_times(i)],[0, abs(y2(i))*0.05],'color',[0 0 0]);
end
title(sprintf('Maximum channel, t=%f, p=%f\n',stat_cmb_1.stat,stat_cmb_1.prob));

set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0 11*4 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['d:\analysis\WANDER\images\ERF_overview_probe.pdf']);
print(fig,'-dpng',['d:\analysis\WANDER\images\ERF_overview_probe.png']);

set(gca,'PlotBoxAspectRatio',[8/3,1,1]);
% saveas(fig,'w:\WANDER\images\FIGb_onset_timecourse.eps','eps');


%% get highest pos. and neg. deflections at SECOND COMPONENT
peak_latency            = [0.280 0.320]; % for PROBELOCKED

peak_latency            = [0.270 0.400]; % for PROBELOCKED

cfg                     = [];
cfg.latency             = peak_latency;
cfg.avgovertime         = 'yes';
cfg.avgovertime         = 'no';

cfg.channel             = 'MEG*3';
% 
% % only right side
% cfg.channel = {'MEG0722+0723', 'MEG0732+0733', 'MEG0912+0913', 'MEG0922+0923', 'MEG0932+0933', 'MEG0942+0943', 'MEG1022+1023', 'MEG1032+1033', 'MEG1042+1043', 'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG1212+1213', 'MEG1222+1223', 'MEG1232+1233', 'MEG1242+1243', 'MEG1312+1313', 'MEG1322+1323', 'MEG1332+1333', 'MEG1342+1343', 'MEG1412+1413', 'MEG1422+1423', 'MEG1432+1433', 'MEG1442+1443', 'MEG2022+2023', 'MEG2032+2033', 'MEG2132+2133', 'MEG2212+2213', 'MEG2222+2223', 'MEG2232+2233', 'MEG2242+2243', 'MEG2312+2313', 'MEG2322+2323', 'MEG2332+2333', 'MEG2342+2343', 'MEG2412+2413', 'MEG2422+2423', 'MEG2432+2433', 'MEG2442+2443', 'MEG2512+2513', 'MEG2522+2523', 'MEG2532+2533', 'MEG2542+2543', 'MEG2612+2613', 'MEG2622+2623', 'MEG2632+2633', 'MEG2642+2643'};

temp                    = ft_selectdata(cfg,ERF_filt_GA);

% max timepoint
[Y_time,I_time]         = max(mean(abs(temp.avg)));
time_max_cmb            = temp.time(I_time);

% max sensor at max timepoint
[Y_sens,I_sens]         = max(temp.avg(:,I_time));
chan_pos_cmb            = ft_channelselection(I_sens, temp.label);

cfg.channel             = 'MEG*1';
% 
% % only right side
% cfg.channel             = {'MEG0721', 'MEG0731', 'MEG0911', 'MEG0921', 'MEG0931', 'MEG0941', 'MEG1021', 'MEG1031', 'MEG1041', 'MEG1111', 'MEG1121', 'MEG1131', 'MEG1141', 'MEG1211', 'MEG1221', 'MEG1231', 'MEG1241', 'MEG1311', 'MEG1321', 'MEG1331', 'MEG1341', 'MEG1411', 'MEG1421', 'MEG1431', 'MEG1441', 'MEG2021', 'MEG2031', 'MEG2131', 'MEG2211', 'MEG2221', 'MEG2231', 'MEG2241', 'MEG2311', 'MEG2321', 'MEG2331', 'MEG2341', 'MEG2411', 'MEG2421', 'MEG2431', 'MEG2441', 'MEG2511', 'MEG2521', 'MEG2531', 'MEG2541', 'MEG2611', 'MEG2621', 'MEG2631', 'MEG2641'};
% 
temp                    = ft_selectdata(cfg,ERF_filt_GA);

% max/min timepoint
[Y_time,I_time]         = max(mean(abs(temp.avg)));
time_max_mag            = temp.time(I_time);

% max sensor at max timepoint
[Y_sens,I_sens]         = max(temp.avg(:,I_time));
chan_pos_mag            = ft_channelselection(I_sens, temp.label);

% min sensor at max timepoint
[Y_sens,I_sens]         = min(temp.avg(:,I_time));
chan_neg_mag            = ft_channelselection(I_sens, temp.label);

% select data only of max min sensors
cfg = [];
cfg.parameter                   = {'avg','std'};
cfg.avgoverchan                 = 'yes';
cfg.channel                     = chan_neg_mag;
ERF_low_filt_sel_neg_GA_mag     = ft_selectdata(cfg,ERF_low_filt_GA);
ERF_high_filt_sel_neg_GA_mag    = ft_selectdata(cfg,ERF_high_filt_GA);
cfg.channel                     = chan_pos_mag;
ERF_low_filt_sel_pos_GA_mag     = ft_selectdata(cfg,ERF_low_filt_GA);
ERF_high_filt_sel_pos_GA_mag    = ft_selectdata(cfg,ERF_high_filt_GA);
cfg.channel                     = chan_pos_cmb;
ERF_low_filt_sel_pos_GA_cmb     = ft_selectdata(cfg,ERF_low_filt_GA);
ERF_high_filt_sel_pos_GA_cmb    = ft_selectdata(cfg,ERF_high_filt_GA);

%% statistics SECOND COMPONENT

neigh_cmb = load('neuromag306cmb_neighb');
neigh_mag = load('neuromag306mag_neighb');
pad = 0.01;

% statistics
Nsub                    = length(slist);
cfg                     = [];
cfg.latency             = [time_max_mag-pad time_max_mag+pad];
cfg.avgovertime         = 'yes';
cfg.avgoverchan         = 'yes';
cfg.method              = 'analytic';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.minnbchan           = 0;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.05;
cfg.numrandomization    = 2000;
cfg.neighbours          = neigh_mag.neighbours; %  ft_prepare_neighbours(cfg_neighb, grad);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

cfg.channel             = chan_pos_mag;
stat_pos_mag_2          = ft_timelockstatistics(cfg,ERF_low_filt{slist},ERF_high_filt{slist});

cfg.channel             = chan_neg_mag;
stat_neg_mag_2          = ft_timelockstatistics(cfg,ERF_low_filt{slist},ERF_high_filt{slist});

cfg.channel             = chan_pos_cmb;
stat_cmb_2              = ft_timelockstatistics(cfg,ERF_low_filt{slist},ERF_high_filt{slist});


%% plot overview with ERF pos, neg and topo SECOND COMPONENT
fig = figure;
cfg                     = [];
cfg.layout              = 'neuromag306mag';
cfg.xlim                = [time_max_mag-pad time_max_mag+pad];
cfg.zlim                = 'maxabs';
cfg.highlightchannel    = chan_neg_mag;
cfg.highlight           = 'yes';
cfg.highlightsymbol     = '.';
cfg.highlightcolor      = [0 0 0];
cfg.highlightsize       = 20;
cfg.marker              = 'off';
cfg.comment             = 'no';
cfg.colorbar            = 'yes';
cfg.colormap            = colormap(jet(1024));
cfg.gridscale           = 64;
ft_topoplotER(cfg,ERF_filt_GA);

% set(gca,'PlotBoxAspectRatio',[8/3,1,1]);
saveas(fig,'w:\WANDER\images\FIGd_offset_topo_2nd.eps','eps');
% 
% fig = figure;
% t = 1;
% stim_times = -0.5:1/16:0;
% 
% subplot(3,2,3);
% axis([-0.5,1,-10e-14,20e-14]); ax = axis; hold;
% patch([time_max_mag-pad time_max_mag+pad time_max_mag+pad time_max_mag-pad],[ax(3),ax(3) ax(4),ax(4)],[0 0 0 ],'FaceAlpha',.2,'LineStyle','None');
% patch([ERF_low_filt_sel_pos_GA_mag.time, ERF_low_filt_sel_pos_GA_mag.time(end:-1:1)],[ERF_low_filt_sel_pos_GA_mag.avg+ERF_low_filt_sel_pos_GA_mag.std ./sqrt(22)*t, ERF_low_filt_sel_pos_GA_mag.avg(end:-1:1)-ERF_low_filt_sel_pos_GA_mag.std(end:-1:1) ./sqrt(21)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
% patch([ERF_high_filt_sel_pos_GA_mag.time, ERF_high_filt_sel_pos_GA_mag.time(end:-1:1)],[ERF_high_filt_sel_pos_GA_mag.avg+ERF_high_filt_sel_pos_GA_mag.std ./sqrt(22)*t, ERF_high_filt_sel_pos_GA_mag.avg(end:-1:1)-ERF_high_filt_sel_pos_GA_mag.std(end:-1:1) ./sqrt(21)*t],[1 0 0 ],'FaceAlpha',.5,'LineStyle','None');
% 
% y1 = ones(size(stim_times))*ax(3);
% y2 = ones(size(stim_times))*ax(4);
% for i = 1 : length(stim_times)
%     line([stim_times(i),stim_times(i)],[y1(i), y1(i) + abs(y1(i))*0.1],'color',[0 0 0]);
% end
% title(sprintf('Positive channel, t=%f, p=%f\n',stat_pos_mag_2.stat,stat_pos_mag_2.prob));


%% Figure 1D
pad = 0.01;
fig = figure;
axis([-0.5,1,-20e-14,15e-14]); 
ax = axis; hold;
patch([time_max_mag-pad time_max_mag+pad time_max_mag+pad time_max_mag-pad],[ax(3),ax(3) ax(4),ax(4)],[0 0 0],'LineStyle','None');
patch([ERF_low_filt_sel_neg_GA_mag.time, ERF_low_filt_sel_neg_GA_mag.time(end:-1:1)],[ERF_low_filt_sel_neg_GA_mag.avg+ERF_low_filt_sel_neg_GA_mag.std ./sqrt(22)*t,     ERF_low_filt_sel_neg_GA_mag.avg(end:-1:1)-ERF_low_filt_sel_neg_GA_mag.std(end:-1:1)   ./sqrt(22)*t],[0 0.3 0],'LineStyle','None');
patch([ERF_high_filt_sel_neg_GA_mag.time, ERF_high_filt_sel_neg_GA_mag.time(end:-1:1)],[ERF_high_filt_sel_neg_GA_mag.avg+ERF_high_filt_sel_neg_GA_mag.std ./sqrt(22)*t, ERF_high_filt_sel_neg_GA_mag.avg(end:-1:1)-ERF_high_filt_sel_neg_GA_mag.std(end:-1:1) ./sqrt(22)*t],[1 0 0],'LineStyle','None');

line(ERF_low_filt_sel_neg_GA_mag.time, ERF_low_filt_sel_neg_GA_mag.avg);
line(ERF_low_filt_sel_neg_GA_mag.time, ERF_high_filt_sel_neg_GA_mag.avg);

y1 = ones(size(stim_times))*ax(3);
y2 = ones(size(stim_times))*ax(4);
for i = 1 : length(stim_times)
    line([stim_times(i),stim_times(i)],[y1(i), y1(i) + abs(y1(i))*0.1],'color',[0 0 0]);
end
% title(sprintf('Negative channel, t=%f, p=%f\n',stat_neg_mag_2.stat,stat_neg_mag_2.prob));
set(gca,'PlotBoxAspectRatio',[8/3,1,1]);
savefig(fig,'w:\WANDER\images\poster\FIGd_offset_timecourse');

%%







subplot(3,2,4);
axis([-0.5,1,0,8e-12]); ax = axis; hold;
patch([time_max_cmb-pad time_max_cmb+pad time_max_cmb+pad time_max_cmb-pad],[ax(3),ax(3) ax(4),ax(4)],[0 0 0 ],'FaceAlpha',.2,'LineStyle','None');
patch([ERF_low_filt_sel_pos_GA_cmb.time,  ERF_low_filt_sel_pos_GA_cmb.time(end:-1:1)], [ERF_low_filt_sel_pos_GA_cmb.avg  + ERF_low_filt_sel_pos_GA_cmb.std  ./sqrt(22)*t, ERF_low_filt_sel_pos_GA_cmb.avg(end:-1:1)-ERF_low_filt_sel_pos_GA_cmb.std(end:-1:1) ./sqrt(21)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
patch([ERF_high_filt_sel_pos_GA_cmb.time, ERF_high_filt_sel_pos_GA_cmb.time(end:-1:1)],[ERF_high_filt_sel_pos_GA_cmb.avg + ERF_high_filt_sel_pos_GA_cmb.std ./sqrt(22)*t, ERF_high_filt_sel_pos_GA_cmb.avg(end:-1:1)-ERF_high_filt_sel_pos_GA_cmb.std(end:-1:1) ./sqrt(21)*t],[1 0 0 ],'FaceAlpha',.5,'LineStyle','None');
y1 = ones(size(stim_times))*ax(3);
y2 = ones(size(stim_times))*ax(4);
for i = 1 : length(stim_times)
    line([stim_times(i),stim_times(i)],[0, abs(y2(i))*0.05],'color',[0 0 0]);
end
title(sprintf('Maximum channel, t=%f, p=%f\n',stat_cmb_1.stat,stat_cmb_1.prob));

set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0 11*4 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['d:\analysis\WANDER\images\ERF_overview_probe_second.pdf']);
print(fig,'-dpng',['d:\analysis\WANDER\images\ERF_overview_probe_second.png']);



%% Calculate TFR at SSSF frequency
for isubject = slist
    cfg                         = [];
    cfg.pad                     = 'nextpow2';
    cfg.channel                 = 'all';
    cfg.method                  = 'mtmconvol';
    cfg.foi                     = [1:30];
    cfg.taper                   = 'hanning';
    cfg.t_ftimwin               = ones(size(cfg.foi));
    cfg.toi                     = -9.5:0.05:0.45;
    
    TFR{isubject}           = ft_freqanalysis(cfg, ERF{isubject});
    TFR_high{isubject}      = ft_freqanalysis(cfg, ERF_high{isubject});
    TFR_low{isubject}       = ft_freqanalysis(cfg, ERF_low{isubject});
end

% combine planar
for isubject = slist
    TFR{isubject}      = ft_combineplanar([],TFR{isubject});
    TFR_high{isubject} = ft_combineplanar([],TFR_high{isubject});
    TFR_low{isubject}  = ft_combineplanar([],TFR_low{isubject});
end

% difference per subject
for isubject = slist
    TFR_diff{isubject}                  =  TFR_low{isubject};
    TFR_diff{isubject}.powspctrm        =  (TFR_low{isubject}.powspctrm - TFR_high{isubject}.powspctrm)  ./ (TFR_low{isubject}.powspctrm + TFR_high{isubject}.powspctrm);
    TFR_low_norm{isubject}              =  TFR_low{isubject};
    TFR_low_norm{isubject}.powspctrm    = (TFR_low{isubject}.powspctrm) ./ (TFR_low{isubject}.powspctrm + TFR_high{isubject}.powspctrm);    
    TFR_high_norm{isubject}             =  TFR_low{isubject};
    TFR_high_norm{isubject}.powspctrm   = (TFR_high{isubject}.powspctrm) ./ (TFR_low{isubject}.powspctrm + TFR_high{isubject}.powspctrm);        
end

% average over subjects
TFR_GA                  = ft_freqgrandaverage([],TFR{slist});
TFR_high_GA             = ft_freqgrandaverage([],TFR_high{slist});
TFR_low_GA              = ft_freqgrandaverage([],TFR_low{slist});
TFR_diff_GA             = ft_freqgrandaverage([],TFR_diff{slist});


% get maximum sensor
peak_latency_TFR        = [-9.5 -0.5];

cfg                     = [];
cfg.latency             = peak_latency_TFR;
cfg.avgovertime         = 'yes';
cfg.channel             = 'MEG*1';
cfg.frequency           = [16 16];
cfg.avgoverfreq         = 'yes';

temp                    = ft_selectdata(cfg,TFR_GA);
[Y I]                   = sort(temp.powspctrm,'ascend');
chan_TFR_mag            = ft_channelselection(I(end), temp.label);

cfg.channel             = 'MEG*3';
temp                    = ft_selectdata(cfg,TFR_GA);
[Y I]                   = sort(temp.powspctrm,'ascend');
chan_TFR_cmb            = ft_channelselection(I(end), temp.label);

% grand averages for plotting 16Hz
i = 1; clear TFR_ind*
for isubject = slist
    cfg                     = [];
    cfg.frequency           = [16 16];
    cfg.avgoverfreq         = 'yes';
    cfg.channel             = chan_TFR_mag;
    cfg.latency             = [-9.5 0.45];

    temp                    = ft_selectdata(cfg,TFR{isubject});
    TFR_ind_mag(i,:)        = squeeze(temp.powspctrm);
    temp                    = ft_selectdata(cfg,TFR_high{isubject});
    TFR_ind_high_mag(i,:)   = squeeze(temp.powspctrm);
    temp                    = ft_selectdata(cfg,TFR_low{isubject});
    TFR_ind_low_mag(i,:)    = squeeze(temp.powspctrm);
    temp                    = ft_selectdata(cfg,TFR_diff{isubject});
    TFR_ind_diff_mag(i,:)   = squeeze(temp.powspctrm);
    
    cfg.channel             = chan_TFR_cmb;
    temp                    = ft_selectdata(cfg,TFR{isubject});
    TFR_ind_cmb(i,:)        = squeeze(temp.powspctrm);
    temp                    = ft_selectdata(cfg,TFR_high{isubject});
    TFR_ind_high_cmb(i,:)   = squeeze(temp.powspctrm);
    temp                    = ft_selectdata(cfg,TFR_low{isubject});
    TFR_ind_low_cmb(i,:)    = squeeze(temp.powspctrm);    
    temp                    = ft_selectdata(cfg,TFR_diff{isubject});
    TFR_ind_diff_cmb(i,:)   = squeeze(temp.powspctrm);    
    i = i + 1;
end

cfg = [];
cfg.foilim                      = [16 16];
cfg.channel                     = chan_TFR_cmb;
TFR_16hz_cmb_GA                 = ft_freqgrandaverage(cfg,TFR{slist});
TFR_16hz_high_cmb_GA            = ft_freqgrandaverage(cfg,TFR_high{slist});
TFR_16hz_low_cmb_GA             = ft_freqgrandaverage(cfg,TFR_low{slist});
TFR_16hz_diff_cmb_GA            = ft_freqgrandaverage(cfg,TFR_diff{slist});

cfg.channel                     = chan_TFR_mag;
TFR_16hz_mag_GA                 = ft_freqgrandaverage(cfg,TFR{slist});
TFR_16hz_high_mag_GA            = ft_freqgrandaverage(cfg,TFR_high{slist});
TFR_16hz_low_mag_GA             = ft_freqgrandaverage(cfg,TFR_low{slist});
TFR_16hz_diff_mag_GA            = ft_freqgrandaverage(cfg,TFR_diff{slist});

% the rest is just to remove the first timepoints which are nan
TFR_16hz_cmb_GA.time            = temp.time;
TFR_16hz_high_cmb_GA.time       = temp.time;
TFR_16hz_low_cmb_GA.time        = temp.time;
TFR_16hz_mag_GA.time            = temp.time;
TFR_16hz_high_mag_GA.time       = temp.time;
TFR_16hz_low_mag_GA.time        = temp.time;

TFR_16hz_mag_GA.powspctrm       = nanmean(TFR_ind_mag);
TFR_16hz_high_mag_GA.powspctrm  = nanmean(TFR_ind_high_mag);
TFR_16hz_low_mag_GA.powspctrm   = nanmean(TFR_ind_low_mag);
TFR_16hz_diff_mag_GA.powspctrm  = nanmean(TFR_ind_diff_mag);

TFR_16hz_cmb_GA.powspctrm       = nanmean(TFR_ind_cmb);
TFR_16hz_high_cmb_GA.powspctrm  = nanmean(TFR_ind_high_cmb);
TFR_16hz_low_cmb_GA.powspctrm   = nanmean(TFR_ind_low_cmb);
TFR_16hz_diff_cmb_GA.powspctrm  = nanmean(TFR_ind_diff_cmb);

TFR_16hz_cmb_GA.std             = nanstd(TFR_ind_cmb,1);
TFR_16hz_high_cmb_GA.std        = nanstd(TFR_ind_high_cmb,1);
TFR_16hz_low_cmb_GA.std         = nanstd(TFR_ind_low_cmb,1);
TFR_16hz_diff_cmb_GA.std        = nanstd(TFR_ind_diff_cmb,1);

TFR_16hz_mag_GA.std             = nanstd(TFR_ind_mag,1);
TFR_16hz_high_mag_GA.std        = nanstd(TFR_ind_high_mag,1);
TFR_16hz_low_mag_GA.std         = nanstd(TFR_ind_low_mag,1);
TFR_16hz_diff_mag_GA.std        = nanstd(TFR_ind_diff_mag,1);

%% statistics TFR

neigh_cmb = load('neuromag306cmb_neighb');
neigh_mag = load('neuromag306mag_neighb');

% statistics
Nsub                    = length(slist);
cfg                     = [];
cfg.latency             = [-9.5 -0.5];
cfg.avgovertime         = 'yes';
cfg.avgoverchan         = 'yes';
cfg.frequency           = [16 16];

cfg.method              = 'analytic';
cfg.statistic           = 'ft_statfun_depsamplesT';
% cfg.correctm            = 'cluster';
% cfg.clusteralpha        = 0.05;
% cfg.clusterstatistic    = 'maxsum';
% cfg.minnbchan           = 1;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.minnbchan           = 1;
cfg.alpha               = 0.05;
cfg.numrandomization    = 2000;
cfg.neighbours          = neigh_mag.neighbours; %  ft_prepare_neighbours(cfg_neighb, grad);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

cfg.channel             = chan_TFR_mag;
stat_TFR_mag            = ft_freqstatistics(cfg,TFR_low_norm{slist},TFR_high_norm{slist});

cfg.channel             = chan_TFR_cmb;
stat_TFR_cmb            = ft_freqstatistics(cfg,TFR_low_norm{slist},TFR_high_norm{slist});

% 
% figure; 
% patch([TFR_16hz_high_mag_GA.time, TFR_16hz_high_mag_GA.time(end:-1:1)], [TFR_16hz_diff_cmb_GA.powspctrm+TFR_16hz_diff_cmb_GA.std ./sqrt(22)*t, TFR_16hz_diff_cmb_GA.powspctrm(end:-1:1)- TFR_16hz_diff_cmb_GA.std(end:-1:1) ./sqrt(22)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
% 
% 
% figure; 

%% Plot overview SSEF TFR results
fig = figure;

subplot(4,2,1);
cfg                     = [];
cfg.layout              = 'neuromag306mag';
cfg.xlim                = peak_latency_TFR;
cfg.zlim                = 'maxabs';
cfg.highlightchannel    = [chan_TFR_mag];
cfg.highlight           = 'yes';
cfg.highlightsymbol     = '.';
cfg.highlightcolor      = [0 0 0];
cfg.highlightsize       = 20;
cfg.marker              = 'off';
cfg.comment             = 'no';
cfg.colorbar            = 'yes';
cfg.ylim                = [16 16];
ft_topoplotTFR(cfg,TFR_GA);

subplot(4,2,2);
cfg                     = [];
cfg.layout              = 'neuromag306cmb';
cfg.xlim                = peak_latency_TFR;
cfg.zlim                = 'maxabs';
cfg.highlightchannel    = chan_TFR_cmb;
cfg.highlight           = 'yes';
cfg.highlightsymbol     = 'pentagram';
cfg.highlightcolor      = [0 0 0];
cfg.highlightsize       = 16;
cfg.marker              = 'off';
cfg.comment             = 'no';
cfg.colorbar            = 'yes';
cfg.ylim                = [16 16];
ft_topoplotTFR(cfg,TFR_GA);

subplot(4,2,3);
cfg                     = [];
cfg.layout              = 'neuromag306cmb';
cfg.xlim                = [-3 0.5];
cfg.zlim                = 'maxabs';
cfg.channel             = chan_TFR_mag;
cfg.baseline            = [0.4 0.5];

cfg.baselinetype        = 'relative';
cfg.interactive         = 'no';
cfg.ylim                = [1 30];
ft_singleplotTFR(cfg,TFR_GA);

subplot(4,2,4);

%% FIGURE
fig = figure;

cfg                     = [];
cfg.layout              = 'neuromag306cmb';
cfg.xlim                = [-10 0.5];
cfg.zlim                = 'maxabs';
cfg.channel             = chan_TFR_cmb;
cfg.baseline            = [0.4 0.5];
cfg.baselinetype        = 'relative';
cfg.interactive         = 'yes';
cfg.ylim                = [1 30];
ft_singleplotTFR(cfg,TFR_GA);
title('');
set(gca,'PlotBoxAspectRatio',[10,1,1]);
saveas(fig,'w:\WANDER\images\FIGe_offset_tfr.eps','eps');

for isubject = slist
    cfg = [];
    cfg.frequency               = [16 16];
    cfg.latency                 = [-10 0];
    cfg.channel                 = chan_TFR_cmb;
    cfg.avgovertime             = 'yes';
    TFR_low_sel{isubject}  = ft_selectdata(cfg,TFR_low{isubject});
    TFR_high_sel{isubject} = ft_selectdata(cfg,TFR_high{isubject});
    TFR_diff_sel{isubject}  = ft_selectdata(cfg,TFR_diff{isubject});
end

cfg = [];
cfg.keepindividual = 'yes';
TFR_low_sel_avg  = ft_freqgrandaverage(cfg,TFR_low_sel{slist});
TFR_high_sel_avg = ft_freqgrandaverage(cfg,TFR_high_sel{slist});
TFR_diff_sel_avg = ft_freqgrandaverage(cfg,TFR_diff_sel{slist});

figure; hold;
errorbar([mean(TFR_low_sel_avg.powspctrm); mean(TFR_high_sel_avg.powspctrm)],[std(TFR_low_sel_avg.powspctrm) std(TFR_high_sel_avg.powspctrm)]./ sqrt(22));
bar(mean(TFR_low_sel_avg.powspctrm));

% plot bar difference
fig = figure; hold;
errorbar(mean(TFR_diff_sel_avg.powspctrm)*100,std(TFR_diff_sel_avg.powspctrm)*100 ./ sqrt(22));
bar(mean(TFR_diff_sel_avg.powspctrm)*100);
ylim([0,14]);
xlim([0 2]);
saveas(fig,'w:\WANDER\images\FIGf_TFR_diff_bar.eps','eps');


%%

subplot(4,2,5);
figure;

stim_times = [-9.5:1/16:0];
y1 = ones(size(stim_times))*0;
y2 = ones(size(stim_times))*1.3e-27;
for i = 1 : length(stim_times)
    line([stim_times(i),stim_times(i)],[y1', y2'.*0.05],'color',[0 0 0]);
end
ax = axis;
patch([peak_latency_TFR, peak_latency_TFR(2:-1:1)],[0,0 1.3e-27,1.3e-27],[0 0 0 ],'FaceAlpha',.2,'LineStyle','None');
patch([TFR_16hz_high_mag_GA.time, TFR_16hz_high_mag_GA.time(end:-1:1)], [TFR_16hz_high_mag_GA.powspctrm+TFR_16hz_high_mag_GA.std ./sqrt(22)*t,    TFR_16hz_high_mag_GA.powspctrm(end:-1:1)- TFR_16hz_high_mag_GA.std(end:-1:1) ./sqrt(22)*t],[1 0 0 ],'FaceAlpha',.5,'LineStyle','None');
patch([TFR_16hz_low_mag_GA.time,  TFR_16hz_low_mag_GA.time(end:-1:1)],  [TFR_16hz_low_mag_GA.powspctrm+TFR_16hz_low_mag_GA.std ./sqrt(22)*t,      TFR_16hz_low_mag_GA.powspctrm(end:-1:1)-  TFR_16hz_low_mag_GA.std(end:-1:1)  ./sqrt(22)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
axis tight

% subplot(4,2,6);

stim_times = [-9.5:1/16:0];
y1 = ones(size(stim_times))*0;
y2 = ones(size(stim_times))*1.7e-24;
for i = 1 : length(stim_times)
    line([stim_times(i),stim_times(i)],[y1', y2'.*0.05],'color',[0 0 0]);
end
ax = axis;
patch([peak_latency_TFR, peak_latency_TFR(2:-1:1)],[0,0 1.7e-24,1.7e-24],[0 0 0 ],'FaceAlpha',.2,'LineStyle','None');
patch([TFR_16hz_high_cmb_GA.time, TFR_16hz_high_cmb_GA.time(end:-1:1)], [TFR_16hz_high_cmb_GA.powspctrm+TFR_16hz_high_cmb_GA.std ./sqrt(22)*t,    TFR_16hz_high_cmb_GA.powspctrm(end:-1:1)- TFR_16hz_high_cmb_GA.std(end:-1:1) ./sqrt(22)*t],[1 0 0 ],'FaceAlpha',.5,'LineStyle','None');
patch([TFR_16hz_low_cmb_GA.time,  TFR_16hz_low_cmb_GA.time(end:-1:1)],  [TFR_16hz_low_cmb_GA.powspctrm+TFR_16hz_low_cmb_GA.std ./sqrt(22)*t,      TFR_16hz_low_cmb_GA.powspctrm(end:-1:1)-  TFR_16hz_low_cmb_GA.std(end:-1:1)  ./sqrt(22)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
axis tight

subplot(4,2,7); hold;
line(TFR_16hz_high_mag_GA.time,zeros(size(TFR_16hz_diff_mag_GA.powspctrm,2)),'color','k','linestyle','--');
patch([TFR_16hz_high_mag_GA.time, TFR_16hz_high_mag_GA.time(end:-1:1)], [TFR_16hz_diff_mag_GA.powspctrm+TFR_16hz_diff_mag_GA.std ./sqrt(22)*t, TFR_16hz_diff_mag_GA.powspctrm(end:-1:1)- TFR_16hz_diff_mag_GA.std(end:-1:1) ./sqrt(22)*t],[0 0 0 ],'FaceAlpha',.5,'LineStyle','None');
plot(TFR_16hz_high_mag_GA.time,TFR_16hz_diff_mag_GA.powspctrm,'k');
scatter(TFR_16hz_high_mag_GA.time(stat_TFR_mag.mask),TFR_16hz_diff_mag_GA.powspctrm(stat_TFR_mag.mask),8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
axis tight

subplot(4,2,8); hold;
line(TFR_16hz_high_mag_GA.time,zeros(size(TFR_16hz_diff_mag_GA.powspctrm,2)),'color','k','linestyle','--');
patch([TFR_16hz_high_mag_GA.time, TFR_16hz_high_mag_GA.time(end:-1:1)], [TFR_16hz_diff_cmb_GA.powspctrm+TFR_16hz_diff_cmb_GA.std ./sqrt(22)*t, TFR_16hz_diff_cmb_GA.powspctrm(end:-1:1)- TFR_16hz_diff_cmb_GA.std(end:-1:1) ./sqrt(22)*t],[0 0 0 ],'FaceAlpha',.5,'LineStyle','None');
plot(TFR_16hz_high_mag_GA.time,TFR_16hz_diff_cmb_GA.powspctrm,'k');
scatter(TFR_16hz_high_cmb_GA.time(stat_TFR_cmb.mask),TFR_16hz_diff_cmb_GA.powspctrm(stat_TFR_cmb.mask),8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
axis tight

% set(fig,'PaperSize',[11*2 8.5*2 ]);
% set(fig,'PaperPosition',[0 0 11*4 8.5*2]);
% set(fig,'PaperOrientation','landscape');
% set(fig,'Position',[50 50 1200 800]);
% print(fig,'-dpdf',['d:\analysis\WANDER\images\SSTFR_overview_probe.pdf']);
% print(fig,'-dpng',['d:\analysis\WANDER\images\SSTFR_overview_probe.png']);
% 
% 
% 

% grand average and difference steady-state over time
fig = figure;
subplot(2,1,1);

cfg                     = [];
cfg.layout              = 'neuromag306cmb';
cfg.xlim                = [-10 0.5];
cfg.zlim                = 'maxabs';
cfg.channel             = chan_TFR_cmb;
cfg.baseline            = [0.4 0.5];
cfg.baselinetype        = 'relative';
cfg.interactive         = 'yes';
cfg.ylim                = [1 30];
cfg.zlim                = [-50 50];
cfg.colorbar = 'yes';
ft_singleplotTFR(cfg,TFR_GA);
title('');

subplot(2,1,2);
hold;
line(TFR_16hz_high_mag_GA.time,zeros(size(TFR_16hz_diff_mag_GA.powspctrm,2)),'color','k','linestyle','--');
patch([TFR_16hz_high_mag_GA.time, TFR_16hz_high_mag_GA.time(end:-1:1)], [TFR_16hz_diff_cmb_GA.powspctrm+TFR_16hz_diff_cmb_GA.std ./sqrt(22)*t, TFR_16hz_diff_cmb_GA.powspctrm(end:-1:1)- TFR_16hz_diff_cmb_GA.std(end:-1:1) ./sqrt(22)*t],[0 0 0 ],'FaceAlpha',.5,'LineStyle','None');
plot(TFR_16hz_high_mag_GA.time,TFR_16hz_diff_cmb_GA.powspctrm,'k');
scatter(TFR_16hz_high_cmb_GA.time(stat_TFR_cmb.mask),TFR_16hz_diff_cmb_GA.powspctrm(stat_TFR_cmb.mask),8,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);
axis tight;
ylim([-0.25,0.25]);

set(gca,'PlotBoxAspectRatio',[10,1,1]);
saveas(fig,'w:\WANDER\images\TFR_SSEF.eps','eps');
print(fig,'-dpdf',['w:\WANDER\images\TFR_SSEF.pdf']);

fig = figure;

cfg                     = [];
cfg.layout              = 'neuromag306cmb';
cfg.xlim                = peak_latency_TFR;
cfg.zlim                = 'maxabs';
cfg.highlightchannel    = chan_TFR_cmb;
cfg.highlight           = 'yes';
cfg.highlightsymbol     = '.';
cfg.highlightcolor      = [0 0 0];
cfg.highlightsize       = 20;
cfg.marker              = 'off';
cfg.comment             = 'no';
cfg.colorbar            = 'no';
% cfg.baseline            = [0.4 0.5];
% cfg.baselinetype        = 'relative';
cfg.ylim                = [16 16];
cfg.gridscale           = 64;
% TFR_GA = rmfield(TFR_GA,'cfg');
ft_topoplotTFR(cfg,TFR_GA);

savefig('w:\WANDER\images\TFR_SSEF_topo');


saveas(fig,'w:\WANDER\images\TFR_SSEF_topo.eps','eps');
print(fig,'-dpdf',['w:\WANDER\images\TFR_SSEF_topo.pdf']);
