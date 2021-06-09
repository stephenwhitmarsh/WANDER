function WANDER_ERF_probe_GA(force,timing,rootpath)

if ~(strcmp(timing,'cue') || strcmp(timing,'probe'))
    fprintf('Use "cue" or "probe" as third argument.\n');
    return
end

timing = 'cue';
clear ERF_cue*
slist = [1:20 22:26];

for isubject = slist
    fprintf('Loading subject %d \n',isubject);    
    temp = load(['d:\analysis\WANDER\data\ERF\s' num2str(isubject) '_' timing '.mat'],'ERF_probe');
    ERF_probe{isubject} = temp.ERF_probe;
    
    cfg. 
end
clear temp
ERF_probe_GA = ft_timelockgrandaverage([],ERF_probe{slist});
 
for isubject = slist
    fprintf('Loading subject %d \n',isubject);    
    temp = load(['d:\analysis\WANDER\data\ERF\s' num2str(isubject) '_' timing '.mat'],'ERF_probe_high');
    ERF_probe_high{isubject} = temp.ERF_probe_high;
end
ERF_probe_high_GA = ft_timelockgrandaverage([],ERF_probe_high{slist});

for isubject = slist
    fprintf('Loading subject %d \n',isubject);    
    temp = load(['d:\analysis\WANDER\data\ERF\s' num2str(isubject) '_' timing '.mat'],'ERF_probe_low');
    ERF_probe_low{isubject} = temp.ERF_probe_low;
end
ERF_probe_low_GA = ft_timelockgrandaverage([],ERF_probe_low{slist});


cfg         = [];
cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEG*1';
cfg.xlim = [-3 0];
cfg.ylim = 'maxabs';
cfg.marker = 'off';
cfg.comment = 'no';
ft_topoplotER(cfg,ERF_probe_low_GA,ERF_probe_high_GA);

% 
% ERF_cue_diff_GA = ERF_cue_low_GA;
% ERF_cue_diff_GA.avg = ERF_cue_low_GA.avg(:,1:20000) - ERF_cue_high_GA.avg(:,1:20000);
% ERF_cue_diff_GA.time = ERF_cue_diff_GA.time(1:20000);

% extract median split info
trialinfo_low = nan(length(slist),7);
trialinfo_high = nan(length(slist),7);
for isubject = slist
    fprintf('Loading subject %d \n',isubject);    
    trialinfo_high(isubject,:) = ERF_probe_high{isubject}.trialinfo;
    trialinfo_low(isubject,:) = ERF_probe_low{isubject}.trialinfo;  
end
nr_ratings_high = sum(trialinfo_high,2);
nr_ratings_low  = sum(trialinfo_low,2);
nr_ratings_high - nr_ratings_low;

nr_ratings_high_avg = nanmean(nr_ratings_high);
nr_ratings_low_avg = nanmean(nr_ratings_low);
nr_ratings_high_std = nanstd(nr_ratings_high);
nr_ratings_low_std = nanstd(nr_ratings_low);
errorbar([nr_ratings_high_avg,nr_ratings_low_avg],[nr_ratings_high_std,nr_ratings_low_std]);
xlim([0 3]);

% plot average rating distribution between high and low
trialinfo_low_avg = nanmean(trialinfo_low);
trialinfo_high_avg = nanmean(trialinfo_high);
trialinfo_low_std = nanstd(trialinfo_low);
trialinfo_high_std = nanstd(trialinfo_high);
figure; bar([trialinfo_low_avg; trialinfo_high_avg]',0.9); legend({'high','low'});
figure; errorbar([trialinfo_low_avg; trialinfo_high_avg]',[trialinfo_low_std/5; trialinfo_high_std/5]'); legend({'high (SEM)','low (SEM)'});

% shorten trials to have same length over subjects
maxlength = 0;
minlength = 30000;
for isubject = slist
    if size(ERF_probe{isubject}.time,2) > maxlength
        maxlength = size(ERF_probe{isubject}.time,2);
    end
    if size(ERF_probe{isubject}.time,2) < minlength
        minlength = size(ERF_probe{isubject}.time,2);
    end
end

for isubject = slist
    fprintf('Cutting subject %d to size \n',isubject);    
    cfg = [];
    cfg.latency = [-minlength/1000 0];
    ERF_probe_low{isubject}     = ft_selectdata(cfg,ERF_probe_low{isubject});
    ERF_probe_high{isubject}    = ft_selectdata(cfg,ERF_probe_high{isubject});
    ERF_probe{isubject}         = ft_selectdata(cfg,ERF_probe{isubject});
end

% Need to cut off end of trial first !!!!!!!!!!!!!!
% Calculate FFT
for isubject = slist
    cfg                         = [];
    cfg.pad                     = 'nextpow2';
    cfg.channel                 = 'all';
    cfg.method                  = 'mtmfft';
    cfg.foi                     = [1:100];
    cfg.taper                   = 'hanning';
    FFT_probe{isubject}         = ft_freqanalysis(cfg, ERF_probe{isubject});
    FFT_probe_high{isubject}    = ft_freqanalysis(cfg, ERF_probe_high{isubject});
    FFT_probe_low{isubject}     = ft_freqanalysis(cfg, ERF_probe_low{isubject});
    
    % combine planar 
    FFT_probe{isubject}      = ft_combineplanar([],FFT_probe{isubject});
    FFT_probe_high{isubject} = ft_combineplanar([],FFT_probe_high{isubject});
    FFT_probe_low{isubject}  = ft_combineplanar([],FFT_probe_low{isubject});       
end

    
FFT_probe_GA        = ft_freqgrandaverage([],FFT_probe{slist});
FFT_probe_high_GA   = ft_freqgrandaverage([],FFT_probe_high{slist});
FFT_probe_low_GA    = ft_freqgrandaverage([],FFT_probe_low{slist});

% plot topo FFT
cfg         = [];
cfg.layout  = 'neuromag306cmb';
cfg.channel = 'MEG*3';
% cfg.xlim    = [-20 0];
% cfg.zlim    = 'maxabs';
% cfg.marker  = 'off';
cfg.comment = 'no';
figure; ft_singleplotER(cfg,FFT_probe_low_GA,FFT_probe_high_GA);


% Calculate TFR at SSSF frequency
for isubject = slist
    cfg                         = [];
    cfg.pad                     = 'nextpow2';
    cfg.channel                 = 'all';
    cfg.method                  = 'mtmconvol';
    cfg.foi                     = 16;
    cfg.taper                   = 'hanning';
    cfg.t_ftimwin               = 1;
    cfg.toi                     = -minlength:0.25:0;
    TFR_probe{isubject}         = ft_freqanalysis(cfg, ERF_probe{isubject});
    TFR_probe_high{isubject}    = ft_freqanalysis(cfg, ERF_probe_high{isubject});
    TFR_probe_low{isubject}     = ft_freqanalysis(cfg, ERF_probe_low{isubject});
    
    % combine planar 
    TFR_probe{isubject}      = ft_combineplanar([],TFR_probe{isubject});
    TFR_probe_high{isubject} = ft_combineplanar([],TFR_probe_high{isubject});
    TFR_probe_low{isubject}  = ft_combineplanar([],TFR_probe_low{isubject});    
end


% difference per subject
for isubject = slist
    TFR_probe_diff{isubject} =  TFR_probe_low{isubject};
    TFR_probe_diff{isubject}.powspctrm = (TFR_probe_high{isubject}.powspctrm - TFR_probe_low{isubject}.powspctrm) ./ (TFR_probe_low{isubject}.powspctrm - TFR_probe_high{isubject}.powspctrm);
end

% average over subjects
TFR_probe_GA                  = ft_freqgrandaverage([],TFR_probe{slist});
TFR_probe_high_GA             = ft_freqgrandaverage([],TFR_probe_high{slist});
TFR_probe_low_GA              = ft_freqgrandaverage([],TFR_probe_low{slist});
TFR_probe_diff_GA             = TFR_probe_high_GA;
TFR_probe_diff_GA.powspctrm   = TFR_probe_low_GA.powspctrm - TFR_probe_high_GA.powspctrm;

% plot topos average
cfg         = [];
cfg.layout  = 'neuromag306cmb';
cfg.xlim    = [-20 0];
cfg.zlim    = 'maxabs';
% cfg.marker  = 'off';
cfg.comment = 'no';
figure; ft_topoplotTFR(cfg,TFR_probe_diff_GA);

% determine top sensors for SSSF
for isubject = slist
    cfg = [];
    cfg.latency = [-10 -2];
    cfg.channel = 'MEG*3';
    cfg.avgovertime = 'yes';
    temp = ft_selectdata(cfg,TFR_probe{isubject});
    [y indx] = sort(temp.powspctrm,'descend');
    for isens = 1 : 9
        TFR_probe{isubject}.maxsens{isens} = temp.label{indx(isens)};
    end
    
    cfg = [];
    cfg.channel                             = TFR_probe{isubject}.maxsens; 
    cfg.avgoverchan                         = 'yes';
    TFR_probe_sel{isubject}                   = ft_selectdata(cfg,TFR_probe{isubject});
    TFR_probe_sel{isubject}.label{1}          = 'max';
    TFR_probe_low_sel{isubject}               = ft_selectdata(cfg,TFR_probe_low{isubject});
    TFR_probe_low_sel{isubject}.label{1}      = 'max';
    TFR_probe_high_sel{isubject}              = ft_selectdata(cfg,TFR_probe_high{isubject});
    TFR_probe_high_sel{isubject}.label{1}     = 'max';
    TFR_probe_diff_sel{isubject}              = TFR_probe_low_sel{isubject};
    TFR_probe_diff_sel{isubject}.powspctrm    = TFR_probe_low_sel{isubject}.powspctrm ./ (TFR_probe_high_sel{isubject}.powspctrm + TFR_probe_low_sel{isubject}.powspctrm);
%     TFR_cue_diff_sel{isubject}.powspctrm    = TFR_cue_low_sel{isubject}.powspctrm - TFR_cue_high_sel{isubject}.powspctrm;
end

% % add maximum sensor as separate sensor
% for isubject = slist
%     cfg = [];
%     cfg.channel = TFR_probe{isubject}.maxsens;
%     cfg.avgoverchan = 'yes';
%     ERF_probe_maxsens = ft_selectdata(cfg,ERF_probe{isubject});
%     ERF_probe_maxsens.label = 'max';
%     ERF_probe{isubject} = test = ft_appenddata(ERF_probe{isubject},ERF_probe_maxsens);
% end

cfg = [];
cfg.keepindividual      = 'yes';
TFR_probe_sel_GA          = ft_freqgrandaverage(cfg,TFR_probe_sel{slist});
TFR_probe_high_sel_GA     = ft_freqgrandaverage(cfg,TFR_probe_high_sel{slist});
TFR_probe_low_sel_GA      = ft_freqgrandaverage(cfg,TFR_probe_low_sel{slist});
TFR_probe_diff_sel_GA     = ft_freqgrandaverage(cfg,TFR_probe_diff_sel{slist});

TFR_probe_sel_GA.avg      = squeeze(mean(TFR_probe_sel_GA.powspctrm,1));
TFR_probe_sel_GA.std      = squeeze(std(TFR_probe_sel_GA.powspctrm,1));
TFR_probe_high_sel_GA.avg = squeeze(mean(TFR_probe_high_sel_GA.powspctrm,1));
TFR_probe_high_sel_GA.std = squeeze(std(TFR_probe_high_sel_GA.powspctrm,1));
TFR_probe_low_sel_GA.avg  = squeeze(mean(TFR_probe_low_sel_GA.powspctrm,1));
TFR_probe_low_sel_GA.std  = squeeze(std(TFR_probe_low_sel_GA.powspctrm,1));
TFR_probe_diff_sel_GA.avg = squeeze(mean(TFR_probe_diff_sel_GA.powspctrm,1));
TFR_probe_diff_sel_GA.std = squeeze(std(TFR_probe_diff_sel_GA.powspctrm,1));

% ttest2(squeeze(TFR_probe_high_sel_GA.powspctrm),squeeze(TFR_probe_high_sel_GA.powspctrm))


figure; hold;
errorbar(TFR_probe_high_sel_GA.time,        squeeze(TFR_probe_sel_GA.avg),      squeeze(TFR_probe_sel_GA.std));
errorbar(TFR_probe_high_sel_GA.time+0.001,  squeeze(TFR_probe_high_sel_GA.avg), squeeze(TFR_probe_high_sel_GA.std));
errorbar(TFR_probe_high_sel_GA.time+0.002,  squeeze(TFR_probe_low_sel_GA.avg),  squeeze(TFR_probe_low_sel_GA.std));

figure; hold;
errorbar(TFR_probe_diff_sel_GA.time,        squeeze(TFR_probe_diff_sel_GA.avg),      squeeze(TFR_probe_diff_sel_GA.std/5));

figure; hold;
errorbar(TFR_probe_diff_sel_GA.time,        squeeze(TFR_probe_diff_sel_GA.avg),      squeeze(TFR_probe_diff_sel_GA.std/5));

figure; hold;
% plot(TFR_probe_high_sel_GA.time,TFR_probe_sel_GA.avg);
plot(TFR_probe_high_sel_GA.time,TFR_probe_high_sel_GA.avg);
plot(TFR_probe_high_sel_GA.time,TFR_probe_low_sel_GA.avg);
plot(TFR_probe_diff_sel_GA.time,squeeze(stat.mask*4e-25));

figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    plot(TFR_probe_diff_sel{isubject}.time,squeeze(TFR_probe_diff_sel{isubject}.powspctrm));
    title(num2str(isubject));
    ylim([0 1]);
    i = i + 1;
end


% plot individual topos at SSSF frequency
i = 1;
figure;
for isubject = slist
    subplot(5,5,i);
    cfg         = [];
    cfg.layout  = 'neuromag306cmb';
    cfg.xlim = [-10 0];
    cfg.zlim = 'maxabs';
    cfg.highlight           = 'on';
    cfg.highlightchannel    = TFR_probe{isubject}.maxsens;
    cfg.highlightsymbol     = '.';
    cfg.highlighsize        = 13;
    cfg.highlightcolor      = [1 1 1];
    cfg.marker = 'off';
    cfg.comment = 'no';
    ft_topoplotTFR(cfg,TFR_probe{isubject});
    title(num2str(isubject));
    i = i + 1;
end

% plot individual topos of difference in SSSF power
i = 1;
figure;
for isubject = slist
    subplot(5,5,i);
    cfg         = [];
    cfg.layout  = 'neuromag306cmb';
    cfg.xlim = [-10 0];
    cfg.zlim = 'maxabs';
    cfg.marker = 'off';
    cfg.comment = 'no';
    cfg.highlight           = 'on';    
    cfg.highlightchannel    = TFR_probe{isubject}.maxsens;
    cfg.highlightsymbol     = '.';  
    cfg.highlighsize        = 13;
    cfg.highlightcolor      = [1 1 1];
    ft_topoplotTFR(cfg,TFR_probe_diff{isubject});
    title(num2str(isubject));
    i = i + 1;
end


load('D:\analysis\WANDER\scripts\grad','grad');
neigh = load('neuromag306cmb_neighb');

cfg = [];
cfg.grad = grad;
cfg.neighbours          = neigh.neighbours;
ft_neighbourplot(cfg);

% statistics
Nsub                    = length(slist);
cfg                     = [];
cfg.channel             = 'max';
% cfg.latency             = [15 21];
% cfg.latency             = [10 11];
cfg.avgoverfreq         = 'yes';
cfg.avgovertime         = 'no';
cfg.avgoverchan         = 'no';
cfg.parameter           = 'powspctrm';
cfg.method              = 'analytic';
cfg.statistic           = 'ft_statfun_depsamplesT';
% cfg.correctm            = 'cluster';
% cfg.clusteralpha        = 0.05;
% cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.05;
% cfg.numrandomization    = 2000;
cfg.neighbours          = neigh.neighbours; %  ft_prepare_neighbours(cfg_neighb, grad);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_freqstatistics(cfg,TFR_probe_high_sel{slist},TFR_probe_low_sel{slist});




























i = 1;
figure;
for isubject = slist
    subplot(5,5,i);
    cfg         = [];
    cfg.layout  = 'neuromag306mag';
    %     cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143'};
    cfg.channel = {'MEG1111', 'MEG1121', 'MEG1131', 'MEG1141'};
    %     cfg.channel = 'MEG*1';
    cfg.xlim = [-0.5 1];
    cfg.zlim = 'maxabs';
    ft_singleplotER(cfg,ERF_cue{isubject},ERF_cue_low{isubject},ERF_cue_high{isubject});
    title(num2str(isubject));
    i = i + 1;
end

% plot ERF grand-average
cfg         = [];
cfg.layout  = 'neuromag306mag';
cfg.channel = {'MEG1111', 'MEG1121', 'MEG1131', 'MEG1141'};
cfg.channel = 'MEG*1';
cfg.xlim = [-1 3];
cfg.zlim = 'maxabs';
figure; ft_singleplotER(cfg,ERF_cue_GA,ERF_cue_low_GA,ERF_cue_high_GA);

% Calculate TFR
cfg                         = [];
cfg.pad                     = 'nextpow2';
cfg.channel                 = 'all';
cfg.method                  = 'mtmconvol';
cfg.foi                     = 1:1:30;
cfg.taper                   = 'hanning';
cfg.t_ftimwin               = ones(size(cfg.foi));
cfg.toi                     = -1:0.05:20;
TFR_cue_GA                  = ft_freqanalysis(cfg, ERF_cue_GA);
TFR_cue_high_GA             = ft_freqanalysis(cfg, ERF_cue_high_GA);
TFR_cue_low_GA              = ft_freqanalysis(cfg, ERF_cue_low_GA);

% combine planar
TFR_cue_GA      = ft_combineplanar([],TFR_cue_GA);
TFR_cue_high_GA = ft_combineplanar([],TFR_cue_high_GA);
TFR_cue_low_GA  = ft_combineplanar([],TFR_cue_low_GA);

TFR_cue_diff_GA = TFR_cue_GA;
TFR_cue_diff_GA.powspctrm = TFR_cue_low_GA.powspctrm - TFR_cue_high_GA.powspctrm;

% plot TFR
fig             = figure;
cfg             = [];
cfg.layout      = 'neuromag306cmb';
% cfg.layout      = 'neuromag306mag';
% cfg.zlim        = [0.85 1.15];
cfg.baselinetype = 'relative';
cfg.zlim = 'maxabs';
cfg.baseline = [-1 0];
cfg.xlim        = [-1 10];
figure; ft_toptTFR(cfg,TFR_cue_GA);

figure; ft_singleplotTFR(cfg,TFR_cue_diff_GA);



% plot difference
cfg         = [];
cfg.layout  = 'neuromag306mag';
%     cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143'};
cfg.channel = 'MEG*1';
%     cfg.ylim    = [5 30];
cfg.xlim = [-1 10];
% % cfg.baselinetype = 'relative';
cfg.zlim = 'maxabs';
% cfg.baseline = [-10 1];
figure; ft_singleplotER(cfg,ERF_cue_diff_GA);


    
    
% statistics
cfg = [];
cfg.channel     = 'all';
% cfg.latency     = [-5 -1];
cfg.avgovertime = 'no';
cfg.frequency   = [8 14];
cfg.avgoverfreq = 'yes';
cfg.avgoverchan = 'no';
cfg.parameter   = 'powspctrm';
cfg.method           = 'montecarlo';
cfg.method           = 'analytic';
cfg.statistic           = 'ft_statfun_depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05;
% cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb, TFR{1});

Nsub = 25;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat = ft_freqstatistics(cfg,TFR_comb_high{[1:20 22:end]},TFR_comb_low{[1:20 22:end]});


cfg             = [];
cfg.layout      = 'neuromag306cmb';
figure; ft_clusterplot(cfg,stat)

MI_diff = MI_high_GA{1};
MI_diff.avg = MI_high_GA{1}.avg - MI_low_GA{1}.avg;


    %
    cfg         = [];
    cfg.layout  = 'neuromag306mag';
%     cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143'};
    cfg.channel = 'MEG*1';
%     cfg.ylim    = [5 30];
    cfg.xlim = [-1 3];
    % % cfg.baselinetype = 'relative';
    cfg.zlim = 'absmax';
    % cfg.baseline = [-10 1];
    figure; ft_singleplotER(cfg,ERF_cue_GA);



% TFR multiplot average
fig             = figure;
cfg             = [];
cfg.layout      = 'neuromag306cmb';
% cfg.layout      = 'neuromag306mag';
cfg.zlim        = [0.85 1.15];
cfg.xlim        = [-10 0];
figure; ft_singleplotTFR(cfg,TFR_comb_diff_GA);

cfg.parameter = 'stat';
figure; ft_singleplotTFR(cfg,stat);

% 9, 10, 18 is atypical, 11,12, 25 somewhat
figure; ft_singleplotTFR(cfg,TFR_comb_diff{2});


for i = [1:20 22:26]
    subplot(5,5,i);
    cfg             = [];
    cfg.layout      = 'neuromag306cmb';
    cfg.zlim        = 'maxabs';
    
    ft_singleplotTFR(cfg,TFR_comb_diff_GA);
end

% ax = axis;
% perc = squeeze(TFR_corr.nr) ./ max(TFR_corr.nr) * ax(4);
% plot(TFR_corr.time,perc,'r','linewidth',2);
% print(fig,'-dpdf',['d:\analysis\WANDER\images\TFR\s' num2str(isubject) '_' timing '_TFR_corr_cue.pdf']);
% print(fig,'-dpng',['d:\analysis\WANDER\images\TFR\s' num2str(isubject) '_' timing '_TFR_corr_cue.png']);


% Calculate TFR
cfg                         = [];
cfg.pad                     = 'nextpow2';
cfg.channel                 = 'all';
cfg.method                  = 'mtmconvol';
cfg.foi                     = 1:1:30;
cfg.taper                   = 'hanning';
cfg.t_ftimwin               = ones(size(cfg.foi));

cfg.toi                     = -1:0.05:30;
TFR_cue                     = ft_freqanalysis(cfg, ERF_cue_GA);

TFR_cue = ft_combineplanar([],TFR_cue);
    
% TFR multiplot average
fig             = figure;
cfg             = [];
cfg.layout      = 'neuromag306cmb';
% cfg.layout      = 'neuromag306mag';
% cfg.zlim        = [0.85 1.15];
cfg.xlim        = [-1 10];
figure; ft_singleplotTFR(cfg,TFR_cue);
