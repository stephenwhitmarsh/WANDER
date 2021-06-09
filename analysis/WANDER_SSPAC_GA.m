function WANDER_SSPAC_GA(isubject,force,timing,rootpath)

slist = [1:20 22:26]; %task
slist = [1:15 17:20 22:26]; %task

for isubject = slist
    temp = load(['i:\analysis\WANDER\data\SSPAC\s' num2str(isubject) '_SSPAC_cue.mat']);
    
    MI_all{isubject}        = temp.MI;
    MI_high{isubject}       = temp.MI_high;
    MI_low{isubject}        = temp.MI_low;
    
    MI_diff{isubject}       = MI_high{isubject};
    MI_diff{isubject}.avg   = (MI_low{isubject}.avg - MI_high{isubject}.avg) ./ (MI_low{isubject}.avg + MI_high{isubject}.avg);
   
    PAC_all{isubject}       = temp.PAC;
    PAC_high{isubject}      = temp.PAC_high;
    PAC_low{isubject}       = temp.PAC_low;
    
    PAC_all_norm{isubject}              = PAC_all{isubject};
    PAC_all_norm{isubject}.powspctrm    = PAC_all{isubject}.powspctrm ./  mean(PAC_all{isubject}.powspctrm,3);      
    PAC_high_norm{isubject}             = PAC_high{isubject};
    PAC_high_norm{isubject}.powspctrm   = PAC_high_norm{isubject}.powspctrm ./  mean(PAC_high_norm{isubject}.powspctrm,3);  
    PAC_low_norm{isubject}              = PAC_low{isubject};    
    PAC_low_norm{isubject}.powspctrm    = PAC_low_norm{isubject}.powspctrm ./  mean(PAC_low_norm{isubject}.powspctrm,3);
    
    PAC_diff{isubject}                  = PAC_high{isubject};
    PAC_diff{isubject}.powspctrm        = (PAC_low{isubject}.powspctrm - PAC_high{isubject}.powspctrm) ./ (PAC_low{isubject}.powspctrm + PAC_high{isubject}.powspctrm);
    PAC_norm_diff{isubject}             = PAC_high{isubject};
    PAC_norm_diff{isubject}.powspctrm   = (PAC_low_norm{isubject}.powspctrm - PAC_high_norm{isubject}.powspctrm);
    
    PAC_diff_abs{isubject}                  = PAC_high{isubject};
    PAC_diff_abs{isubject}.powspctrm        = (PAC_low{isubject}.powspctrm ./ PAC_high{isubject}.powspctrm);
end

MI_all_GA = ft_timelockgrandaverage([],MI_all{slist});
MI_high_GA = ft_timelockgrandaverage([],MI_high{slist});
MI_low_GA = ft_timelockgrandaverage([],MI_low{slist});
MI_diff_GA = ft_timelockgrandaverage([],MI_diff{slist});

PAC_all_GA = ft_freqgrandaverage([],PAC_all{slist});
PAC_all_norm_GA = ft_freqgrandaverage([],PAC_all_norm{slist});
PAC_high_GA = ft_freqgrandaverage([],PAC_high{slist});
PAC_low_GA = ft_freqgrandaverage([],PAC_low{slist});
PAC_diff_abs_GA = ft_freqgrandaverage([],PAC_diff_abs{slist});
PAC_norm_diff_GA = ft_freqgrandaverage([],PAC_norm_diff{slist});


cfg                 = [];
cfg.layout          = 'neuromag306cmb.lay';
cfg.channel         = 'MEG*3';
cfg.channel         = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};
% cfg.zlim            = 'maxabs';
% cfg.zlim = [0.992 1.008];
figure; ft_singleplotTFR(cfg, PAC_all_norm_GA);
figure; ft_singleplotTFR(cfg, PAC_diff_GA);

figure; ft_singleplotTFR(cfg, PAC_norm_diff_GA);

% cfg.ylim            = 'maxabs';
figure; ft_singleplotER(cfg, MI_diff{1});
figure; ft_singleplotER(cfg, MI_all_GA);
figure; ft_singleplotER(cfg, MI_diff_GA);
figure; ft_singleplotER(cfg, MI_low_GA, MI_high_GA);


% plot PAC per subject
figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    cfg                 = [];
    cfg.layout          = 'neuromag306cmb.lay';
    cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};
    cfg.xlim            = [1 18];
    cfg.ylim            = [30 100];
    ft_singleplotTFR(cfg, PAC_all_norm{isubject});
    title(num2str(isubject));
    i = i + 1;   
end

% plot MI over freq per subject
figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    cfg                 = [];
    cfg.layout          = 'neuromag306cmb.lay';
    cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};
%     cfg.xlim            = [1 18];
%     cfg.ylim            = [30 100];
    ft_singleplotER(cfg, MI_all{isubject});
    title(num2str(isubject));
    i = i + 1;   
end

% plot MI topos over freq per subject
figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    cfg                 = [];
    cfg.layout          = 'neuromag306cmb.lay';
%     cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};
    cfg.xlim            = [30 100];
%     cfg.ylim            = [30 100];
    ft_topoplotER(cfg, MI_all{isubject});
    title(num2str(isubject));
    i = i + 1;   
end

% plot PAC difference per subject
figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    cfg                 = [];
    cfg.layout          = 'neuromag306cmb.lay';
    cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};
    cfg.xlim            = [1 18];
    cfg.ylim            = [30 100];
    ft_singleplotTFR(cfg, PAC_diff{isubject});
    title(num2str(isubject));
    i = i + 1;   
end

% plot MI high and low per subject
figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    cfg                 = [];
    cfg.layout          = 'neuromag306cmb.lay';
    cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};
%     cfg.xlim            = [1 18];
%     cfg.ylim            = [30 100];
    ft_singleplotER(cfg, MI_low{isubject}, MI_high{isubject});    
    title(num2str(isubject));
    i = i + 1;   
end

% plot MI difference per subject
figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    cfg                 = [];
    cfg.layout          = 'neuromag306cmb.lay';
    cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};
%     cfg.xlim            = [1 18];
%     cfg.ylim            = [30 100];
    ft_singleplotER(cfg, MI_diff{isubject});    
    title(num2str(isubject));
    i = i + 1;   
end


% plot histogram 
PAC_sel_hist = nan(length(slist),18);
for isubject = slist
    cfg = [];
    cfg.avgoverfreq = 'yes';
    cfg.avgoverchan = 'yes';
    cfg.frequency = [35 60];
    cfg.channel             = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG1312+1313', 'MEG1342+1343', 'MEG2212+2213', 'MEG2222+2223', 'MEG2412+2413'};
    PAC_sel{isubject} = ft_selectdata(cfg,PAC{isubject});
    PAC_sel_hist(isubject,:) = squeeze(PAC_sel{isubject}.powspctrm);
    
    PAC_sel_high{isubject} = ft_selectdata(cfg,PAC_high{isubject});
    PAC_sel_hist_high(isubject,:) = squeeze(PAC_sel_high{isubject}.powspctrm);
    
    PAC_sel_low{isubject} = ft_selectdata(cfg,PAC_low{isubject});
    PAC_sel_hist_low(isubject,:) = squeeze(PAC_sel_low{isubject}.powspctrm);    
    
    PAC_sel_diff{isubject} = ft_selectdata(cfg,PAC_diff{isubject});
    PAC_sel_hist_diff(isubject,:) = squeeze(PAC_sel_diff{isubject}.powspctrm);        
    
end
PAC_sel_hist_avg = nanmean(PAC_sel_hist,1);
PAC_sel_hist_std = nanstd(PAC_sel_hist,1);
PAC_sel_hist_avg_high = nanmean(PAC_sel_hist_high,1);
PAC_sel_hist_std_high = nanstd(PAC_sel_hist_high,1);
PAC_sel_hist_avg_low = nanmean(PAC_sel_hist_low,1);
PAC_sel_hist_std_low = nanstd(PAC_sel_hist_low,1);
PAC_sel_hist_avg_diff = nanmean(PAC_sel_hist_diff,1);
PAC_sel_hist_std_diff = nanstd(PAC_sel_hist_diff,1);
x = linspace(-pi,pi,18);

figure; hold;
errorbar(x,PAC_sel_hist_avg,PAC_sel_hist_std,'r.')
bar(x,PAC_sel_hist_avg); 
axis tight; 
ylim([0.98 1.02]);
legend({'std','power 35-60Hz'});

figure; hold;
errorbar(x,PAC_sel_hist_avg_diff,PAC_sel_hist_std_diff,'r.')
bar(x,PAC_sel_hist_avg_diff); 
axis tight; 
ylim([0.98 1.02]);
legend({'std','power 35-60Hz'});

figure; hold;
errorbar([x; x+0.05]',[PAC_sel_hist_avg_low; PAC_sel_hist_avg_high]',[PAC_sel_hist_std_high; PAC_sel_hist_std_low]','r.');
bar([x; x+0.05]',[PAC_sel_hist_avg_low; PAC_sel_hist_avg_high]',5); 
legend({'std','std','high 35-60Hz','low 35-60Hz'});
axis tight; 
ylim([0.95 1.05]);
hold;




figure; plot([PAC_sel_hist_avg_low; PAC_sel_hist_avg_high]')

% statistics
load('D:\analysis\WANDER\scripts\grad','grad');
neigh_mag = load('neuromag306mag_neighb_last');
neigh_cmb = load('neuromag306cmb_neighb.mat');

Nsub                    = length(slist);
cfg                     = [];
cfg.channel             = 'MEG*3';
    cfg.channel = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG2212+2213', 'MEG2222+2223'};
% cfg.latency             = [9 18];
cfg.frequency           = [40 50];
cfg.avgovertime         = 'no';
cfg.avgoverchan         = 'yes';
cfg.avgoverfreq         = 'yes';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.05;
cfg.numrandomization    = 2000;
cfg.neighbours          = neigh_cmb.neighbours; %  ft_prepare_neighbours(cfg_neighb, grad);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
% 
stat = ft_timelockstatistics(cfg,MI_high{slist},MI_low{slist});

PAC_dummy = PAC;
for isubject = slist
    PAC_dummy{isubject}.powspctrm = ones(size(PAC_dummy{isubject}.powspctrm));
end

stat = ft_freqstatistics(cfg,PAC{slist},PAC_dummy{slist});

    cfg                 = [];
    cfg.layout          = 'neuromag306cmb.lay';
%     cfg.channel             = {'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG1312+1313', 'MEG1342+1343', 'MEG2212+2213', 'MEG2222+2223', 'MEG2412+2413'};
    cfg.xlim            = [1 18];
    cfg.ylim            = [30 100];
%     ft_singleplotTFR(cfg, PAC{isubject});
    ft_singleplotTFR(cfg, PAC_diff{isubject});
    title(num2str(isubject));
    i = i + 1;   
find(stat.negclusterslabelmat == 1)

