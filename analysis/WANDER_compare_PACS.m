isubject = 11;

temp = load(['W:\WANDER\data\PFR\s' num2str(isubject) '_TFR_phaselocked_V1b.mat']);
FFT_high_V1         = temp.FFT_high;
FFT_low_V1          = temp.FFT_low;
count_orig_V1       = temp.bincount_orig;
count_orig_high_V1  = temp.high_bincount_orig;
count_orig_low_V1   = temp.low_bincount_orig;

cfg = [];
cfg.avgovertime = 'yes';
% FFT_high_V1         = ft_selectdata(cfg,FFT_high_V1);
% FFT_low_V1          = ft_selectdata(cfg,FFT_low_V1);
FFT_high_V1.dimord  = 'chan_freq';
FFT_low_V1.dimord   = 'chan_freq';
FFT_high_V1             = rmfield(FFT_high_V1,'time');
FFT_low_V1              = rmfield(FFT_low_V1,'time');
FFT_high_V1.powspctrm   = squeeze(nanmean(FFT_high_V1.powspctrm,3));
FFT_low_V1.powspctrm    = squeeze(nanmean(FFT_low_V1.powspctrm,3));

PAC_V1              = temp.PPP;
PAC_high_V1         = temp.PPP_high;
PAC_low_V1          = temp.PPP_low;
PAC_V1.dimord       = 'chan_freq_time';
PAC_high_V1.dimord  = 'chan_freq_time';
PAC_low_V1.dimord   = 'chan_freq_time';

temp = load(['W:\WANDER\data\PFR\s' num2str(isubject) '_TFR_phaselocked_V2b.mat']);
FFT_high_V2         = temp.FFT_avg_high;
FFT_low_V2          = temp.FFT_avg_low;

count_orig_V2       = sum(temp.PAC_orig.nrobs,1);
count_orig_high_V2  = sum(temp.PAC_orig_high.nrobs,1);
count_orig_low_V2   = sum(temp.PAC_orig_low.nrobs,1);

% Plot nr of observations per bin
figure;
subplot(1,2,1); hold;
plot(1:18,count_orig_V1); 
plot(1:18,count_orig_high_V1);
plot(1:18,count_orig_low_V1);
axis tight; title('version 1 (TFR)');

subplot(1,2,2); hold;
plot(1:18,count_orig_V2); 
plot(1:18,count_orig_high_V2);
plot(1:18,count_orig_low_V2);
axis tight; title('version 2 (FFT)');


PAC_V2              = temp.PAC_orig;
PAC_high_V2         = temp.PAC_orig_high; 
PAC_low_V2          = temp.PAC_orig_low;
PAC_high_V2.freq    = [1:30];
PAC_low_V2.freq     = [1:30];
PAC_high_V2.dimord  = 'chan_freq_time';
PAC_low_V2.dimord   = 'chan_freq_time';

PAC_V1_norm                 = PAC_V1;
PAC_V1_norm.powspctrm       = PAC_V1_norm.powspctrm ./ nanmean(PAC_V1_norm.powspctrm,3);
PAC_high_V1_norm            = PAC_high_V1;
PAC_high_V1_norm.powspctrm  = PAC_high_V1_norm.powspctrm ./ nanmean(PAC_high_V1_norm.powspctrm,3);
PAC_low_V1_norm             = PAC_low_V1;
PAC_low_V1_norm.powspctrm   = PAC_low_V1_norm.powspctrm  ./ nanmean(PAC_low_V1_norm.powspctrm,3);

PAC_V2_norm                 = PAC_V2;
PAC_V2_norm.powspctrm       = PAC_V2_norm.powspctrm ./ nanmean(PAC_V2_norm.powspctrm,3);
PAC_high_V2_norm            = PAC_high_V2;
PAC_high_V2_norm.powspctrm  = PAC_high_V2_norm.powspctrm ./ nanmean(PAC_high_V2_norm.powspctrm,3);
PAC_low_V2_norm             = PAC_low_V2;
PAC_low_V2_norm.powspctrm   = PAC_low_V2_norm.powspctrm  ./ nanmean(PAC_low_V2_norm.powspctrm,3);

cfg             = [];
cfg.layout      = 'neuromag306cmb';
cfg.channel     = 'MEG*3';
cfg.interactive = 'yes';
cfg.zlim        = 'maxabs';
cfg.channel     = {'MEG2332+2333', 'MEG2342+2343', 'MEG2512+2513'};

cfg.zlim = [-1.2e-23 1.2e-23];
figure;  ft_singleplotTFR(cfg,PAC_high_V1);
figure;  ft_singleplotTFR(cfg,PAC_high_V2);
figure;  ft_singleplotTFR(cfg,PAC_low_V1);
figure;  ft_singleplotTFR(cfg,PAC_low_V2);

cfg.zlim = [0.8 1.2];
figure;  ft_singleplotTFR(cfg,PAC_high_V1_norm);
figure;  ft_singleplotTFR(cfg,PAC_low_V1_norm);
figure;  ft_singleplotTFR(cfg,PAC_high_V2_norm);
figure;  ft_singleplotTFR(cfg,PAC_low_V2_norm);

cfg.layout      = 'neuromag306mag';
cfg.channel     = 'MEG*1';

figure; ft_singleplotER(cfg,FFT_high_V1,FFT_high_V2);
figure; ft_topoplotER(cfg,FFT_high_V1,FFT_high_V2)
figure; ft_singleplotER(cfg,FFT_low_V1,FFT_low_V2);

% PAC per channel version 1
figure;
for ichan = 1 : 100
    subplot(10,10,ichan); hold;
    plot(squeeze(PAC_high_V1_norm.powspctrm(ichan,10,:)));
    plot(squeeze(PAC_high_V1_norm.powspctrm(ichan,11,:)));
    ylim([.5 1.5]);
end

% PAC per channel version 1
figure;
for ichan = 1 : 100
    subplot(10,10,ichan); hold;
    plot(squeeze(PAC_high_V2_norm.powspctrm(ichan,10,:)));
    plot(squeeze(PAC_high_V2_norm.powspctrm(ichan,11,:)));
    ylim([.5 1.5]);
end

% PAC compared 
figure;
for ichan = 1 : 100
    subplot(10,10,ichan); hold;
    plot(squeeze(PAC_high_V1_norm.powspctrm(ichan,10,:)));
    plot(squeeze(PAC_high_V2_norm.powspctrm(ichan,10,:)));
    ylim([.5 1.5]);
end

figure;  hold;
plot(squeeze(PAC_low_V1_norm.powspctrm(ichan,10,:))); 
plot(squeeze(PAC_low_V1_norm.powspctrm(ichan,11,:))); 
ylim([.5 1.5]);

figure;  hold;
plot(squeeze(PAC_low_V2_norm.powspctrm(ichan,10,:))); 
plot(squeeze(PAC_low_V2_norm.powspctrm(ichan,11,:))); 
ylim([.5 1.5]);

figure;  hold;
plot(squeeze(PAC_V1_norm.powspctrm(1,10,:)));
plot(squeeze(PAC_V2_norm.powspctrm(1,10,:)));



%%%
