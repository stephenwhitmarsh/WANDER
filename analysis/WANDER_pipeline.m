% --- TASK DATA
% subject06_run3 was maxfiltered without xscan, but with autobad:  maxfilter_cenir -f s06run3.fif -st 2000 -corr 0.98 -trans meanhp.fif -autobad on
% subject21_run3 was maxfiltered without xscan, but with autobad, but with sss instead of tsss === VERY BAD, WITH MANY JUMPS THROUGHTOUT RECORDING

% subject05 has noisy EGG sensors, so picked maximum sensor by hand: nr. 3
% subject26 has noisy EGG sensors, so picked maximum sensor by hand: nr. 6
% subject23_run1 has ECG left
% subject05_run1 has an ICA component removed that looks like tactile response
% --- RESTING STATE
% subject11 has an EOG component that includes strong posterior topo
% subject26 has an ECG component that includes a bit of a central lateralized component
% subject03 has noisy MEG data in resting state due to muscle activity
% subject06 has tonic muscle artifacs, esepcially in run 2
% subject06 has clear EOG remaining in run 3, although ICA looks good. Rejected them manually as artefacts
% subject07 has many (frontal) bursts of EMG in all runs, many are manually selected
% subject07 has ECG left in run 3, manually selected as artefacts
% subject08 has some EOG left in run 1
% subject08 has some ECG left in run 2
% subject09 has a lot of highfrequency brain activity (beta or higher)
% rejected an extra EOG ICA component in subject22 resting state, to no avail since still had to select many EOG artefacts by hand


addpath('D:\fieldtrip\fieldtrip.git\trunk');
addpath('D:/analysis/WANDER/scripts/');
ft_defaults

% load filenames
[dataset_task, dataset_rs] = WANDER_subjectinfo;

restingstate = 0;
rootpath     = 1;
isubject     = 1;
force        = 1;

% preprocess EGG (identical for cue or probelocked analysis)
WANDER_FFT_EGG(isubject,0,rootpath,restingstate)
WANDER_EGG_layout(isubject,0)
WANDER_filter_EGG(isubject,force,rootpath,restingstate)
WANDER_filter_EGG_broadband(isubject,force,rootpath,restingstate)
WANDER_artefact_detection_EGG(isubject,1,rootpath,restingstate)

% behaviour
WANDER_behaviour(isubject,1);
WANDER_behaviour_average;
WANDER_phase_ratings;

% preprocess MEG
WANDER_epoch_MEG(isubject,force,rootpath,restingstate)
WANDER_artefact_detection_MEG(isubject,force,rootpath,restingstate)
WANDER_blink_detection(isubject,force,rootpath,restingstate)
WANDER_ICA(isubject,force,rootpath,restingstate)
WANDER_artefact_detection_MEG_after_ICA(isubject,force,rootpath,restingstate)

% cuelocked freqanalysis
WANDER_TFR(isubject,0,'cue',rootpath,restingstate)
WANDER_TFR_corr(isubject,0,'cue')
WANDER_TFR_average(isubject,0,'cue');

% cuelocked timelocked analysis
WANDER_ERF(isubject,0,'cue')

% redefine trials to probelocked
WANDER_redefine_MEG_to_probe(isubject,0)
WANDER_redefine_EGG_to_probe(isubject,0)

% can now do probe-locked analysis
WANDER_TFR(isubject,1,'probe')
WANDER_FFT_MEG(isubject,0,'probe')
WANDER_TFR_phaselocking_balanced(isubject,force,timing,latency,rootpath)
% WANDER_TFR_phaselocking_server_and_rs(isubject,1,'cue','all',rootpath)
WANDER_TFR_corr(isubject,1,'probe');

% averages over subjects
WANDER_TFR_phaselocking_GA(isubject,0,'probe','all')
WANDER_TFR_corr_GA;
WANDER_TFR_median_split_GA(1,'probe');

WANDER_FFT_MEG_alpha_GA(isubject,1,[],1,1);


% source analysis
WANDER_grid(isubject,force)
WANDER_common_filter(isubject,rootpath,force)
WANDER_beamformer(isubject,rootpath,force)
WANDER_source_MI(isubject,rootpath,force)





% LOOPS

restingstate    = 0;
rootpath        = 1;
force           = 0;
latency         = 'all';
timing          = 'cue';

isubject        = 1;
slist = [1:20 22:26];
slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

for isubject = slist
%     try
        % for isubject = [26:-1:14 12:-1:1]
%         WANDER_filter_EGG(isubject,force,rootpath,restingstate);
        
        %     WANDER_epoch_MEG(isubject,1,rootpath,restingstate)
        %     WANDER_artefact_detection_MEG(isubject,1,rootpath,restingstate)
        %     WANDER_blink_detection(isubject,0,rootpath,restingstate)
        %     WANDER_ICA(isubject,1,rootpath,restingstate);
        %     WANDER_FFT_MEG_alpha(isubject,0,[],1,1);
        %     WANDER_TFR_phaselocking_server_and_rs_freqshift(isubject,force,timing,latency,rootpath);
%         WANDER_artefact_detection_MEG_after_ICA(isubject,force,rootpath,restingstate)
        % WANDER_filter_EGG_broadband(isubject,force,rootpath,restingstate)
%         WANDER_artefact_detection_EGG(isubject,1,rootpath,restingstate)
% %          WANDER_TFR_phaselocking_balanced(isubject,force,timing,latency,rootpath)
%             close all
%                 WANDER_ERF(isubject,force,timing,rootpath)
%     catch
%         fprintf('something went wrong with subject %d ',isubject);
%     end
% WANDER_SSPAC(isubject,force,timing,latency,rootpath);
% WANDER_Heart(isubject,force,timing,rootpath);
% WANDER_source_MI(isubject,rootpath,force);

% close all
% WANDER_common_filter(isubject,rootpath,force)
WANDER_grid(isubject,1,1)
end

i = 1;
for isubject = [9:12 14:26]
   
   [FFT_mag{i} FFT_grad{i}] = WANDER_FFT_MEG_alpha(isubject,0,[],1,1);
   
   cfg              = [];
   cfg.frequency    = FFT_mag{i}.max_freq;
   cfg.avgoverfreq  = 'yes';
   FFT_alpha{i}     = ft_selectdata(cfg,FFT_mag{i});
   FFT_grad{i}      = ft_combineplanar([],FFT_grad{i});
   i = i + 1;  
end
% save('trans_FFT','FFT_mag','FFT_grad');

FFT_alpha_GA    = ft_freqgrandaverage([],FFT_alpha{:});
FFT_mag_GA          = ft_freqgrandaverage([],FFT_mag{:});
FFT_grad_GA          = ft_freqgrandaverage([],FFT_grad{:});

% TFR multiplot average
cfg         = [];
cfg.layout  = 'neuromag306all';
cfg.parameter = 'powspctrm';
figure; ft_topoplotER(cfg,FFT_alpha_GA);
figure; ft_topoplotER(cfg,FFT_mag_GA);
figure; ft_topoplotER(cfg,FFT_mag{20});

cfg.layout  = 'neuromag306cmb';
figure; ft_topoplotER(cfg,FFT_grad_GA);





% 
% subplot(2,1,1);
% cfg         = [];
% cfg.layout  = 'neuromag306all';
% cfg.parameter = 'powspctrm';
% cfg.xlim = [FFT_alpha_GA.max_freq FFT_MEG.max_freq];
% cfg.marker = 'off';
% cfg.highlight = 'on';
% cfg.highlightchannel = FFT_MEG.max_chan;
% cfg.highlightsymbol = 'pentagram';
% cfg.highlightsize = 14;
% cfg.comment = 'no';
% ft_topoplotER(cfg,FFT_MEG);
