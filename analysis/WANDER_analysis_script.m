% Note that everywhere LOW/HIGH refers to ratings of distraction, so LOW = HIGH ATTENTION, and HIGH = LOW ATTENTION

% addpath /shared/projects/project_wander/WANDER/scripts/
% addpath /shared/projects/project_wander/WANDER/scripts/
% addpath /shared/projects/project_wander/WANDER/scripts/
% addpath /shared/projects/project_wander/WANDER/fieldtrip/
% addpath D:\WANDER\scripts
% addpath D:\WANDER\fieldtrip

if isunix
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/fieldtrip
    addpath /network/lustre/iss01/charpier/analyses/stephen.whitmarsh/WANDER/scripts/
    ft_defaults    
end

dbstop if error

% List of subjects without more than 2 SD of artefacts
slist = [1:5 8:13 15:20 22:26];  

% Get age of subjects
% WANDER_age

% Force recalculation?
force = 0;


% Run scripts from PC (1) or cluster (0)?
rootpath = 0;
timing = 'cue';
restingstate = 0;


for isubject = slist
    
    % Behavioural analysis
    % WANDER_behaviour(isubject,force,rootpath)
    
    % Overview behaviour
    % WANDER_behaviour_average
    
    % Clean up data
    % WANDER_artefact_detection_MEG(isubject,force,rootpath,restingstate)
    % WANDER_ICA(isubject,force,rootpath,restingstate)
    % WANDER_artefact_detection_MEG_after_ICA(isubject,force,rootpath,restingstate)
    
    % Get invidual alpha freq
    % WANDER_FFT_MEG_alpha(isubject,force,timing,rootpath,restingstate)% so
    
    % check potiion of head in helmet
    % WANDER_check_sensor_alignment
    
    % Calculate common filter for individual alpha freq
    % WANDER_common_filter_DICS_individual_onlyalpha(isubject,rootpath,force)
    
    % Source analysis median split
    % WANDER_source_individual_onlyalpha_GRAD(isubject,force,rootpath)
    
    % Source grand average, and create ROI
    % WANDER_source_individual_onlyalpha_GRAD_GA 
    
    % Source expand to ratings based on previous ROI
    % WANDER_source_individual_onlyalpha_GRAD_ratings_GA 

    % Eyelink analysis
    % WANDER_EYE_GA(isubject,force,rootpath);
        
    % SSEF create timelock average of single stim
    % WANDER_timelock_SSEF_noEGG(isubject,rootpath,force) % This needs high-memory cluster

    % SSEF dipole fitting and plotting
    % WANDER_SSEF_GA_topo_and_dipoles
    
    % source data, keeping individual trials
    % WANDER_source_individual_onlyalpha_GRAD_keeptrials(isubject,force,rootpath)
    
    % create GLM for each subject, combining measures
    % WANDER_combine_measures_GLM
    
end



for isubject = slist
    
    [FFT{isubject}, FFT_high{isubject}, FFT_low{isubject}, ~, ~, ~] = WANDER_source_individual_onlyalpha_GRAD(isubject,force,rootpath);
    f(isubject)                    = FFT{isubject}.freq;
    FFT_low{isubject}.freq         = 10;
    FFT_high{isubject}.freq        = 10;
    FFT_low{isubject}.powspctrm    = FFT_low{isubject}.powspctrm';
    FFT_high{isubject}.powspctrm   = FFT_high{isubject}.powspctrm';
    
    
    [dataset_task, ~] = WANDER_subjectinfo(rootpath);
    hdr = ft_read_header(dataset_task{isubject,1});
    FFT_low{isubject}.grad = hdr.grad;
    FFT_high{isubject}.grad = hdr.grad;
% 
%     % for some reason dimord of powscptrm is flipped
%     FFT_low{isubject}.powspctrm = FFT_low{isubject}.powspctrm';
%     FFT_high{isubject}.powspctrm = FFT_high{isubject}.powspctrm';

    FFT_low_comb{isubject}         = ft_combineplanar([], FFT_low{isubject});
    FFT_high_comb{isubject}        = ft_combineplanar([], FFT_high{isubject});
    
    FFT_diff{isubject}             = FFT_low_comb{isubject};
    FFT_diff{isubject}.powspctrm   = (FFT_low_comb{isubject}.powspctrm - FFT_high_comb{isubject}.powspctrm) ./ (FFT_low_comb{isubject}.powspctrm + FFT_high_comb{isubject}.powspctrm);
    
end

cfg                         = [];
FFT_high_GA                 = ft_freqgrandaverage(cfg,FFT_high_comb{slist}); 
FFT_low_GA                  = ft_freqgrandaverage(cfg,FFT_low_comb{slist}); 
FFT_diff_GA                 = ft_freqgrandaverage(cfg,FFT_diff{slist}); 

figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    i=i+1;
    cfg = [];
    cfg.zlim                = 'maxabs';
    cfg.layout              = 'neuromag306cmb';
    cfg.comment             = 'no';
    ft_topoplotER(cfg,FFT_diff{isubject} );
    title(num2str(isubject));
end

figure;
cfg = [];
cfg.zlim                    = 'maxabs';
cfg.layout                  = 'neuromag306cmb';
cfg.comment                 = 'no';
ft_topoplotER(cfg,FFT_diff_GA);


% statistics

load grad 

cfg           = [];
% cfg.channel   = 'MEG*1';
cfg.method    = 'triangulation';
cfg.grad      = grad;
cfg.feedback  = 'yes';
neighbours    = ft_prepare_neighbours(cfg);


cfg = [];
% cfg.channel         = 'MEG*1';
cfg.statistic           = 'depsamplesT';
cfg.ivar                = 1;
cfg.design              = [ones(1, length(slist)) ones(1, length(slist))*2];
cfg.design(2,:)         = [1:length(slist) 1:length(slist)];
cfg.ivar                = 1;
cfg.uvar                = 2;
cfg.method              = 'montecarlo';
cfg.correctm            = 'yes';
cfg.clusteralpha        = 0.025;
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 1;
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.alpha               = 0.025;
cfg.numrandomization    = 4000;
stat                    = ft_freqstatistics(cfg, FFT_low_comb{slist}, FFT_high_comb{slist});

figure;
cfg = [];
cfg.zlim                    = 'maxabs';
cfg.layout                  = 'neuromag306cmb';
cfg.comment                 = 'no';
cfg.highlight               = 'on';
cfg.highlightcolor             = [1 1 1];
cfg.highlightchannel        = find(stat.mask);
cfg.parameter = 'stat';
ft_topoplotER(cfg,stat);

figure;
for i = 1 : length(slist)
    subplot(5,5,i);
    cfg = [];
    cfg.zlim                    = 'maxabs';
    cfg.layout                  = 'neuromag306cmb';
    cfg.comment                 = 'no';
    ft_topoplotER(cfg,FFT_diff{slist(i)});
    title(slist(i));
end

