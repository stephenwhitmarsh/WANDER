slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

template_mri            = ft_read_mri('D:\fieldtrip\fieldtrip.git\trunk\external\spm8\templates\T1.nii');
template_mri.coordsys   = 'spm';

cfg                     = [];
cfg.output              = 'brain';
cfg.template            = 'D:\fieldtrip\fieldtrip.git\trunk\external\spm8\templates\T1.nii';
template_seg            = ft_volumesegment(cfg, template_mri);

[mri_segmented, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(1,1,0);
template_mri.inside     = template_seg.brain;

for isubject = slist
  
    [source_MI_matrix_all{isubject}, source_MI_matrix_high{isubject}, source_MI_matrix_low{isubject}, source_MI_matrix_diff{isubject}, source_TFR_all{isubject}, source_TFR_high{isubject}, source_TFR_low{isubject}, source_TFR_diff{isubject}] = WANDER_source_MI(isubject,1,0);
    
    source_MI_matrix_all{isubject}.pos          = template_grid.pos;
    source_MI_matrix_all{isubject}.pos          = source_MI_matrix_all{isubject}.pos(template_grid.inside,:);
    source_MI_matrix_diff{isubject}.pos         = template_grid.pos;
    source_MI_matrix_diff{isubject}.pos         = source_MI_matrix_diff{isubject}.pos(template_grid.inside,:);
    source_MI_matrix_high{isubject}.pos         = template_grid.pos;
    source_MI_matrix_high{isubject}.pos         = source_MI_matrix_high{isubject}.pos(template_grid.inside,:);
    source_MI_matrix_low{isubject}.pos          = template_grid.pos;
    source_MI_matrix_low{isubject}.pos          = source_MI_matrix_low{isubject}.pos(template_grid.inside,:);
    
    source_TFR_all{isubject}.pos                = template_grid.pos;
    source_TFR_all{isubject}.pos                = source_TFR_all{isubject}.pos(template_grid.inside,:);
    source_TFR_diff{isubject}.pos               = template_grid.pos;
    source_TFR_diff{isubject}.pos               = source_TFR_diff{isubject}.pos(template_grid.inside,:);
    source_TFR_high{isubject}.pos               = template_grid.pos;
    source_TFR_high{isubject}.pos               = source_TFR_high{isubject}.pos(template_grid.inside,:);
    source_TFR_low{isubject}.pos                = template_grid.pos;
    source_TFR_low{isubject}.pos                = source_TFR_low{isubject}.pos(template_grid.inside,:);
       
    source_TFR_high_norm{isubject}              = source_TFR_high{isubject};
    source_TFR_high_norm{isubject}.pow          = source_TFR_high{isubject}.pow ./ (source_TFR_low{isubject}.pow + source_TFR_high{isubject}.pow);
    source_TFR_low_norm{isubject}               = source_TFR_low{isubject};
    source_TFR_low_norm{isubject}.pow           = source_TFR_low{isubject}.pow  ./ (source_TFR_low{isubject}.pow + source_TFR_high{isubject}.pow);
           
    cfg = []; % MAYBE TRY TO AVERAGE MI, I.E. LATER
    cfg.avgoverfreq = 'yes';
    source_MI_matrix_all{isubject}              = ft_selectdata(cfg,source_MI_matrix_all{isubject});
    source_MI_matrix_diff{isubject}             = ft_selectdata(cfg,source_MI_matrix_diff{isubject});
    source_MI_matrix_high{isubject}             = ft_selectdata(cfg,source_MI_matrix_high{isubject});
    source_MI_matrix_low{isubject}              = ft_selectdata(cfg,source_MI_matrix_low{isubject});
    
    source_TFR_all{isubject}                    = ft_selectdata(cfg,source_TFR_all{isubject});
    source_TFR_diff{isubject}                   = ft_selectdata(cfg,source_TFR_diff{isubject});
    source_TFR_high{isubject}                   = ft_selectdata(cfg,source_TFR_high{isubject});
    source_TFR_low{isubject}                    = ft_selectdata(cfg,source_TFR_low{isubject});   
    
    source_TFR_all{isubject}.powdimord          = 'pos';
    source_TFR_diff{isubject}.powdimord         = 'pos';
    source_TFR_high{isubject}.powdimord         = 'pos';
    source_TFR_low{isubject}.powdimord          = 'pos';
    
    source_MI_matrix_all{isubject}.MIdimord     = 'pos';
    source_MI_matrix_diff{isubject}.MIdimord    = 'pos';
    source_MI_matrix_high{isubject}.MIdimord    = 'pos';
    source_MI_matrix_low{isubject}.MIdimord     = 'pos';
    
    source_MI_matrix_all{isubject}.dim          = template_grid.dim;
    source_MI_matrix_diff{isubject}.dim         = template_grid.dim;
    source_MI_matrix_high{isubject}.dim         = template_grid.dim;
    source_MI_matrix_low{isubject}.dim          = template_grid.dim;
    
    source_TFR_all{isubject}.dim                = template_grid.dim;
    source_TFR_diff{isubject}.dim               = template_grid.dim;
    source_TFR_high{isubject}.dim               = template_grid.dim;
    source_TFR_low{isubject}.dim                = template_grid.dim;    
    
    source_MI_matrix_high_norm{isubject}        = source_MI_matrix_high{isubject};
    source_MI_matrix_high_norm{isubject}.MI     = source_MI_matrix_high{isubject}.MI ./ (source_MI_matrix_low{isubject}.MI + source_MI_matrix_high{isubject}.MI);
    source_MI_matrix_low_norm{isubject}         = source_MI_matrix_low{isubject};
    source_MI_matrix_low_norm{isubject}.MI      = source_MI_matrix_low{isubject}.MI ./ (source_MI_matrix_low{isubject}.MI + source_MI_matrix_high{isubject}.MI);  
    
    source_TFR_high_norm{isubject}              = source_TFR_high{isubject};
    source_TFR_high_norm{isubject}.pow          = source_TFR_high{isubject}.pow ./ (source_TFR_low{isubject}.pow + source_TFR_high{isubject}.pow);
    source_TFR_low_norm{isubject}               = source_TFR_low{isubject};
    source_TFR_low_norm{isubject}.pow           = source_TFR_low{isubject}.pow ./ (source_TFR_low{isubject}.pow + source_TFR_high{isubject}.pow);  

end

cfg = [];
cfg.parameter               = 'pow';
source_TFR_diff_avg         = ft_sourcegrandaverage(cfg,source_TFR_diff{slist});
source_TFR_all_avg          = ft_sourcegrandaverage(cfg,source_TFR_all{slist});
cfg.parameter               = 'MI';
source_MI_matrix_diff_avg   = ft_sourcegrandaverage(cfg,source_MI_matrix_diff{slist});

% sparse to full represetnation
temp = source_TFR_all_avg.pow;
source_TFR_all_avg.statdimord                                = 'pos';
source_TFR_all_avg.pos                                       = template_grid.pos;
source_TFR_all_avg.pow                                       = zeros(prod(source_TFR_all_avg.dim),1);
source_TFR_all_avg.pow(template_grid.inside)                 = temp;

temp = source_TFR_diff_avg.pow;
source_TFR_diff_avg.statdimord                                = 'pos';
source_TFR_diff_avg.pos                                       = template_grid.pos;
source_TFR_diff_avg.pow                                       = zeros(prod(source_TFR_diff_avg.dim),1);
source_TFR_diff_avg.pow(template_grid.inside)                 = temp;

temp = source_MI_matrix_diff_avg.MI;
source_MI_matrix_diff_avg.statdimord                                = 'pos';
source_MI_matrix_diff_avg.pos                                       = template_grid.pos;
source_MI_matrix_diff_avg.MI                                       = zeros(prod(source_MI_matrix_diff_avg.dim),1);
source_MI_matrix_diff_avg.MI(template_grid.inside)                 = temp;

cfg = [];
cfg.parameter                       = {'pow'};
source_TFR_all_avg                  = ft_sourceinterpolate(cfg, source_TFR_all_avg, template_mri);
source_TFR_all_avg.inside           = template_mri.inside;
source_TFR_diff_avg                 = ft_sourceinterpolate(cfg, source_TFR_diff_avg, template_mri);
source_TFR_diff_avg.inside          = template_mri.inside;
cfg.parameter                       = {'MI'};
source_MI_matrix_diff_avg           = ft_sourceinterpolate(cfg, source_MI_matrix_diff_avg, template_mri);
source_MI_matrix_diff_avg.inside    = template_mri.inside;


% plot averages
cfg = [];
cfg.method        = 'ortho';
cfg.funparameter  = 'pow';
ft_sourceplot(cfg, source_TFR_all_avg);
ft_sourceplot(cfg, source_TFR_diff_avg);
cfg.funparameter  = 'MI';
ft_sourceplot(cfg, source_MI_matrix_diff_avg);


cfg.method         = 'surface';
cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
cfg.surfinflated   = 'surface_inflated_both.mat'; % Cortical sheet from canonical MNI brain
cfg.funparameter  = 'pow';
ft_sourceplot(cfg, source_TFR_all_avg);
ft_sourceplot(cfg, source_TFR_diff_avg);
material dull


cfg.funparameter  = 'MI';
ft_sourceplot(cfg, source_MI_matrix_diff_avg);





% cluster corrected analysis
cfg = [];
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 10000;
cfg.alpha               = 0.05;
cfg.clusteralpha        = 0.01; 
cfg.tail                = 0;
cfg.design(1,:)         = [1:length(slist) 1:length(slist)];
cfg.design(2,:)         = [ones(size(slist)) ones(size(slist))*2];
cfg.uvar                = 1; 
cfg.ivar                = 2;
cfg.parameter           = 'MI';
stat_MI                 = ft_sourcestatistics(cfg, source_MI_matrix_low_norm{slist}, source_MI_matrix_high_norm{slist});

cfg.parameter           = 'pow';
stat_TFR                = ft_sourcestatistics(cfg, source_TFR_low_norm{slist}, source_TFR_high_norm{slist});


% convert to full to allow interpolation
stat_MI_full                                            = keepfields(stat_MI,{'label','stat'});
stat_MI_full.statdimord                                 = 'pos';
stat_MI_full.dimord                                     = 'pos';
stat_MI_full.dim                                        = template_grid.dim;
stat_MI_full.pos                                        = template_grid.pos;
stat_MI_full.stat                                       = zeros(prod(stat_MI_full.dim),1);
stat_MI_full.stat(template_grid.inside)                 = stat_MI.stat;
stat_MI_full.mask                                       = zeros(prod(stat_MI_full.dim),1);
stat_MI_full.mask(template_grid.inside)                 = stat_MI.mask;
stat_MI_full.negclusterslabelmat                        = zeros(prod(stat_MI_full.dim),1);
stat_MI_full.negclusterslabelmat(template_grid.inside)  = stat_MI.negclusterslabelmat;

cfg = [];
cfg.parameter       = {'stat','mask'};
stat_MI_int         = ft_sourceinterpolate(cfg, stat_MI_full, template_mri);
stat_MI_int.inside  = template_mri.inside;


% convert to full to allow interpolation
clear stat_TFR_full
stat_TFR_full                                            = keepfields(stat_TFR,{'label','stat'});
stat_TFR_full.statdimord                                 = 'pos';
stat_TFR_full.dimord                                     = 'pos';
stat_TFR_full.dim                                        = template_grid.dim;
stat_TFR_full.pos                                        = template_grid.pos;
stat_TFR_full.stat                                       = zeros(prod(stat_TFR_full.dim),1);
stat_TFR_full.stat(template_grid.inside)                 = stat_TFR.stat;
stat_TFR_full.negclusterslabelmat                        = zeros(prod(stat_TFR_full.dim),1);
stat_TFR_full.mask                                       = zeros(prod(stat_MI_full.dim),1);
stat_TFR_full.mask(template_grid.inside)                 = stat_TFR.mask;
stat_TFR_full.negclusterslabelmat                        = zeros(prod(stat_TFR_full.dim),1);
stat_TFR_full.negclusterslabelmat(template_grid.inside)  = stat_TFR.negclusterslabelmat;

cfg = [];
cfg.parameter       = {'stat','mask','negclusterslabelmat'};
stat_TFR_int         = ft_sourceinterpolate(cfg, stat_TFR_full, template_mri);
stat_TFR_int.inside  = template_mri.inside;


cfg = [];
cfg.method        = 'ortho';
cfg.funcolorlim   = [-5 5];
cfg.funparameter  = 'stat';
% cfg.maskparameter = 'mask';
% cfg.funcolorlim   = [0.0 1.2];
% cfg.opacitylim    = [0 1];
cfg.opacitymap    = 'rampup';
cfg.maskparameter = 'mask';

ft_sourceplot(cfg, stat_TFR_int);
ft_sourceplot(cfg, stat_MI_int);


cfg.method         = 'surface';
cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
cfg.surfinflated   = 'surface_inflated_both.mat'; % Cortical sheet from canonical MNI brain

ft_sourceplot(cfg, stat_TFR_int);
ft_sourceplot(cfg, stat_MI_int);

material dull



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% cluster uncorrected %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cluster UNcorrected analysis
cfg = [];
% cfg.dim         = 3;
cfg.method              = 'analytic';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.alpha               = 0.05; % note that this only implies single-sided testing
cfg.tail                = 0;
cfg.design(1,:)         = [1:length(slist) 1:length(slist)];
cfg.design(2,:)         = [ones(size(slist)) ones(size(slist))*2];
cfg.uvar                = 1; % row of design matrix that contains unit variable (in this case: trials)
cfg.ivar                = 2; % row of design matrix that contains independent variable (the conditions)
cfg.parameter           = 'MI';
stat_MI                 = ft_sourcestatistics(cfg, source_MI_matrix_low_norm{slist}, source_MI_matrix_high_norm{slist});

cfg.parameter           = 'pow';
stat_TFR                = ft_sourcestatistics(cfg, source_TFR_low_norm{slist}, source_TFR_high_norm{slist});


% convert to full to allow interpolation
stat_MI_full                                            = keepfields(stat_MI,{'label','stat'});
stat_MI_full.statdimord                                 = 'pos';
stat_MI_full.dimord                                     = 'pos';
stat_MI_full.dim                                        = template_grid.dim;
stat_MI_full.pos                                        = template_grid.pos;
stat_MI_full.stat                                       = zeros(prod(stat_MI_full.dim),1);
stat_MI_full.stat(template_grid.inside)                 = stat_MI.stat;
stat_MI_full.mask                                       = zeros(prod(stat_MI_full.dim),1);
stat_MI_full.mask(template_grid.inside)                 = stat_MI.mask;

cfg = [];
cfg.parameter       = {'stat','mask'};
stat_MI_int         = ft_sourceinterpolate(cfg, stat_MI_full, template_mri);
stat_MI_int.inside  = template_mri.inside;


% convert to full to allow interpolation
clear stat_TFR_full
stat_TFR_full                                            = keepfields(stat_TFR,{'label','stat'});
stat_TFR_full.statdimord                                 = 'pos';
stat_TFR_full.dimord                                     = 'pos';
stat_TFR_full.dim                                        = template_grid.dim;
stat_TFR_full.pos                                        = template_grid.pos;
stat_TFR_full.stat                                       = zeros(prod(stat_TFR_full.dim),1);
stat_TFR_full.stat(template_grid.inside)                 = stat_TFR.stat;
stat_TFR_full.mask                                       = zeros(prod(stat_MI_full.dim),1);
stat_TFR_full.mask(template_grid.inside)                 = stat_TFR.mask;

cfg = [];
cfg.parameter       = {'stat','mask'};
stat_TFR_int         = ft_sourceinterpolate(cfg, stat_TFR_full, template_mri);
stat_TFR_int.inside  = template_mri.inside;


cfg = [];
cfg.method        = 'ortho';
cfg.funcolorlim   = [-5 5];
cfg.funparameter  = 'stat';
% cfg.maskparameter = 'mask';
% cfg.funcolorlim   = [0.0 1.2];
cfg.opacitylim    = [0 1];
cfg.opacitymap    = 'rampup';
cfg.maskparameter = 'mask';

ft_sourceplot(cfg, stat_MI_int);

cfg.method         = 'surface';
cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
cfg.surfinflated   = 'surface_inflated_both.mat'; % Cortical sheet from canonical MNI brain

ft_sourceplot(cfg, stat_MI_int);

