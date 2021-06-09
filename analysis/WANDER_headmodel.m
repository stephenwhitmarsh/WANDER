function headmodel = WANDER_headmodel(isubject)


if nargin == 0
    disp('Not enough input arguments');
    return;
end

WANDER_subjectinfo

% load subject specific info
subjectdata     = eval(isubject);
shape           = ft_read_headshape(dataset{isubject,1},'unit','cm');
grad            = ft_read_sens(dataset{isubject,1},'senstype','meg');
mri_orig        = ft_read_mri('d:/analysis/WANDER/data/Test08/mm140137/MM140137_MPR.nii');

figure;
subplot(2,2,1);
ft_plot_headshape(shape);
ft_plot_sens(grad, 'style', '*b');
view([1 0 0])
subplot(2,2,2);
ft_plot_headshape(shape);
ft_plot_sens(grad, 'style', '*b');
view([1 1 0])
subplot(2,2,3);
ft_plot_headshape(shape);
ft_plot_sens(grad, 'style', '*b');
view([0 0 1])


figure;
cfg = [];
ft_sourceplot(cfg, mri_orig);


% coregister

cfg = [];
cfg.method          = 'interactive';
cfg.coordsys        = 'neuromag';
[mri_realigned]    = ft_volumerealign(cfg, mri_orig);

mri_real

cfg = [];
cfg.method                      = 'headshape';
cfg.coordsys                    = 'neuromag';
cfg.headshape.headshape         = shape;
cfg.headshape.scalpsmooth       = 4; % (default = 2)
cfg.headshape.scalpthreshold    = 0.08; % (default = 0.1)

[mri_realigned2]                = ft_volumerealign(cfg, mri_realigned);










template_mri    = ft_read_mri(template);

% segment the template brain
cfg                     = [];
% cfg.coordinates         = 'spm';
template_seg            = ft_volumesegment(cfg, template_mri);



try
    load([outputdir filesep 'template_hdm']);
    load([outputdir filesep 'template_scm']);
    disp('template data found...');
catch
    
    disp('template data not found, making it now...');
    
    load grad % just any will do, will not need it further down the line
    
    % segment the template brain
    cfg                     = [];
    cfg.coordinates         = 'spm';
    template_seg            = ft_volumesegment(cfg, template_mri);
    
    % construct a volume conduction model (i.e. head model) from the template
    % only for using it as a boundary for making the grid
    cfg                     = [];
    cfg.method              = 'singleshell';
    template_hdm            = ft_prepare_headmodel(cfg, template_seg);
    
    % make template grid
    cfg                     = [];
    cfg.grid.resolution     = 0.5;
    cfg.inwardshift         = -1.5;
    template_scm            = ft_prepare_sourcemodel(cfg, template_hdm, grad);
    
    % check
    cfg                     = [];
    cfg.vol                 = template_hdm;
    cfg.grid                = template_scm;
    cfg.plotsensors         = 'no';
    ft_headmodelplot(cfg);
    saveas(gcf,[outputdir filesep 'template_headmodelplot' ],'png');
    
    % save template_grid
    save([outputdir filesep 'template_hdm'], 'template_hdm');
    save([outputdir filesep 'template_scm'], 'template_scm');
end