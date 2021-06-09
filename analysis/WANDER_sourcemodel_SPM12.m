function [mri_segmented_brain, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_sourcemodel_SPM12(isubject,rootpath,force)

if nargin == 0
    disp('Not enough input arguments');
    return;
end

if rootpath == 1
    outputdir_data      = 'W:\WANDER\data\beamformer';
else
    outputdir_data      = '/shared/projects/project_wander/WANDER/data/beamformer';  
end

load([outputdir_data filesep 'template_hdm']);
load([outputdir_data filesep 'template_scm']);
load([outputdir_data filesep 'template_grid']);

fname_grid = [outputdir_data filesep num2str(isubject) '_grid.mat' ];
load(fname_grid);

% segment the MRI - brain
cfg                     = [];
cfg.spmversion          = 'spm12';
cfg.spmmethod           = 'new';
cfg.output              = {'brain','scalp'};
mri_segmented_brain     = ft_volumesegment(cfg, mri_realigned);

% show segmentation
%     figure;
%     mri = imagesc(flipud(squeeze(mri_realigned.anatomy(150,:,:))'));
%     alphaval = zeros(size(mri.CData));
%     hold;
%     brain = image(flipud(squeeze(mri_segmented_brain.gray(150,:,:))')*600);
%     alphaval(find(brain.CData)) = 0.5;
%     set(brain, 'AlphaData', alphaval);
%     axis square
%     saveas(gcf,[outputdir_image filesep num2str(isubject) '_segmentation' ],'png');

% construct the volume conductor model (i.e. head model) for each subject
cfg                     = [];
cfg.method              = 'singleshell';
cfg.spmversion          = 'spm12';
headmodel               = ft_prepare_headmodel(cfg, mri_segmented_brain);

%     if isubject == 8 || isubject == 12
%         mri_realigned.anatomy(mri_segmented_brain.softtissue > 0.2) = 0;
%     end

% create the subject specific grid, using the template grid that has just been created
cfg                     = [];
cfg.spmversion          = 'spm12';
cfg.grid.warpmni        = 'yes';
cfg.grid.template       = template_grid;
cfg.grid.nonlinear      = 'yes';
cfg.mri                 = mri_realigned;
cfg.grid.unit           = 'mm';
subject_grid            = ft_prepare_sourcemodel(cfg);

%     % make a figure of the single subject headmodel, and grid positions
%     figure; hold on;
%     ft_plot_vol(headmodel, 'edgecolor', 'none', 'facecolor', 'brain','facealpha', 0.4);
%     ft_plot_mesh(subject_grid.pos(subject_grid.inside,:));
%     az = 90;
%     el = 0;
%     view(az, el);
%     saveas(gcf,[outputdir_image filesep num2str(isubject) '_headmodelplot' ],'png');

save(fname_grid, 'mri_segmented_brain', 'headmodel', 'subject_grid', 'template_grid','mri_realigned');
