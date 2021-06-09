function [mri_segmented_brain, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(isubject,rootpath,force)

if nargin == 0
    disp('Not enough input arguments');
    return;
end

if rootpath == 1
    datapath = 'D:';
else
    datapath = '/shared/projects/project_wander';
end

% slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD
% addpath('D:\fieldtrip-master_21092017\fieldtrip-master');
% ft_defaults

% if rootpath == 1
template_file       = [datapath,filesep,'WANDER',filesep,'fieldtrip',filesep,'template',filesep,'anatomy',filesep,'single_subj_T1.nii'];
outputdir_data      = [datapath,filesep,'WANDER',filesep,'data',filesep,'beamformer'];
outputdir_image     = [datapath,filesep,'WANDER',filesep,'images',filesep,'beamformer'];
fname_MRI{1}        = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test01',filesep,'2016_05_09_WANDER_16_05_01_04',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160509160503378.MR.dic'];
fname_MRI{2}        = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test02',filesep,'Sujet90',filesep,'Anat',filesep,'s_S11_t1mpr_SAG_NSel_S176.nii'];
fname_MRI{3}        = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test03',filesep,'2015_02_27_THRESHOLD_S06',filesep,'S11_t1mpr_SAG_NSel_S224_0_8iso',filesep,'20150227164905829.MR.dic'];
fname_MRI{4}        = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test04',filesep,'S04_3DT1',filesep,'20140116173427161.MR.dic'];
fname_MRI{5}        = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test05',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20151027135859918.MR.dic'];
fname_MRI{6}        = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test06',filesep,'S08_mp2rage_1iso_ipat2_T1_Images',filesep,'20150129154324466.MR.dic']; % NO T1, but not used anyway
fname_MRI{7}        = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test07',filesep,'2016_05_13_WANDER_16_05_01_07',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160513125117013.MR.dic'];
fname_MRI{8}        = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test08',filesep,'mm140137',filesep,'MM140137_MPR_CROPPED.nii']; % grid too large on top, probably due to wrong scalp segmentation and therefor normalization
fname_MRI{9}        = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test09',filesep,'S04_3DT1',filesep,'20140115105944084.MR.dic']; % grid a bit too small, especially frontally
fname_MRI{10}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test10',filesep,'Anat',filesep,'s_S11_t1mpr_SAG_NSel_S176.nii'];
fname_MRI{11}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test11',filesep,'2016_06_02_WANDER_16_05_01_11',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160602112315979.MR.dic'];
fname_MRI{12}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test12',filesep,'15-05-01-40',filesep,'S04_mp2rage_1iso_ipat2_INV1',filesep,'20151124121653179.MR.dic']; % grid seems too large on top
fname_MRI{13}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test13',filesep,'2016_05_24_WANDER_16_05_01_13_E2',filesep,'S03_t1_mpr_sag_0_8iso',filesep,'20160524132446395.MR.dic'];
fname_MRI{14}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test14',filesep,'S04_3DT1',filesep,'20140120101403044.MR.dic'];
fname_MRI{15}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test15',filesep,'2016_06_06_WANDER_16_05_01_15',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160606173031429.MR.dic'];
fname_MRI{16}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test16',filesep,'15-05-01-14',filesep,'S04_mp2rage_1iso_ipat2_INV1',filesep,'20151117105942397.MR.dic']; % grid is not well-fitting
fname_MRI{17}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test17',filesep,'2016_06_03_WANDER_16_05_01_17',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160603125354990.MR.dic'];
fname_MRI{18}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test18',filesep,'2016_06_03_WANDER_16_05_01_18',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160603172447161.MR.dic'];
fname_MRI{19}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test19',filesep,'2016_06_10_WANDER_16_05_01_19',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160610171054717.MR.dic'];
fname_MRI{20}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test20',filesep,'2016_06_10_WANDER_16_05_01_20',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160610173819059.MR.dic'];
fname_MRI{21}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test21',filesep,'2016_06_10_WANDER_16_05_01_70',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160610124230474.MR.dic'];
fname_MRI{22}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test22',filesep,'2016_06_23_WANDER_16_05_01_34',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160623130006747.MR.dic'];
fname_MRI{23}       = template_file;
fname_MRI{24}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test24',filesep,'2016_06_22_WANDER_16_05_01_28',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160622135330623.MR.dic'];
fname_MRI{25}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test25',filesep,'2016_07_18_WANDER_16_05_01_61',filesep,'S02_t1_mpr_sag_0_8iso',filesep,'20160718102304814.MR.dic'];
fname_MRI{26}       = [datapath,filesep,'WANDER',filesep,'data',filesep,'MRI',filesep,'Test26',filesep,'S02_t1mpr_SAG_NSel_S176',filesep,'20150624164211176.MR.dic'];
% else
%
%     template_file       = '/shared/projects/project_wander/WANDER/fieldtrip\external\spm8\templates\T1.nii';
%     outputdir_data      = '/shared/projects/project_wander/WANDER/data/beamformer/';
%     outputdir_image     = '/shared/projects/project_wander/WANDER/data/beamformer/';
%
%     fname_MRI{1}        = '/shared/projects/project_wander/WANDER/data/MRI\Test01\2016_05_09_WANDER_16_05_01_04\S02_t1_mpr_sag_0_8iso\20160509160503378.MR.dic';
%     fname_MRI{2}        = '/shared/projects/project_wander/WANDER/data/MRI\Test02\Sujet90\Anat\s_S11_t1mpr_SAG_NSel_S176.nii';
%     fname_MRI{3}        = '/shared/projects/project_wander/WANDER/data/Test03\2015_02_27_THRESHOLD_S06\S11_t1mpr_SAG_NSel_S224_0_8iso\20150227164905829.MR.dic';
%     fname_MRI{4}        = '/shared/projects/project_wander/WANDER/data/Test04\S04_3DT1\20140116173427161.MR.dic';
%     % fname_MRI{4}        = 'W:\WANDER\fieldtrip\template\anatomy\single_subj_T1_1mm.nii';
%     fname_MRI{5}        = '/shared/projects/project_wander/WANDER/data/MRI\Test05\S02_t1_mpr_sag_0_8iso\20151027135859918.MR.dic';
%     fname_MRI{6}        = '/shared/projects/project_wander/WANDER/data/MRI\Test06\S08_mp2rage_1iso_ipat2_T1_Images\20150129154324466.MR.dic'; % NO T1, but not used anyway
%     fname_MRI{7}        = '/shared/projects/project_wander/WANDER/data/MRI\Test07\2016_05_13_WANDER_16_05_01_07\S02_t1_mpr_sag_0_8iso\20160513125117013.MR.dic';
%     fname_MRI{8}        = '/shared/projects/project_wander/WANDER/data/MRI\Test08\mm140137\MM140137_MPR_CROPPED.nii'; % grid too large on top, probably due to wrong scalp segmentation and therefor normalization
%     fname_MRI{9}        = '/shared/projects/project_wander/WANDER/data/MRI\Test09\S04_3DT1\20140115105944084.MR.dic'; % grid a bit too small, especially frontally
%     fname_MRI{10}       = '/shared/projects/project_wander/WANDER/data/MRI\Test10\Anat\s_S11_t1mpr_SAG_NSel_S176.nii';
%     fname_MRI{11}       = '/shared/projects/project_wander/WANDER/data/MRI\Test11\2016_06_02_WANDER_16_05_01_11\S02_t1_mpr_sag_0_8iso\20160602112315979.MR.dic';
%     fname_MRI{12}       = '/shared/projects/project_wander/WANDER/data/MRI\Test12\15-05-01-40\S04_mp2rage_1iso_ipat2_INV1\20151124121653179.MR.dic'; % grid seems too large on top
%     fname_MRI{13}       = '/shared/projects/project_wander/WANDER/data/MRI\Test13\2016_05_24_WANDER_16_05_01_13_E2\S03_t1_mpr_sag_0_8iso\20160524132446395.MR.dic';
%     fname_MRI{14}       = '/shared/projects/project_wander/WANDER/data/MRI\Test14\S04_3DT1\20140120101403044.MR.dic';
%     fname_MRI{15}       = '/shared/projects/project_wander/WANDER/data/MRI\Test15\2016_06_06_WANDER_16_05_01_15\S02_t1_mpr_sag_0_8iso\20160606173031429.MR.dic';
%     fname_MRI{16}       = '/shared/projects/project_wander/WANDER/data/MRI\Test16\15-05-01-14\S04_mp2rage_1iso_ipat2_INV1\20151117105942397.MR.dic'; % grid is not well-fitting
%     fname_MRI{17}       = '/shared/projects/project_wander/WANDER/data/MRI\Test17\2016_06_03_WANDER_16_05_01_17\S02_t1_mpr_sag_0_8iso\20160603125354990.MR.dic';
%     fname_MRI{18}       = '/shared/projects/project_wander/WANDER/data/MRI\Test18\2016_06_03_WANDER_16_05_01_18\S02_t1_mpr_sag_0_8iso\20160603172447161.MR.dic';
%     fname_MRI{19}       = '/shared/projects/project_wander/WANDER/data/MRI\Test19\2016_06_10_WANDER_16_05_01_19\S02_t1_mpr_sag_0_8iso\20160610171054717.MR.dic';
%     fname_MRI{20}       = '/shared/projects/project_wander/WANDER/data/MRI\Test20\2016_06_10_WANDER_16_05_01_20\S02_t1_mpr_sag_0_8iso\20160610173819059.MR.dic';
%     fname_MRI{21}       = '/shared/projects/project_wander/WANDER/data/MRI\Test21\2016_06_10_WANDER_16_05_01_70\S02_t1_mpr_sag_0_8iso\20160610124230474.MR.dic';
%     fname_MRI{22}       = '/shared/projects/project_wander/WANDER/data/MRI\Test22\2016_06_23_WANDER_16_05_01_34\S02_t1_mpr_sag_0_8iso\20160623130006747.MR.dic';
%     fname_MRI{23}       = 'D:\fieldtrip\fieldtrip.git\trunk\template\anatomy\single_subj_T1_1mm.nii';
%     fname_MRI{24}       = '/shared/projects/project_wander/WANDER/data/MRI\Test24\2016_06_22_WANDER_16_05_01_28\S02_t1_mpr_sag_0_8iso\20160622135330623.MR.dic';
%     fname_MRI{25}       = '/shared/projects/project_wander/WANDER/data/MRI\Test25\2016_07_18_WANDER_16_05_01_61\S02_t1_mpr_sag_0_8iso\20160718102304814.MR.dic';
%     fname_MRI{26}       = '/shared/projects/project_wander/WANDER/data/MRI\Test26\S02_t1mpr_SAG_NSel_S176\20150624164211176.MR.dic';
%
% end


fname_grid          = [outputdir_data filesep num2str(isubject) '_grid.mat' ];
[dataset_task, ~]   = WANDER_subjectinfo(rootpath);

try
    load([outputdir_data filesep 'template_hdm']);
    load([outputdir_data filesep 'template_scm']);
    load([outputdir_data filesep 'template_grid']);
    disp('template headmodel found...');
catch
    
    % read template data
    disp('template data not found, making it now...');
    template_mri            = ft_read_mri(template_file);
    template_mri.coordsys   = 'spm';
    
    % segment the template brain
    cfg                     = [];
    cfg.template            = template_file;
    template_seg            = ft_volumesegment(cfg, template_mri);
    
    % construct a volume conduction model (i.e. head model) from the template
    % only for using it as a boundary for making the grid
    cfg                     = [];
    cfg.method              = 'singleshell';
    template_hdm            = ft_prepare_headmodel(cfg, template_seg);
    
    % make template grid
    cfg                     = [];
    cfg.grid.resolution     = 5;
    cfg.inwardshift         = -20;
    template_scm            = ft_prepare_sourcemodel(cfg, template_hdm);
    
    % check
    figure
    hold on
    ft_plot_vol(template_hdm, 'facecolor', 'cortex', 'edgecolor', 'none'); alpha 0.5; camlight;
    ft_plot_mesh(template_scm.pos(template_scm.inside,:));
    
    % construct the dipole grid in the template brain coordinates
    % the negative inwardshift means an outward shift of the brain surface for inside/outside detection
    cfg = [];
    cfg.grid.resolution     = 5; % used to be 10
    cfg.grid.tight          = 'yes';
    cfg.inwardshift         = -20; % used to be 15
    cfg.headmodel           = template_hdm;
    template_grid           = ft_prepare_sourcemodel(cfg);
    
    % save template_grid
    save([outputdir_data filesep 'template_hdm'], 'template_hdm');
    save([outputdir_data filesep 'template_scm'], 'template_scm');
    save([outputdir_data filesep 'template_grid'], 'template_grid');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
% creating subject grids %
%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(fname_grid,'file') && force ~= 1
    
    load(fname_grid);
    disp('subject headmodel found...');
    
else
    
    mri_orig                = ft_read_mri(fname_MRI{isubject});
    headshape               = ft_read_headshape(dataset_task{isubject,1},'unit','mm');
    
    % reslice
    mri_orig                = ft_volumereslice([],mri_orig);
    
    % coregister - reslicing took care of radiological mapping
    cfg = [];
    cfg.method              = 'interactive';
    cfg.coordsys            = 'neuromag';
    cfg.spmversion          = 'spm12';
    [mri_realigned]         = ft_volumerealign(cfg, mri_orig);
    
    cfg = [];
    cfg.method              = 'headshape';
    cfg.coordsys            = 'neuromag';
    cfg.headshape.icp       = 'yes';
    cfg.spmversion          = 'spm12';
    cfg.headshape.headshape = headshape;
    [mri_realigned]         = ft_volumerealign(cfg, mri_realigned);
    
    cfg = [];
    cfg.spmversion          = 'spm12';
    cfg.method              = 'headshape';
    cfg.coordsys            = 'neuromag';
    cfg.headshape.icp       = 'yes';
    cfg.headshape.headshape = headshape;
    [mri_realigned]         = ft_volumerealign(cfg, mri_realigned);
    
    % segment the MRI - brain
    cfg                     = [];
    cfg.spmversion          = 'spm12';
    cfg.spmmethod           = 'new';
    %     cfg.output              = 'brain';
    cfg.output              = {'brain','skull','scalp'};
    %     if isubject == 8
    %         cfg.brainthreshold      = 0.7; % higher is more strict
    %     end
    %     if isubject == 2 || isubject == 4 || isubject == 9
    %         cfg.brainthreshold      = 0.9; % higher is more strict
    %     end
    %     if isubject == 10 || isubject == 11
    %         cfg.brainthreshold      = 0.99; % higher is more strict
    %     end
    mri_segmented_brain     = ft_volumesegment(cfg, mri_realigned);
    
    % show segmentation
    figure;
    mri = imagesc(flipud(squeeze(mri_realigned.anatomy(150,:,:))'));
    alphaval = zeros(size(mri.CData));
    hold;
    brain = image(flipud(squeeze(mri_segmented_brain.brain(150,:,:))')*600);
    alphaval(find(brain.CData)) = 0.5;
    set(brain, 'AlphaData', alphaval);
    axis square
    % %
    %     figure;
    %     mri = imagesc(flipud(squeeze(mri_segmented_brain.air(150,:,:)>0.1)'));
    %     alphaval = zeros(size(mri.CData));
    %     axis square
    
    %
    %
    %     cfg = [];
    %     cfg.funparameter = 'softtissue';
    %     ft_sourceplot(cfg,mri_segmented_brain)
    %
    %
    % construct the volume conductor model (i.e. head model) for each subject
    cfg                     = [];
    cfg.method              = 'singleshell';
    cfg.spmversion          = 'spm12';
    headmodel               = ft_prepare_headmodel(cfg, mri_segmented_brain);
    
    %     if isubject == 8 || isubject == 12
    %         mri_realigned.anatomy(mri_segmented_brain.softtissue > 0.2) = 0;
    %     end
    %
    
    % create the subject specific grid, using the template grid that has just been created
    cfg                     = [];
    cfg.spmversion          = 'spm12';
    cfg.grid.warpmni        = 'yes';
    cfg.grid.template       = template_grid;
    cfg.grid.nonlinear      = 'yes';
    %     if isubject == 12
    %         cfg.grid.nonlinear = 'no';
    %     end
    cfg.mri                 = mri_realigned;
    cfg.grid.unit           = 'mm';
    subject_grid            = ft_prepare_sourcemodel(cfg);
    
    make a figure of the single subject headmodel, and grid positions
    figure; hold on;
    ft_plot_vol(headmodel, 'edgecolor', 'none', 'facecolor', 'brain','facealpha', 0.4);
    ft_plot_mesh(subject_grid.pos(subject_grid.inside,:));
    saveas(gcf,[outputdir_image filesep num2str(isubject) '_headmodelplot' ],'png');
    
    save(fname_grid, 'mri_segmented_brain', 'headmodel', 'subject_grid', 'template_grid','mri_realigned');
end


disp('Done..');
