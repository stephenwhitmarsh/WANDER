slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

force = 0;
rootpath = 1;
restingstate = 0;

for isubject = slist  
    disp(['subject ', num2str(isubject)]);    
    fname_SSEF_phase = ['w:\WANDER\data\SSEF\SSEF_phase_' num2str(isubject) '.mat'];

    % load directly from file
    temp = load(fname_SSEF_phase,'data_MEG_stim_avg');
    data_MEG_stim_avg{isubject} = temp.data_MEG_stim_avg;

    % load headmodel
    [mri_segmented{isubject}, headmodel{isubject}, subject_grid{isubject}, template_grid, mri_realigned{isubject}] = WANDER_grid(isubject,rootpath,0);
    
    % load sensor information
    [dataset_task, ~] = WANDER_subjectinfo;
    hdr{isubject} = ft_read_header(dataset_task{isubject,1});
%     
%     cfg = [];
%     cfg.latency = [0.04 0.06];
%     cfg.avgovertime = 'yes';
%     data_MEG_stim_avg{isubject} = ft_selectdata(cfg,data_MEG_stim_avg{isubject});
    
    % dipole fit
    cfg = [];
    cfg.numdipoles = 1;
    cfg.grid        = subject_grid{isubject};
%     cfg.grid.pos    = template_grid.pos;   
    cfg.headmodel   = headmodel{isubject};
    cfg.grad        = hdr{isubject}.grad;
    cfg.channel     = {'MEG*2','MEG*3'};
    cfg.model = 'regional';
    cfg.nonlinear   = 'no';
    source_mag{isubject} = ft_dipolefitting(cfg, data_MEG_stim_avg{isubject});
end

% read template MRI
% template_mri            = ft_read_mri('D:\fieldtrip\fieldtrip.git\trunk\external\spm8\templates\T1.nii');
% template_mri.coordsys   = 'spm';

template_mri            = ft_read_mri('D:\analysis\WANDER\scripts\colin27_t1_tal_lin.nii');
% template_mri.coordsys   = 'spm';
% 
% cfg = [];
% cfg.resolution = 1;
% cfg.xrange = [-100 100];
% cfg.yrange = [-110 140];
% cfg.zrange = [-80 120];
% template_resliced = ft_volumereslice(cfg, template_mri);
% template_mri = ft_convert_units(template_mri, 'mm');

for isubject = slist
    
%     % read template MRI
%     template_mri            = ft_read_mri('D:\fieldtrip\fieldtrip.git\trunk\external\spm8\templates\T1.nii');
%     template_mri.coordsys   = 'spm';
%     template_mri.transform  = inv(mri_realigned{isubject}.transform);
%     
    pos_subject{isubject} = source_mag{isubject}.dip.pos;
    ori_subject{isubject} = source_mag{isubject}.dip.mom;
%     cfg = [];
%     cfg.location = pos_subject{isubject};
%     ft_sourceplot(cfg, mri_realigned{isubject});
%     
    index = find(subject_grid{isubject}.pos(:,1) == source_mag{isubject}.dip.pos(:,1) & subject_grid{isubject}.pos(:,2) == source_mag{isubject}.dip.pos(:,2) & subject_grid{isubject}.pos(:,3) == source_mag{isubject}.dip.pos(:,3));
    pos_template{isubject} = template_grid.pos(index,:);
    ori_template{isubject} = source_mag{isubject}.dip.mom;
    
%     cfg = [];
%     cfg.location = pos_template{isubject}; 
%     ft_sourceplot(cfg, template_mri);
%     
 
end


i = 1;
fig = figure;
for isubject = slist
    subplot(5,5,i); hold;
    ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos_template{isubject}, 'plotmarker' ,pos_template{isubject},'markersize',30,'markercolor',[1 1 1])

    view(0,90);
    i = i + 1;
%     title(num2str(isubject));
end
print(fig,'-dpng','-r300',['d:\analysis\WANDER\images\SSSEF_source_template.png']);

i = 1;
figure;
for isubject = slist
    subplot(5,5,i);
    ft_plot_slice(mri_realigned{isubject}.anatomy, 'transform', mri_realigned{isubject}.transform, 'location', pos_subject{isubject}, 'plotmarker' ,pos_subject{isubject},'markersize',30)
    view(0,90);
    i = i + 1;
end
print(gcf,'-dpng','-r300',['d:\analysis\WANDER\images\SSSEF_source_subject.png']);



ft_plot_dipole(pos_template, ori_template, 'color', 'r','diameter',1)

ft_plot_slice(mri_realigned{isubject}.anatomy, 'transform', mri_realigned{isubject}.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_realigned{isubject}.anatomy, 'transform', mri_realigned{isubject}.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)


ft_plot_dipole(source_mag.dip.pos(2,:), mean(source_mag.dip.mom(4:6,:),2), 'color', 'r')
ft_plot_dipole(source_planar.dip.pos(1,:), mean(source_planar.dip.mom(1:3,:),2), 'color', 'g')
ft_plot_dipole(source_planar.dip.pos(2,:), mean(source_planar.dip.mom(4:6,:),2), 'color', 'g')
pos = mean(source_mag.dip.pos,1);
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [1 0 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 1 0], 'resolution', 0.1)
ft_plot_slice(mri_resliced_cm.anatomy, 'transform', mri_resliced_cm.transform, 'location', pos, 'orientation', [0 0 1], 'resolution', 0.1)
ft_plot_crosshair(pos, 'color', [1 1 1]/2);
axis tight
axis off
view(12, -10)