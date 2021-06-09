
slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD
glist = [1 0 0 1 0 0 0 1 1 0 1 1 0 1 0 0 1 0 1 0 1 1]; 
% read template MRI
template_mri            = ft_read_mri('D:\fieldtrip\fieldtrip.git\trunk\external\spm8\templates\T1.nii');
template_mri.coordsys   = 'spm';

% segment template MRI for visualization (new inside)
cfg                     = [];
cfg.output              = 'brain';
cfg.template            = 'D:\fieldtrip\fieldtrip.git\trunk\external\spm8\templates\T1.nii';
template_seg            = ft_volumesegment(cfg, template_mri);

% get template grid
[mri_segmented, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(2,1,0);

for isubject = slist
    
    % load data magnetometers
    temp = load(['W:\WANDER\data\beamformer\s' num2str(isubject) '_TFR_phaselocked_V2b_source_commonfilter_MAG.mat'],'source*','MI*');
    
    % average over 10-11Hz
    FFT_MAG{isubject}                   = temp.source10_avg;
    FFT_MAG{isubject}.pow               = (temp.source10_avg.pow + temp.source11_avg.pow) ./ 2;
    FFT_MAG{isubject}.dimord            = 'pos';
    FFT_low_MAG{isubject}               = temp.source10_avg_low;
    FFT_low_MAG{isubject}.pow           = (temp.source10_avg_low.pow + temp.source11_avg_low.pow) ./ 2;
    FFT_low_MAG{isubject}.dimord        = 'pos';
    FFT_high_MAG{isubject}              = temp.source10_avg_high;
    FFT_high_MAG{isubject}.pow          = (temp.source10_avg_high.pow + temp.source11_avg_high.pow) ./ 2;
    FFT_high_MAG{isubject}.dimord       = 'pos';
    
    % replace individual pos by template
    FFT_MAG{isubject}.pos               = template_grid.pos;
    FFT_low_MAG{isubject}.pos           = template_grid.pos;
    FFT_high_MAG{isubject}.pos          = template_grid.pos;

    % difference
    FFT_diff_abs_MAG{isubject}          = FFT_MAG{isubject};
    FFT_diff_abs_MAG{isubject}.pow      = (FFT_low_MAG{isubject}.pow - FFT_high_MAG{isubject}.pow);
    FFT_diff_rel_MAG{isubject}          = FFT_MAG{isubject};
    FFT_diff_rel_MAG{isubject}.pow      = (FFT_low_MAG{isubject}.pow - FFT_high_MAG{isubject}.pow) ./ (FFT_low_MAG{isubject}.pow + FFT_high_MAG{isubject}.pow);
    
    % average over 10-11Hz
    MI_MAG{isubject}                    = temp.MI10;
    MI_MAG{isubject}.MI                 = (temp.MI10.MI + temp.MI11.MI) ./ 2;
    MI_high_MAG{isubject}               = temp.MI10_high;
    MI_high_MAG{isubject}.MI            = (temp.MI10_high.MI + temp.MI11_high.MI) ./ 2;
    MI_low_MAG{isubject}                = temp.MI10_low;
    MI_low_MAG{isubject}.MI             = (temp.MI10_low.MI + temp.MI11_low.MI) ./ 2;
    
    % replace individual pos by template
    MI_MAG{isubject}.pos                = template_grid.pos;
    MI_high_MAG{isubject}.pos           = template_grid.pos;
    MI_low_MAG{isubject}.pos            = template_grid.pos;
    
    % difference
    MI_diff_abs_MAG{isubject}           = MI_MAG{isubject};
    MI_diff_abs_MAG{isubject}.MI        = (MI_low_MAG{isubject}.MI - MI_high_MAG{isubject}.MI);
    MI_diff_rel_MAG{isubject}           = MI_MAG{isubject};
    MI_diff_rel_MAG{isubject}.MI        = (MI_low_MAG{isubject}.MI - MI_high_MAG{isubject}.MI) ./ (MI_low_MAG{isubject}.MI + MI_high_MAG{isubject}.MI);

    % load data gradiometers
    temp = load(['W:\WANDER\data\beamformer\s' num2str(isubject) '_TFR_phaselocked_V2b_source_commonfilter.mat'],'source*','MI*');
    
    % average over 10-11Hz
    FFT_GRAD{isubject}                   = temp.source10_avg;
    FFT_GRAD{isubject}.pow               = (temp.source10_avg.pow + temp.source11_avg.pow) ./ 2;
    FFT_GRAD{isubject}.dimord            = 'pos';
    FFT_low_GRAD{isubject}               = temp.source10_avg_low;
    FFT_low_GRAD{isubject}.pow           = (temp.source10_avg_low.pow + temp.source11_avg_low.pow) ./ 2;
    FFT_low_GRAD{isubject}.dimord        = 'pos';
    FFT_high_GRAD{isubject}              = temp.source10_avg_high;
    FFT_high_GRAD{isubject}.pow          = (temp.source10_avg_high.pow + temp.source11_avg_high.pow) ./ 2;
    FFT_high_GRAD{isubject}.dimord       = 'pos';
    
    % replace individual pos by template
    FFT_GRAD{isubject}.pos               = template_grid.pos;
    FFT_low_GRAD{isubject}.pos           = template_grid.pos;
    FFT_high_GRAD{isubject}.pos          = template_grid.pos;
    
    % difference
    FFT_diff_abs_GRAD{isubject}          = FFT_GRAD{isubject};
    FFT_diff_abs_GRAD{isubject}.pow      = (FFT_low_GRAD{isubject}.pow - FFT_high_GRAD{isubject}.pow);
    FFT_diff_rel_GRAD{isubject}          = FFT_GRAD{isubject};
    FFT_diff_rel_GRAD{isubject}.pow      = (FFT_low_GRAD{isubject}.pow - FFT_high_GRAD{isubject}.pow) ./ (FFT_low_GRAD{isubject}.pow + FFT_high_GRAD{isubject}.pow);
    
    % average over 10-11Hz
    MI_GRAD{isubject}                    = temp.MI10;
    MI_GRAD{isubject}.MI                 = (temp.MI10.MI + temp.MI11.MI) ./ 2;
    MI_high_GRAD{isubject}               = temp.MI10_high;
    MI_high_GRAD{isubject}.MI            = (temp.MI10_high.MI + temp.MI11_high.MI) ./ 2;
    MI_low_GRAD{isubject}                = temp.MI10_low;
    MI_low_GRAD{isubject}.MI             = (temp.MI10_low.MI + temp.MI11_low.MI) ./ 2;
    
    % replace individual pos by template
    MI_GRAD{isubject}.pos                = template_grid.pos;
    MI_high_GRAD{isubject}.pos           = template_grid.pos;
    MI_low_GRAD{isubject}.pos            = template_grid.pos;
    
    % difference
    MI_diff_abs_GRAD{isubject}               = MI_GRAD{isubject};
    MI_diff_abs_GRAD{isubject}.MI            = (MI_low_GRAD{isubject}.MI - MI_high_GRAD{isubject}.MI);
    MI_diff_rel_GRAD{isubject}               = MI_GRAD{isubject};
    MI_diff_rel_GRAD{isubject}.MI            = (MI_low_GRAD{isubject}.MI - MI_high_GRAD{isubject}.MI) ./ (MI_low_GRAD{isubject}.MI + MI_high_GRAD{isubject}.MI);
    
    % interpolate to template
    cfg = [];
    cfg.parameter                        = 'pow';
    FFT_MAG{isubject}                    = ft_sourceinterpolate(cfg, FFT_MAG{isubject},  template_mri);
    FFT_GRAD{isubject}                   = ft_sourceinterpolate(cfg, FFT_GRAD{isubject}, template_mri);
    FFT_diff_rel_MAG{isubject}           = ft_sourceinterpolate(cfg, FFT_diff_rel_MAG{isubject},  template_mri);
    FFT_diff_rel_GRAD{isubject}          = ft_sourceinterpolate(cfg, FFT_diff_rel_GRAD{isubject}, template_mri);    
    FFT_diff_abs_MAG{isubject}           = ft_sourceinterpolate(cfg, FFT_diff_abs_MAG{isubject},  template_mri);
    FFT_diff_abs_GRAD{isubject}          = ft_sourceinterpolate(cfg, FFT_diff_abs_GRAD{isubject}, template_mri);
    
    cfg.parameter                        = {'MI'};
    MI_MAG{isubject}                    = ft_sourceinterpolate(cfg, MI_MAG{isubject},  template_mri);
    MI_GRAD{isubject}                   = ft_sourceinterpolate(cfg, MI_GRAD{isubject}, template_mri);
    MI_diff_rel_MAG{isubject}           = ft_sourceinterpolate(cfg, MI_diff_rel_MAG{isubject},  template_mri);
    MI_diff_rel_GRAD{isubject}          = ft_sourceinterpolate(cfg, MI_diff_rel_GRAD{isubject}, template_mri);    
    MI_diff_abs_MAG{isubject}           = ft_sourceinterpolate(cfg, MI_diff_abs_MAG{isubject},  template_mri);
    MI_diff_abs_GRAD{isubject}          = ft_sourceinterpolate(cfg, MI_diff_abs_GRAD{isubject}, template_mri);
    
    % fix edge inside brain
    MI_MAG{isubject}.inside             = template_seg.brain;
    MI_high_MAG{isubject}.inside        = template_seg.brain;
    MI_low_MAG{isubject}.inside         = template_seg.brain;   
    MI_diff_rel_MAG{isubject}.inside    = template_seg.brain;   
    MI_diff_abs_MAG{isubject}.inside    = template_seg.brain;   
    FFT_MAG{isubject}.inside            = template_seg.brain;
    FFT_high_MAG{isubject}.inside       = template_seg.brain;
    FFT_low_MAG{isubject}.inside        = template_seg.brain;   
    FFT_diff_rel_MAG{isubject}.inside   = template_seg.brain;   
    FFT_diff_abs_MAG{isubject}.inside   = template_seg.brain;    
    MI_GRAD{isubject}.inside            = template_seg.brain;
    MI_high_GRAD{isubject}.inside       = template_seg.brain;
    MI_low_GRAD{isubject}.inside        = template_seg.brain;   
    MI_diff_rel_GRAD{isubject}.inside   = template_seg.brain;   
    MI_diff_abs_GRAD{isubject}.inside   = template_seg.brain;   
    FFT_GRAD{isubject}.inside           = template_seg.brain;
    FFT_high_GRAD{isubject}.inside      = template_seg.brain;
    FFT_low_GRAD{isubject}.inside       = template_seg.brain;   
    FFT_diff_rel_GRAD{isubject}.inside  = template_seg.brain;   
    FFT_diff_abs_GRAD{isubject}.inside  = template_seg.brain;       
end

location = [0 -70 17];
location = [-28 -29 56];
for isubject = slist

    % plot ortho
    cfg = [];
    cfg.method                          = 'ortho';
    
    cfg.funcolorlim                     = 'auto';
    cfg.funcolormap                     = hot(1000);
%     cfg.location                        = 'max';
    cfg.location                        = location;

    cfg.funparameter                    = 'pow';
   
    ft_sourceplot(cfg, FFT_MAG{isubject});
    saveas(gcf,['D:\temp\' num2str(isubject) '_MAG.jpg']);
    f1 = imread(['D:\temp\' num2str(isubject) '_MAG.jpg']);
    
    cfg.funcolorlim                     = 'maxabs';
    cfg.funcolormap                     = jet(1000);
%     cfg.location                        = 'min';
    
    cfg.funparameter                    = 'pow';
    
    ft_sourceplot(cfg, FFT_diff_abs_MAG{isubject});
    saveas(gcf,['D:\temp\' num2str(isubject) '_diff_abs_MAG.jpg']);
    f2 = imread(['D:\temp\' num2str(isubject) '_diff_abs_MAG.jpg']);
    
    ft_sourceplot(cfg, FFT_diff_rel_MAG{isubject});
    saveas(gcf,['D:\temp\' num2str(isubject) '_diff_rel_MAG.jpg']);
    f3 = imread(['D:\temp\' num2str(isubject) '_diff_rel_MAG.jpg']);

    cfg.funcolorlim                     = 'auto';    
    cfg.funcolormap                     = hot(1000);
%     cfg.location                        = 'max';
    
    cfg.funparameter                    = 'MI';
    
    ft_sourceplot(cfg, MI_MAG{isubject});
    saveas(gcf,['D:\temp\MI_' num2str(isubject) '_MAG.jpg']);
    f4 = imread(['D:\temp\MI_' num2str(isubject) '_MAG.jpg']);
    
    cfg.funcolorlim                     = 'maxabs';
    cfg.funcolormap                     = jet(1000);
%     cfg.location                        = 'min';
    
    cfg.funparameter                    = 'MI';
    
    ft_sourceplot(cfg, MI_diff_abs_MAG{isubject});
    saveas(gcf,['D:\temp\MI_' num2str(isubject) '_diff_abs_MAG.jpg']);
    f5 = imread(['D:\temp\MI_' num2str(isubject) '_diff_abs_MAG.jpg']);
    
    ft_sourceplot(cfg, MI_diff_rel_MAG{isubject});
    saveas(gcf,['D:\temp\MI_' num2str(isubject) '_diff_rel_MAG.jpg']);
    f6 = imread(['D:\temp\MI_' num2str(isubject) '_diff_rel_MAG.jpg']);
    
    cfg.funcolorlim                     = 'auto';    
    cfg.funcolormap                     = hot(1000);
%     cfg.location                        = 'max';
    
    cfg.funparameter                    = 'pow';
    
    ft_sourceplot(cfg, FFT_GRAD{isubject});
    saveas(gcf,['D:\temp\' num2str(isubject) '_GRAD.jpg']);
    f7 = imread(['D:\temp\' num2str(isubject) '_GRAD.jpg']);

    cfg.funcolorlim                     = 'maxabs';
    cfg.funcolormap                     = jet(1000);    
%     cfg.location                        = 'min';
    
    cfg.funparameter                    = 'pow';
        
    ft_sourceplot(cfg, FFT_diff_abs_GRAD{isubject});
    saveas(gcf,['D:\temp\' num2str(isubject) '_diff_abs_GRAD.jpg']);
    f8 = imread(['D:\temp\' num2str(isubject) '_diff_abs_GRAD.jpg']);
    
    ft_sourceplot(cfg, FFT_diff_rel_GRAD{isubject});
    saveas(gcf,['D:\temp\' num2str(isubject) '_diff_rel_GRAD.jpg']);
    f9 = imread(['D:\temp\' num2str(isubject) '_diff_rel_GRAD.jpg']);
    
    cfg.funcolorlim                     = 'auto';
    cfg.funcolormap                     = hot(1000);
%     cfg.location                        = 'max';
    
    cfg.funparameter                    = 'MI';
    
    ft_sourceplot(cfg, MI_GRAD{isubject});
    saveas(gcf,['D:\temp\MI_' num2str(isubject) '_GRAD.jpg']);
    f10 = imread(['D:\temp\MI_' num2str(isubject) '_GRAD.jpg']);
    
    cfg.funcolorlim                     = 'maxabs';
    cfg.funcolormap                     = jet(1000);
%     cfg.location                        = 'min';
    
    cfg.funparameter                    = 'MI';

    ft_sourceplot(cfg, MI_diff_abs_GRAD{isubject});
    saveas(gcf,['D:\temp\MI_' num2str(isubject) '_diff_abs_GRAD.jpg']);
    f11 = imread(['D:\temp\MI_' num2str(isubject) '_diff_abs_GRAD.jpg']);
    
    ft_sourceplot(cfg, MI_diff_rel_GRAD{isubject});
    saveas(gcf,['D:\temp\MI_' num2str(isubject) '_diff_rel_GRAD.jpg']);
    f12 = imread(['D:\temp\MI_' num2str(isubject) '_diff_rel_GRAD.jpg']);

    figure;
    subplot(2,3,1); image(f1); title('A mag');          axis square tight off
    subplot(2,3,2); image(f2); title('A abs. mag');     axis square tight off
    subplot(2,3,3); image(f3); title('A rel. mag');     axis square tight off
    subplot(2,3,4); image(f4); title('MI mag');         axis square tight off
    subplot(2,3,5); image(f5); title('MI abs. mag');    axis square tight off
    subplot(2,3,6); image(f6); title('MI rel. mag');    axis square tight off    
    print(['D:\temp\summary_' num2str(isubject) '_MAG'],'-dpng','-r1200')

    figure;
    subplot(2,3,1); image(f7);  title('A grad');        axis square tight off
    subplot(2,3,2); image(f8);  title('A abs. grad');   axis square tight off
    subplot(2,3,3); image(f9);  title('A rel. grad');   axis square tight off
    subplot(2,3,4); image(f10); title('MI grad');       axis square tight off
    subplot(2,3,5); image(f11); title('MI abs. grad');  axis square tight off
    subplot(2,3,6); image(f12); title('MI rel. grad');  axis square tight off
    print(['D:\temp\summary_' num2str(isubject) '_GRAD'],'-dpng','-r1200')
  
    close all
    
end

% GA
mm = ft_read_mri('D:\analysis\WANDER\scripts\MNI_template\mni_icbm152_t1_tal_nlin_asym_09c.nii');



cfg = [];
cfg.parameter               = 'MI';
MI_MAG_GA                   = ft_sourcegrandaverage(cfg,MI_MAG{slist});
MI_GRAD_GA                  = ft_sourcegrandaverage(cfg,MI_GRAD{slist});
MI_MAG_GA.dimord            = 'pos';
MI_GRAD_GA.dimord           = 'pos';

cfg = [];
cfg.parameter               = 'pow';
FFT_MAG_GA                  = ft_sourcegrandaverage(cfg,FFT_MAG{slist});
FFT_GRAD_GA                 = ft_sourcegrandaverage(cfg,FFT_GRAD{slist});
FFT_MAG_GA.dimord           = 'pos';
FFT_GRAD_GA.dimord          = 'pos';

cfg = [];
cfg.parameter                = 'MI';
MI_MAG_GA                    = ft_sourceinterpolate(cfg, MI_MAG_GA,  mm);
MI_GRAD_GA                   = ft_sourceinterpolate(cfg, MI_GRAD_GA,  mm);

cfg = [];
cfg.parameter               = 'pow';
FFT_MAG_GA                  = ft_sourceinterpolate(cfg, FFT_MAG_GA,  mm);
FFT_GRAD_GA                 = ft_sourceinterpolate(cfg, FFT_GRAD_GA,  mm);


cfg = [];
cfg.parameter               = 'pow';
FFT_MAG_GA                  = ft_sourceinterpolate(cfg, FFT_MAG_GA,  template_mri);
FFT_GRAD_GA                 = ft_sourceinterpolate(cfg, FFT_GRAD_GA,  template_mri);

% plot ortho
cfg = [];
cfg.method                          = 'ortho';
cfg.funcolorlim                     = 'auto';
cfg.funcolormap                     = jet(1000);
cfg.location                        = 'max';
cfg.funparameter                    = 'MI';
ft_sourceplot(cfg, MI_MAG_GA);
ft_sourceplot(cfg, MI_GRAD_GA);

cfg.funparameter                    = 'pow';
ft_sourceplot(cfg, FFT_MAG_GA);
ft_sourceplot(cfg, FFT_GRAD_GA);


% calculate max / ALSO ON SINGLE SUBJECT MRIs

atlas       = ft_read_atlas('d:/fieldtrip-master_21092017/fieldtrip-master/template/atlas/aal/ROI_MNI_V4.nii');
% 
% % interpolate to template
% cfg             = [];
% cfg.parameter   = 'tissue';
% cfg.interpmethod = 'nearest';
% atlas_int       = ft_sourceinterpolate(cfg, atlas,  mm);

for isubject = slist
%     FFT_MAG{isubject}.pow(template_seg.brain) = nan;
    [~,pos(isubject)] = max(FFT_GRAD{isubject}.pow);
    
    cfg = [];
    cfg.method                          = 'ortho';
    cfg.funcolorlim                     = 'auto';
    cfg.funcolormap                     = jet(1000);
    cfg.location                        = FFT_GRAD{isubject}.pos(pos(isubject),:);
%     cfg.atlas                           = atlas;
    ft_sourceplot(cfg, mm);
    
    saveas(gcf,['D:\temp\temp' num2str(isubject) '.jpg']);
    close all
end

% get label
for isubject = slist
    tissue{isubject} = atlas.tissue(pos(isubject));
    if tissue{isubject} ~= 0
        label{isubject} = atlas.tissuelabel(tissue{isubject});
    else
        label{isubject} = 'EMPTY';
    end
end

% find calcarine or lingual
for isubject = slist
    label{isubject}
    
    q = strfind(label{isubject},'Calcarine','ForceCellOutput',true);
    if q{1}
        calcerine(isubject) = true;
    else
        calcerine(isubject) = false;
    end
    
    q = strfind(label{isubject},'Lingual','ForceCellOutput',true);
    if q{1}
        lingual(isubject) = true;
    else
        lingual(isubject) = false;
    end
    
    q = strfind(label{isubject},'Cuneus','ForceCellOutput',true);
    if q{1}
        cuneus(isubject) = true;
    else
        cuneus(isubject) = false;
    end  
    
    q = strfind(label{isubject},'Precuneus','ForceCellOutput',true);
    if q{1}
        precuneus(isubject) = true;
    else
        precuneus(isubject) = false;
    end 
end

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD
slistbool = false(1,26);
slistbool(slist) = true;

sum(slistbool)
sum(calcerine)
sum(lingual)
sum(cuneus)
sum(precuneus)
slist = find(calcerine | lingual | cuneus | precuneus);

figure;
i = 1;
for isubject = slist
    f = imread(['D:\temp\temp' num2str(isubject) '.jpg']);
    subplot(5,6,i);
    image(f);
    i = i + 1;
    title(num2str(isubject));
end
    print(['D:\temp\sources_overview.jpg'],'-dpng','-r1200')
    
figure;
ft_plot_ortho(mm.anatomy,'location',[FFT_GRAD{isubject}.pos(pos(isubject),:)])

figure;
ft_plot_ortho(mm.anatomy,'transform',eye(4),'location',[FFT_GRAD{isubject}.pos(pos(isubject),:)])

atlas       = ft_read_atlas('ROI_MNI_V4.nii');
atlas       = ft_read_atlas('d:/fieldtrip-master_21092017/fieldtrip-master/template/atlas/afni/TTatlas+tlrc.BRIK');
atlas       = ft_read_atlas('d:/fieldtrip-master_21092017/fieldtrip-master/template/atlas/aal/ROI_MNI_V4.nii');

    atlas = ft_convert_coordsys(atlas, 'acpc')

    

mm.coordsys = 'mni';

cfg = [];
cfg.funparameter                       = 'pow';
cfg.method                          = 'ortho';
cfg.funcolorlim                     = 'auto';
cfg.funcolormap                     = jet(1000);
cfg.location                        = FFT_GRAD{isubject}.pos(pos(isubject),:);
cfg.atlas                           = atlas;
ft_sourceplot(cfg, FFT_GRAD{isubject});
ft_sourceplot(cfg, mm);

cfg         = [];
cfg.roi     = [FFT_GRAD{isubject}.pos(pos(isubject),:)];
cfg.sphere  = 1;
labels      = ft_volumelookup(cfg,FFT_GRAD{isubject});

