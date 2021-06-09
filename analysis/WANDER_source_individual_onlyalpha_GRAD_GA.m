function WANDER_source_individual_onlyalpha_GRAD_GA(rootpath)

if rootpath == 1
    datapath = 'z:';
else
    datapath = '/shared/projects/project_wander';
end

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

% Do this once with standard template for stats, and with 1mm for image

% t = 2.086;

load([datapath,'/WANDER/scripts/grad.mat'],'grad');
neigh_cmb = load('neuromag306cmb_neighb');
neigh_mag = load('neuromag306mag_neighb_last');

% read template MRI
template_mri            = ft_read_mri([datapath,'/WANDER/fieldtrip/external/spm8/templates/T1.nii']);
template_mri.coordsys   = 'spm';

% segment template MRI for visualization
cfg                     = [];
cfg.output              = 'brain';
cfg.template            = [datapath,'/WANDER/fieldtrip/external/spm8/templates/T1.nii'];
template_seg            = ft_volumesegment(cfg, template_mri);

% get template grid
[mri_segmented, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(1,rootpath,0);  
    
% define inside according to segmentation
template_mri.inside     = template_seg.brain;


% load single subject data
for isubject = slist

    [~, ~, ~, source{isubject}, source_low{isubject}, source_high{isubject}] = WANDER_source_individual_onlyalpha_GRAD(isubject,0,rootpath);
       
    source_high{isubject}.pos   = template_grid.pos;
    source_low{isubject}.pos    = template_grid.pos;
    source{isubject}.pos        = template_grid.pos;
   
    source_high{isubject}.freq  = 10;
    source_low{isubject}.freq   = 10;
    source{isubject}.freq       = 10;
    
    source_high{isubject}.pow   = source_high{isubject}.avg.pow;
    source_high{isubject}       = rmfield(source_high{isubject},'avg');
    source_low{isubject}.pow    = source_low{isubject}.avg.pow;
    source_low{isubject}        = rmfield(source_low{isubject},'avg');

    source_diff{isubject}       = source_low{isubject};
    source_diff{isubject}.pow   = (source_low{isubject}.pow - source_high{isubject}.pow) ./ (source_low{isubject}.pow + source_high{isubject}.pow);
    source_diff{isubject}.pos   = template_grid.pos; 

end


% mri_1mm = ft_read_mri('d:/fieldtrip-master_21092017/fieldtrip-master/template/anatomy/single_subj_T1_1mm.nii');
% mri_1mm.coordsys = 'spm';
% 
% cfg                     = [];
% cfg.output              = 'brain';
% mri_1mm_seg             = ft_volumesegment(cfg, mri_1mm);
% 

for isubject = slist
    
%     try 
%         source_diff{isubject}.pow   = source_diff{isubject}.avg.pow;
%         source_diff{isubject}       = rmfield(source_diff{isubject},'avg');
%     
%         source_high{isubject}       = rmfield(source_high{isubject},'method');
%         source_low{isubject}        = rmfield(source_low{isubject},'method');
%         source_diff{isubject}       = rmfield(source_diff{isubject},'method');
%     catch
%     end
    cfg = [];
    cfg.parameter                       = {'pow'};
%     source_high_int{isubject}           = ft_sourceinterpolate(cfg, source_high{isubject}, template_mri);
%     source_low_int{isubject}            = ft_sourceinterpolate(cfg, source_low{isubject}, template_mri);
    source_diff_int{isubject}            = ft_sourceinterpolate(cfg, source_diff{isubject}, template_mri);
    
    source_high_int{isubject}.inside    = template_seg.brain;
    source_low_int{isubject}.inside     = template_seg.brain;
    source_diff_int{isubject}.inside     = template_seg.brain;

end

cfg = [];
source_diff_int_GA          = ft_sourcegrandaverage(cfg,source_diff_int{slist});

% plot interpolated diff, no mask, no stats
cfg = [];
cfg.method                  = 'ortho';
cfg.funparameter            = 'pow';
ft_sourceplot(cfg, source_diff_int_GA);

%% statistics
cfg = [];
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 4000;
cfg.alpha               = 0.025; % threshold second-level
cfg.clusteralpha        = 0.001; % voxel level
cfg.tail                = 0;
cfg.design(1,:)         = [1:length(slist) 1:length(slist)];
cfg.design(2,:)         = [ones(size(slist)) ones(size(slist))*2];
cfg.uvar                = 1; 
cfg.ivar                = 2;
cfg.parameter           = 'pow';

% stats on non-interpolated data
stat_FFT                = ft_sourcestatistics(cfg, source_low{slist},source_high{slist});
stat_FFT.mask1          = stat_FFT.negclusterslabelmat == 1;
stat_FFT.mask2          = stat_FFT.negclusterslabelmat == 2;
stat_FFT.mask3          = stat_FFT.negclusterslabelmat == 3;

% save results
save(fullfile(datapath,'/WANDER/data/stat/stat_FFT'),'stat_FFT');
 
% stats on interpolated data
int_stat_FFT            = ft_sourcestatistics(cfg, source_low_int{slist},source_high_int{slist});
int_stat_FFT.anatomy    = template_mri.anatomy;
int_stat_FFT.diff       = source_diff_int_GA.pow;
int_stat_FFT.mask1      = int_stat_FFT.negclusterslabelmat == 1;
int_stat_FFT.mask2      = int_stat_FFT.negclusterslabelmat == 2;
int_stat_FFT.mask3      = int_stat_FFT.negclusterslabelmat == 3;
% int_stat_FFT.stat(~template_mri.inside) = nan;

% save results
save(fullfile(datapath,'/WANDER/data/stat/int_stat_FFT'),'int_stat_FFT');

% interpolate stats from non-interpolated data
% cfg = [];
% cfg.parameter           = {'stat','mask','mask1','mask2','mask3'};
% cfg.interpmethod        = 'nearest'; 
% stat_FFT_int            = ft_sourceinterpolate(cfg, stat_FFT, template_mri);
% stat_FFT_int.inside     = template_mri.inside;

% interpolate interpolated stats on higher resoluation
mri_1mm = ft_read_mri([datapath,'/WANDER/fieldtrip/template/anatomy/single_subj_T1_1mm.nii']);

cfg = [];
cfg.parameter           = {'stat','mask','mask1','mask2','mask3','diff'};
% cfg.interpmethod        = 'nearest'; % for making table later
int_stat_FFT_1mm        = ft_sourceinterpolate(cfg, int_stat_FFT, mri_1mm);
% int_stat_FFT_1mm.inside = mri_1mm.anatomy ~= 0;

% 
% cfg                     = [];
% cfg.output              = 'brain';
% mri_1mm_seg             = ft_volumesegment(cfg, mri_1mm);

% % interpolate to template MRI
% cfg = [];
% cfg.parameter           = {'stat','mask'};
% stat_FFT_1mm            = ft_sourceinterpolate(cfg, stat_FFT, mri_1mm);
% stat_FFT_1mm.stat       = stat_FFT_1mm.stat(mri_1mm_seg.brain == 1);


cfg = [];
cfg.method                  = 'ortho';
cfg.funcolormap             = parula(1000);
cfg.funparameter            = 'stat';
% cfg.maskparameter           = 'mask1';
cfg.location = 'min';
cfg.crosshair = 'no';
cfg.interactive = 'yes';
% cfg.colormap = parula(1000);
% cfg.opacitylim = [0 1];
% 
% temp = stat_FFT_int;
% temp.stat(~temp.mask1) = nan;
% ft_sourceplot(cfg, temp);

% temp = stat_FFT_int;
% temp.stat(~temp.mask1) = nan;
% ft_sourceplot(cfg, temp);

% no mask: high res for image
% cfg.funparameter            = 'diff';
% ft_sourceplot(cfg, int_stat_FFT);

% mask: high res for image
temp = int_stat_FFT_1mm;
cfg.maskparameter           = 'mask1';
% 
temp.stat(~temp.mask1) = nan;
% temp.diff(~temp.mask) = nan;
% cfg.funparameter       = 'diff';
ft_sourceplot(cfg, temp);
colormap parula(1000)

saveas(gcf,fullfile(datapath,'source_alpha_ortho'),'fig');
print(gcf,fullfile(datapath,'source_alpha_ortho'),'-dpdf');

h = gcf;
h.Renderer = 'Painters';
print(gcf,fullfile(datapath,'source_alpha_ortho_painters'),'-dpdf');


% stat_FFT_1mm.mask = stat_FFT_1mm.mask * 0.5;
% ft_sourceplot(cfg, stat_FFT_1mm);


% 
% cfg = [];
% cfg.method          = 'surface';
% cfg.surfinflated    = 'surface_inflated_both.mat'; % Cortical sheet from canonical MNI brain
% % cfg.surffile        = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
% % cfg.surffile        = 'surface_pial_right.mat';
% % cfg.opacitymap    = 'stat';
% % cfg.opacitymap    = 'rampup';
% cfg.maskparameter   = 'mask';
% cfg.funparameter    = 'stat';
% % cfg.surfdownsample  = 10;
% ft_sourceplot(cfg, stat_FFT_int);
% ft_sourceplot(cfg, int_stat_FFT);
% view(60,30);
% material dull
% saveas(gcf,'C:\Users\Admin\Dropbox\WANDER_article\source_alpha.pdf');




% view(0,90)
% 
% 
% material dull
% material shiny

%%


%%%%
aal = ft_read_atlas(fullfile(datapath,'/WANDER/aal/ROI_MNI_V5.nii'));

% template_mri = ft_read_mri('D:\fieldtrip\fieldtrip.git\trunk\external\spm8\templates\T1.nii');

% make table
clear target tabsource tabtemp
GA_int = int_stat_FFT;
stat_FFT = int_stat_FFT;
% stat_FFT = stat_FFT_int;



GA_int.anatomy = template_mri.anatomy;
GA_int.negclusterslabelmat = reshape(GA_int.negclusterslabelmat,size(GA_int.anatomy));
GA_int.stat = reshape(GA_int.stat,size(GA_int.anatomy));
GA_int.mask = reshape(GA_int.mask,size(GA_int.anatomy));
GA_int.prob = reshape(GA_int.prob,size(GA_int.anatomy));

tabsource = [];
for isource = 1 : size(unique(stat_FFT.negclusterslabelmat(~isnan(stat_FFT.negclusterslabelmat))),1)
    
    disp(['working on source ',num2str(isource),' of ',num2str(size(unique(stat_FFT.negclusterslabelmat(~isnan(stat_FFT.negclusterslabelmat))),1))]);
    itissue = 0; 
    
    for i = 1 : 120
        
        overlap = false([91,109,91]);
        overlap(GA_int.negclusterslabelmat == isource & aal.tissue == i) = true;
        
        if ~isempty(find(overlap == 1, 1))
            
            itissue = itissue + 1;
            target{isource}.tissuenr(itissue)                       = i;
            target{isource}.sourcecount(itissue)                    = size(find(overlap == 1),1);
            target{isource}.stat_mean(itissue)                      = nanmean(GA_int.stat(overlap));
            target{isource}.stat_min(itissue)                       = min(GA_int.stat(overlap));
            target{isource}.stat_max(itissue)                       = max(GA_int.stat(overlap));
%             target{isource}.stat_firstlevel_mean(itissue)           = nanmean(GA_int.stat_firstlevel(overlap));
%             target{isource}.stat_firstlevel_min(itissue)            = min(GA_int.stat_firstlevel(overlap));
%             target{isource}.stat_firstlevel_max(itissue)            = max(GA_int.stat_firstlevel(overlap));
            
%             target{isource}.stat_mean_correct(itissue)              = nanmean(correct_int.prob(overlap));
%             target{isource}.stat_mean_RT(itissue)                   = nanmean(RT_int.prob(overlap));
%             target{isource}.stat_mean_trialduration(itissue)        = nanmean(trialduration_int.prob(overlap));

%             target{isource}.stat_cluster_correct(itissue)           = nanmean(correct_int.prob(GA_int.negclusterslabelmat == isource));
%             target{isource}.stat_cluster_RT(itissue)                = nanmean(RT_int.prob(GA_int.negclusterslabelmat == isource));
%             target{isource}.stat_cluster_trialduration(itissue)     = nanmean(trialduration_int.prob(GA_int.negclusterslabelmat == isource));
            
            temp = GA_int.stat;
            temp(overlap == false) = NaN;
            [Y,I] = nanmin(temp(:));
            [I1,I2,I3] = ind2sub([91,109,91],I);
            
            target{isource}.stat_loc_min(itissue)                   = I;
            XYZ = ft_warp_apply(template_mri.transform, [I1,I2,I3]) ;

            target{isource}.stat_loc_min1(itissue)                  = XYZ(1);
            target{isource}.stat_loc_min2(itissue)                  = XYZ(2);
            target{isource}.stat_loc_min3(itissue)                  = XYZ(3);
 
            target{isource}.clusterprob(itissue)                    = nanmean(GA_int.prob(GA_int.negclusterslabelmat == isource & aal.tissue == i));
            target{isource}.regioncount(itissue)                    = size(find(aal.tissue == i),1);
            
        end
    end
end

% remove empty cells
temp = target;
clear target
ii = 0;
for i = 1 : size(temp,2)
    if ~isempty(temp{i})
        ii = ii + 1;
        target{ii} = temp{i};
    end
end
if rootpath == 1
    datapath = 'z:';
else
    datapath = '/shared/projects/project_wander';
end

tabsource = [];
for isource = 1 : size(target,2)
    clear clusterprob tissuenr loc loc1 loc2 loc3 stat_min stat_max stat_firstlevel_mean sourcecount sourceperc tissuelabel cluster stat_mean prob_mean regionperc clusterperc clustercount stat_mean_correct stat_mean_RT stat_mean_trialduration stat_cluster_correct stat_cluster_RT stat_cluster_trialduration

    for i = 1 : size(target{isource}.sourcecount,2)
        
        clustercount(i)         = target{isource}.sourcecount(i);
        clusterperc(i)          = target{isource}.sourcecount(i) / sum(target{isource}.sourcecount) * 100;
        regionperc(i)           = target{isource}.sourcecount(i) / target{isource}.regioncount(i) * 100;    
        tissuelabel{i}          = aal.tissuelabel{target{isource}.tissuenr(i)};
        tissuenr(i)             = target{isource}.tissuenr(i);
        cluster(i)              = isource;
        stat_mean(i)            = target{isource}.stat_mean(i);
        stat_min(i)             = target{isource}.stat_min(i);
        stat_max(i)             = target{isource}.stat_max(i);
        loc1(i)                 = target{isource}.stat_loc_min1(i);
        loc2(i)                 = target{isource}.stat_loc_min2(i);
        loc3(i)                 = target{isource}.stat_loc_min3(i);
        loc(i)                  = target{isource}.stat_loc_min(i); 
        clusterprob(i)          = target{isource}.clusterprob(i);
%         stat_firstlevel_mean(i) = target{isource}.stat_firstlevel_mean(i);       
%         
%         stat_mean_correct(i)       = target{isource}.stat_mean_correct(i);
%         stat_mean_RT(i)            = target{isource}.stat_mean_RT(i);
%         stat_mean_trialduration(i) = target{isource}.stat_mean_trialduration(i);
%         
%         stat_cluster_correct(i)       = target{isource}.stat_cluster_correct(i);
%         stat_cluster_RT(i)            = target{isource}.stat_cluster_RT(i);
%         stat_cluster_trialduration(i) = target{isource}.stat_cluster_trialduration(i);
    end
%     
%     cluster = cluster';
%     tissuelabel = tissuelabel';
%     clusterperc = clusterperc';
%     clustercount = clustercount';
%     stat_mean = stat_mean';
%     stat_min = stat_min';
%     stat_max = stat_max';
%     prob_mean = prob_mean';
%     regionperc = regionperc';
%     loc1 = loc1';
%     loc2 = loc2';
%     loc3 = loc3';
%     loc = loc';
%     stat_firstlevel_mean = stat_firstlevel_mean';

    tabtemp = table(cluster',clusterprob',tissuelabel',stat_min',stat_mean',clustercount',regionperc', loc1', loc2', loc3','VariableNames',{'Cluster','Clusterprob','AAL_region','Peak_t','Mean_t','Size','Percent_activation','X','Y','Z'});
%     tabtemp = table(cluster,tissuelabel,regionperc,clusterperc,clustercount,stat_mean,stat_min,stat_max,prob_mean,stat_firstlevel_mean, loc, loc1, loc2, loc3);

    
    tabsource = [tabsource; tabtemp];
        clear tabtemp

end

% prune table/.
tabsource(tabsource.Cluster > 4,:) = [];
tabsource(tabsource.Percent_activation < 1,:) = [];
% tabsource = sortrows(tabsource,{'Cluster','Size'},{'ascend','descend'})
tabsource = sortrows(tabsource,{'Cluster','Mean_t'},{'ascend','ascend'})
% tabsource = sortrows(tabsource,{'regionperc','prob_mean'},{'descend','ascend'})
% tabsource = sortrows(tabsource,{'stat_mean'},{'descend'})

tabsource(:,1:10)

save(fullfile(datapath,'WANDER/data/stat/int_stat_FFT_table'),'tabsource');

load(fullfile(datapath,'WANDER/data/stat/int_stat_FFT_table'),'tabsource');




