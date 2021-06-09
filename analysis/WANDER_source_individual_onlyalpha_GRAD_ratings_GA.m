slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD
% slist = [2     4     5     8     9    11    12    13    16    17    18    19    20    23    24    25    26]; % without those with strange alpha localization

% t = 2.086;

load('D:\analysis\WANDER\scripts\grad','grad');
neigh_cmb = load('neuromag306cmb_neighb');
neigh_mag = load('neuromag306mag_neighb_last');

% addpath('D:\fieldtrip\fieldtrip.git\trunk');
addpath('D:/analysis/WANDER/scripts/');
ft_defaults

% read template MRI
template_mri            = ft_read_mri('D:\fieldtrip\fieldtrip.git\trunk\external\spm8\templates\T1.nii');
template_mri.coordsys   = 'spm';

% segment template MRI for visualization
cfg                     = [];
cfg.output              = 'brain';
cfg.template            = 'D:\fieldtrip\fieldtrip.git\trunk\external\spm8\templates\T1.nii';
template_seg            = ft_volumesegment(cfg, template_mri);


% define inside according to segmentation
template_mri.inside     = template_seg.brain;


mri_1mm = ft_read_mri('d:/fieldtrip-master_21092017/fieldtrip-master/template/anatomy/single_subj_T1_1mm.nii');
mri_1mm.coordsys = 'spm';

cfg                     = [];
cfg.output              = 'brain';
mri_1mm_seg             = ft_volumesegment(cfg, mri_1mm);

% load single subject data
for isubject = slist  
    [~, source{isubject}] = WANDER_source_individual_onlyalpha_GRAD_ratings(isubject,0,1);
    [mri_segmented, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(isubject,1,0);
    
    for irating = 1 : 4     
        cfg = [];
        cfg.parameter                           = {'pow'};
        source{isubject}{irating}.pos           = template_grid.pos;
        source_int{isubject}{irating}           = ft_sourceinterpolate(cfg, source{isubject}{irating}, template_mri);
        source_int{isubject}{irating}.freq      = 10;
    end
%     
% % get template grid
% [mri_segmented, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(isubject,1,0);
% source{isubject}{1}.anatomy = 
% 
% cfg = [];
% cfg.method                  = 'ortho';
% % cfg.funcolormap             = hot(1000);
% cfg.funparameter            = 'pow';
% % cfg.funcolorlim = [-60 -49];
% ft_sourceplot(cfg, source_int{isubject}{1});

end




load('w:/WANDER/data/stat/int_stat_FFT');

roi_ss = find(int_stat_FFT.negclusterslabelmat == 1);
roi_v  = find(int_stat_FFT.negclusterslabelmat == 2);

for irating = 1 : 4
    i = 1;
    for isubject = slist
        pow_ss(irating,i) = nanmean(source_int{isubject}{irating}.pow(roi_ss));
        pow_v(irating,i) = nanmean(source_int{isubject}{irating}.pow(roi_v));
        i = i + 1;
    end
end


pow_ss_avg = mean(zscore(pow_ss)');
pow_ss_std = std(zscore(pow_ss)') / sqrt(22);

pow_v_avg = mean(zscore(pow_v)');
pow_v_std = std(zscore(pow_v)') / sqrt(22);


figure; hold;
plot(pow_ss_avg + pow_ss_std,'k');
plot(pow_ss_avg,'k');
plot(pow_ss_avg - pow_ss_std,'k');


fig = figure; hold;
errorbar(pow_ss_avg,pow_ss_std);
errorbar(pow_v_avg,pow_v_std);
axis([0,5,-1,1]);
xticks([0:5])
xticklabels({'','High','','','Low'});
xlabel('Attention');
ylabel('Normalized alpha (SEM)');
legend({'Somatosensory','Visual'},'location','East');
saveas(fig,'W:/WANDER/images/article/source_alpha_ratings','png');
saveas(fig,'W:/WANDER/images/article/source_alpha_ratings','svg');

[p,t,stats] = anova1(zscore(pow_ss)') ;
[c,m,h,nms] = multcompare(stats);

[p,t,stats] = anova1(zscore(pow_v)') ;
[c,m,h,nms] = multcompare(stats);

figure; hold;
plot(roi_ss_val)
plot(roi_v_val)



for isubject = slist
    %     source_diff{isubject}.pow   = source_diff{isubject}.avg.pow;
    %     source_diff{isubject}       = rmfield(source_diff{isubject},'avg');
    %     source_high{isubject}.pow   = source_high{isubject}.avg.pow;
    %     source_high{isubject}       = rmfield(source_high{isubject},'avg');
    %     source_low{isubject}.pow   = source_low{isubject}.avg.pow;
    %     source_high     = rmfield(source_low{isubject},'avg');
    %     source_high{isubject} = rmfield(source_high{isubject},'method');
    %     source_low{isubject} = rmfield(source_low{isubject},'method');
    source_high_int_rel{isubject}       = source_high_int{isubject};
    source_high_int_rel{isubject}.pow   = source_high_int{isubject}.pow ./ (source_high_int{isubject}.pow + source_high_int{isubject}.pow);
    source_low_int_rel{isubject}        = source_low_int{isubject};
    source_low_int_rel{isubject}.pow    = source_low_int{isubject}.pow ./ (source_high_int{isubject}.pow + source_high_int{isubject}.pow);

end



cfg = [];
% cfg.parameter               = 'pow';
source_diff_GA              = ft_sourcegrandaverage(cfg,source_diff{slist});


% interpolate to template MRI
cfg = [];
cfg.parameter               = {'pow'};
source_diff_GA_int          = ft_sourceinterpolate(cfg, source_diff_GA, template_mri);


cfg = [];
cfg.method                  = 'ortho';
% cfg.funcolormap             = hot(1000);
cfg.funparameter            = 'pow';
% cfg.funcolorlim = [-60 -49];
ft_sourceplot(cfg, source_diff_GA_int);


cfg = [];
cfg.method         = 'surface';
cfg.surffile       = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
% cfg.surfinflated   = 'surface_inflated_both.mat'; % Cortical sheet from canonical MNI brain
cfg.funparameter  = 'pow';
ft_sourceplot(cfg, source_diff_GA_int);
material dull



%% statistics
cfg = [];
cfg.method              = 'montecarlo';
% cfg.method              = 'analytic';

cfg.statistic           = 'ft_statfun_depsamplesT';
% cfg.statistic           = 'ft_statfun_indepsamplesT';
% cfg.clusterthreshold    = 'nonparametric';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 4000;
cfg.alpha               = 0.025;
cfg.clusteralpha        = 0.01;
cfg.tail                = 0;
cfg.design(1,:)         = [1:length(slist) 1:length(slist)];
cfg.design(2,:)         = [ones(size(slist)) ones(size(slist))*2];
cfg.uvar                = 1; 
cfg.ivar                = 2;
cfg.parameter           = 'pow';
% stat_FFT                = ft_sourcestatistics(cfg, source_low{slist},source_high{slist});
int_stat_FFT            = ft_sourcestatistics(cfg, source_low_int{slist},source_high_int{slist});
int_stat_FFT.anatomy    = template_mri.anatomy;
int_stat_FFT.stat(~int_stat_FFT.inside) = nan;

% cfg = [];
% cfg.parameter = {'stat','mask','negclusterslabelmat','prob'};
% cfg.interpmethod = 'nearest'; % for making table later
% stat_FFT_int            = ft_sourceinterpolate(cfg, stat_FFT, template_mri);
% stat_FFT_int.inside     = template_mri.inside;

mri_1mm = ft_read_mri('d:/fieldtrip-master_21092017/fieldtrip-master/template/anatomy/single_subj_T1_1mm.nii');

% interpolate to template MRI
cfg = [];
cfg.parameter               = {'stat','mask'};
int_stat_FFT_1mm          = ft_sourceinterpolate(cfg, int_stat_FFT, mri_1mm);


cfg = [];
cfg.method                  = 'ortho';
% cfg.funcolormap             = hot(1000);
cfg.funparameter            = 'stat';
cfg.maskparameter           = 'mask';
% cfg.funcolorlim = [-60 -49];
% ft_sourceplot(cfg, stat_FFT_int);
cfg.location = 'min';
cfg.crosshair = 'no';
% ft_sourceplot(cfg, int_stat_FFT);
ft_sourceplot(cfg, int_stat_FFT_1mm);


% stat_FFT_int = rmfield(stat_FFT_int,'cfg');

cfg = [];
cfg.method          = 'surface';
% cfg.surfinflated    = 'surface_inflated_both.mat'; % Cortical sheet from canonical MNI brain
cfg.surffile        = 'surface_white_both.mat'; % Cortical sheet from canonical MNI brain
cfg.surffile        = 'surface_pial_right.mat';
% cfg.opacitymap    = 'stat';
% cfg.opacitymap    = 'rampup';
cfg.maskparameter   = 'mask';
cfg.funparameter    = 'stat';
% cfg.surfdownsample  = 10;
ft_sourceplot(cfg, int_stat_FFT);
view(60,30);
material dull
saveas(gcf,'C:\Users\Admin\Dropbox\WANDER_article\source_alpha','fig');


ft_sourceplot(cfg, int_stat_FFT);


% view(0,90)


material dull
material shiny

%%


%%%%
aal = ft_read_atlas('W:\WANDER\aal\ROI_MNI_V5.nii');

template_mri = ft_read_mri('D:\fieldtrip\fieldtrip.git\trunk\external\spm8\templates\T1.nii');

% make table
clear target tabsource tabtemp
GA_int = int_stat_FFT;
stat_FFT = int_stat_FFT
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
% tabsource(tabsource.Percent_activation < 1,:) = [];
tabsource = sortrows(tabsource,{'Cluster','Size'},{'ascend','descend'})
% tabsource = sortrows(tabsource,{'cluster','stat_mean'},{'ascend','ascend'})
% tabsource = sortrows(tabsource,{'regionperc','prob_mean'},{'descend','ascend'})
% tabsource = sortrows(tabsource,{'stat_mean'},{'descend'})

tabsource(:,1:10)






