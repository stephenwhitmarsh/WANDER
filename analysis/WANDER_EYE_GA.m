% List of subjects without more than 2 SD of artefacts
slist = [1:5 8:13 15:20 22:26]; 

% addpath('D:\fieldtrip\fieldtrip.git\trunk');
addpath('D:/analysis/WANDER/scripts/');
addpath D:\analysis\WANDER\scripts\Inpaint_nans

ft_defaults

rootpath = 1;
[dataset_task, dataset_rs] = WANDER_subjectinfo;

for isubject = slist
    fprintf('loading subject %s \n',num2str(isubject));
    fname = ['w:\WANDER\data\EYE\s' num2str(isubject) '.mat'];
    temp = load(fname,'EYE_fix');
    EYE{isubject} = temp.EYE_fix;
    
    cfg = [];
    cfg.resamplefs = 100;
    EYE{isubject} = ft_resampledata(cfg, EYE{isubject});
end



% 
% [~, EYE_avg{isubject}, EYE_avg_low{isubject}, EYE_avg_high{isubject} ,EYE_fix{isubject} ,EYE_avg_fix{isubject}, EYE_avg_fix_low{isubject}, EYE_avg_fix_high{isubject}] = WANDER_EYE(isubject,0,1);
% 
% figure; hist(EYE_fix{isubject}.trialinfo(:,2),[1:8])
% 
% for isubject = slist
% figure;
%     
%     subplot(26,1,isubject)
%     F{isubject} = ceil(5 * tiedrank(EYE{isubject}.EYE_fix.trialinfo(:,2)) / length(EYE{isubject}.EYE_fix.trialinfo(:,2)));
%     hist(F{isubject},4);
%     
%     
% end

% Nr. 10 has a single rating of 8, and none of 7, so just add to 6.
% indx = find(EYE{10}.EYE_fix.trialinfo(:,2) == 8);
% EYE{10}.EYE_fix.trialinfo(indx,2) = 5;
% indx = find(EYE{10}.EYE_fix.trialinfo(:,2) == 6);
% EYE{10}.EYE_fix.trialinfo(indx,2) = 5; 
    
nrbins = 4;
i = 1;
for isubject = slist
    [ratings_binned{isubject}, ratings_old{isubject}, ratings_stat{isubject}] = bin_data(EYE{isubject}.EYE_fix.trialinfo(:,2),nrbins);
    i = i + 1;
end

% normalize data over trials per subject
for isubject = slist
    M = size(EYE{isubject}.EYE_fix.trial,2);
    N = size(EYE{isubject}.EYE_fix.trial{1},2);
    EYE_norm{isubject} = EYE{isubject}.EYE_fix; 
    temp1 = zscore(cell2mat(EYE{isubject}.EYE_fix.trial));
    temp2 = reshape(temp1,N,M)';
    EYE_norm{isubject}.trial = mat2cell(temp2,ones(1,M),N)';
end

figure; hold;
c = hsv(nrbins);
c = c(nrbins:-1:1,:);
for irating = 1:nrbins

    % average over binned ratings
    for isubject = slist
        cfg = [];
        cfg.trials = find(ratings_binned{isubject} == irating);
        cfg.avgoverrpt                          = 'yes';
        rating_avg{irating}{isubject}           = ft_selectdata(cfg,EYE_norm{isubject});
        rating_avg{irating}{isubject}.avg       = rating_avg{irating}{isubject}.trial{1};
        rating_avg{irating}{isubject}.time      = rating_avg{irating}{isubject}.time{1};
        rating_avg{irating}{isubject}.dimord    = 'time';
        rating_avg{irating}{isubject}           = rmfield(rating_avg{irating}{isubject},'trial');

    end
    rating_GA{irating} = ft_timelockgrandaverage([],rating_avg{irating}{slist});
    
    % manual calculate std
    i = 1;
    clear temp
    for isubject = slist
        temp(i,:) = rating_avg{irating}{isubject}.avg;
        i = i + 1;
    end
    rating_GA{irating}.std = nanstd(temp,1) / sqrt(22);
           
%     fill([rating_GA{irating}.time rating_GA{irating}.time(end:-1:1)],[rating_GA{irating}.avg+rating_GA{irating}.std, rating_GA{irating}.avg(end:-1:1)-rating_GA{irating}.std(end:-1:1)],c(irating,:),'linestyle','none');
    %     plot(rating_GA{irating}.avg,'Color',c(irating,:));
end

figure; hold;
for irating = 1 : nrbins
    fill([rating_GA{irating}.time rating_GA{irating}.time(end:-1:1)],[rating_GA{irating}.avg+rating_GA{irating}.std, rating_GA{irating}.avg(end:-1:1)-rating_GA{irating}.std(end:-1:1)],c(irating,:),'linestyle','none');
end
for irating = 1 : nrbins
    plot(rating_GA{irating}.time,rating_GA{irating}.avg,'Color',c(irating,:),'Linewidth',1);
end

% set(gcf, 'renderer', 'painters', 'paperpositionmode', 'auto');
print('-dsvg','w:\WANDER\images\pupil_4ratings.svg');

% 
% figure; hold;
% c = jet(8);
% for irating = 1:8
%     clear t
%     i = 1;
%     for isubject = slist
%         if rating_present(isubject,irating)
%             t{i} = EYE_avg_fix_rating{isubject}{irating};
% %             t{i}.avg = t{i}.avg ./ EYE_avg_fix{isubject}.avg;
%             i = i + 1;
%         end
%     end
%     avg{irating} = ft_timelockgrandaverage(t{:});
%     plot(avg{irating}.avg,'Color',c(irating,:));
% 
% end

 EYE_avg_fix_rating{isubject}{irating}.avg ./ EYE_avg_fix{isubject}.avg;



figure; hold;
for irating = 1:8
    



EYE_avg_fix_rating{rating_present(:,1)}{1}

{1}


for isubject = slist
    
 

    EYE_avg_fix_rel{isubject} = EYE_avg_fix_high{isubject};
    EYE_avg_fix_rel{isubject}.avg = (EYE_avg_fix_low{isubject}.avg - EYE_avg_fix_high{isubject}.avg) ./ (EYE_avg_fix_low{isubject}.avg + EYE_avg_fix_high{isubject}.avg); 
%     
% %     cfg = [];
% %     cfg.lpfilter = 'yes';
%     cfg.lpfreq = 30;
%     EYE_avg_fix_high{isubject} = EYE_avg{isubject};
%     
%     
end

figure;
i = 1;
for isubject = slist
    
%         EYE_avg_fix_low{isubject}.avg = EYE_avg_fix_low{isubject}.avg';
%         EYE_avg_fix_high{isubject}.avg = EYE_avg_fix_high{isubject}.avg';
%         EYE_avg_fix_rel{isubject}.avg = EYE_avg_fix_rel{isubject}.avg';
%     
    cfg = [];
    subplot(5,5,i);
    cfg.channel = 'MISC009';
    ft_singleplotER(cfg,EYE_avg_fix_low{isubject},EYE_avg_fix_high{isubject});
%         ft_singleplotER(cfg,EYE_avg_fix_rel{isubject});

    title(num2str(isubject));
    i = i + 1;
end

cfg = [];
cfg.channel = 'MISC009';
EYE_GA          = ft_timelockgrandaverage(cfg,EYE_avg{slist});
EYE_high_GA     = ft_timelockgrandaverage(cfg,EYE_avg_high{slist});
EYE_low_GA      = ft_timelockgrandaverage(cfg,EYE_avg_low{slist});

EYE_fix_GA      = ft_timelockgrandaverage(cfg,EYE_avg_fix{slist});
EYE_fix_high_GA = ft_timelockgrandaverage(cfg,EYE_avg_fix_high{slist});
EYE_fix_low_GA  = ft_timelockgrandaverage(cfg,EYE_avg_fix_low{slist});
EYE_fix_rel_GA  = ft_timelockgrandaverage(cfg,EYE_avg_fix_rel{slist});


figure;
ft_singleplotER(cfg,EYE_fix_low_GA,EYE_fix_high_GA);


Nsub                    = length(slist);
cfg                     = [];
cfg.channel             = 'MISC009';
% cfg.latency             = [10 11];
cfg.avgovertime         = 'no';
% cfg.avgoverchan         = 'no';
cfg.parameter           = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.clusteralpha        = 0.05;
cfg.clusterstatistic    = 'maxsum';
% cfg.minnbchan           = 1;
cfg.tail                = 1;
cfg.clustertail         = 1;
cfg.alpha               = 0.05;
cfg.numrandomization    = 2000;
% cfg.neighbours          = neigh_cmb.neighbours; 
% cfg.neighbours          = neigh_mag.neighbours; 
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; 
cfg.uvar                = 2;
stat_EYE                 = ft_timelockstatistics(cfg,EYE_avg_fix_low{slist},EYE_avg_fix_high{slist});

figure; plot(stat_EYE.time,stat_EYE.stat);
hold
plot(stat_EYE.time,stat_EYE.mask);

stat_EYE.prob(stat_EYE.posclusterslabelmat==1)
mean(stat_EYE.stat(stat_EYE.posclusterslabelmat==1))

t = stat_EYE.time(stat_EYE.posclusterslabelmat==1)
t(1)-10
t(end)-10

figure; plot(stat_EYE.time-10,[EYE_fix_low_GA.avg; EYE_fix_high_GA.avg]); hold;
plot(stat_EYE.time-10,stat_EYE.mask);

figure; plot(EYE_fix_rel_GA.time,EYE_fix_rel_GA.avg);

% addpath('D:\analysis\WANDER\scripts\plot2svg-master\src');

fig = figure; hold;
t = 2.086;
t = 1;
% patch([EYE_fix_GA.time-10 EYE_fix_GA.time(end:-1:1)-10],[EYE_fix_low_GA.avg+sqrt(EYE_fix_low_GA.var)   ./sqrt(22)*t, EYE_fix_low_GA.avg(end:-1:1)-sqrt(EYE_fix_low_GA.var(end:-1:1)) ./sqrt(22)*t],[0 0.8 0 ],'FaceAlpha',.5,'LineStyle','None');
% patch([EYE_fix_GA.time-10 EYE_fix_GA.time(end:-1:1)-10],[EYE_fix_high_GA.avg+sqrt(EYE_fix_high_GA.var) ./sqrt(22)*t, EYE_fix_high_GA.avg(end:-1:1)-sqrt(EYE_fix_high_GA.var(end:-1:1)) ./sqrt(22)*t],[0.8 0 0 ],'FaceAlpha',.5,'LineStyle','None');


fill([EYE_fix_GA.time-10 EYE_fix_GA.time(end:-1:1)-10],[EYE_fix_low_GA.avg+sqrt(EYE_fix_low_GA.var)   ./sqrt(22)*t, EYE_fix_low_GA.avg(end:-1:1)-sqrt(EYE_fix_low_GA.var(end:-1:1))   ./sqrt(22)*t],[1 1 0 ],'linestyle','none');
fill([EYE_fix_GA.time-10 EYE_fix_GA.time(end:-1:1)-10],[EYE_fix_high_GA.avg+sqrt(EYE_fix_high_GA.var) ./sqrt(22)*t, EYE_fix_high_GA.avg(end:-1:1)-sqrt(EYE_fix_high_GA.var(end:-1:1)) ./sqrt(22)*t],[1 0 1 ],'linestyle','none');
line(EYE_fix_GA.time-10,EYE_fix_low_GA.avg,'color','g');
line(EYE_fix_GA.time-10,EYE_fix_high_GA.avg,'color','r');
axis tight;

set(gca,'PlotBoxAspectRatio',[8/3,1,1]);
savefig(fig,'w:\WANDER\images\FIGb_pupil');

set(gcf, 'renderer', 'painters', 'paperpositionmode', 'auto');
print('-dpdf','w:\WANDER\images\FIGb_pupil.pdf');
% saveas(fig,'w:\WANDER\images\FIGb_pupil.eps','eps');


axis tight;


%% plot percentage change
eye_avg_diff = (mean(abs(EYE_fix_low_GA.avg(1,stat_EYE.mask))) - mean(abs(EYE_fix_high_GA.avg(1,stat_EYE.mask)))) / (mean(abs(EYE_fix_low_GA.avg(1,stat_EYE.mask))) + mean(abs(EYE_fix_high_GA.avg(1,stat_EYE.mask))));



% 
% 
% 
% cfg = [];
% cfg.channel = 'MISC009';
% figure; ft_singleplotER(cfg,EYE_GA);
% figure; ft_singleplotER(cfg,EYE_fix_rel_GA);
% 
% figure; ft_singleplotER(cfg,EYE_high_GA,EYE_low_GA);
% 
% cfg = [];
% cfg.channel = 'MISC009';
% figure; ft_singleplotER(cfg,EYE_fix_GA);
% figure; ft_singleplotER(cfg,EYE_fix_high_GA,EYE_fix_low_GA);
% 
% 

% 
% 
% % maximum length trial
% maxlength = 0;
% for itrial = 1 : size(data_epoch_MEG.trial,2)
%     if size(data_epoch_MEG.trial{itrial},2) > maxlength
%         maxlength = size(data_epoch_MEG.trial{itrial},2);
%     end
% end
% 
% % put data in matrix
% dat = NaN(size(data_epoch_MEG.trial,2),size(data_epoch_MEG.label,1),maxlength);
% for itrial = 1 : size(data_epoch_MEG.trial,2)
%     fprintf('Adding trial %d to data matrix \n',itrial);
%     dat(itrial,:,1:size(data_epoch_MEG.trial{itrial},2)) = data_epoch_MEG.trial{itrial};
% end
% 
% 
% 
% 
% 
% % manually average over all trials, excluding nans
% cfg = [];
% cfg.vartrllength    = 2;
% ERF_cue             = ft_timelockanalysis(cfg,data_epoch_MEG);
% ERF_cue.avg         = squeeze(nanmean(dat,1));
% ERF_cue             = rmfield(ERF_cue,'var');
% 
% % make index according to median split based on ratings
% F = ceil(2 * tiedrank(data_epoch_MEG.trialinfo(:,2)) / length(data_epoch_MEG.trialinfo(:,2)));
% 
% % before equalizing high/low bins
% for irand = 1 : 24
%     
%     rating_split        = ones(size(F));
%     rating_split(F==2)  = 2;
%     low_cnt             = size(find(rating_split == 1),1);
%     high_cnt            = size(find(rating_split == 2),1);
%     
%     rating_split       = ones(size(F));
%     rating_split(F==2) = 2;
%     if high_cnt > low_cnt
%         diff_bar                = min(data_epoch_MEG.trialinfo(rating_split == 2,2)); % rating bin between high/low split
%         nr_to_remove            = high_cnt-low_cnt;
%         discarted               = randsample(find(data_epoch_MEG.trialinfo(:,2) == diff_bar),nr_to_remove);
%         rating_split(discarted) = 0;
%     elseif high_cnt < low_cnt
%         diff_bar                = max(data_epoch_MEG.trialinfo(rating_split == 1,2)); % rating bin between high/low split
%         nr_to_remove            = low_cnt-high_cnt;
%         discarted               = randsample(find(data_epoch_MEG.trialinfo(:,2) == diff_bar),nr_to_remove);
%         rating_split(discarted) = 0;
%     else
%         discarted = [];
%     end
%     fprintf('I made %d low ratings, and %d high ratings in randomization %d \n',sum(rating_split==1),sum(rating_split==2),irand);
%     
%     temp_low(irand,:,:)  = squeeze(nanmean(dat(rating_split==1,:,:),1));
%     temp_high(irand,:,:) = squeeze(nanmean(dat(rating_split==2,:,:),1));
%     
%     % average over randomizations
%     ERF_cue_high        = ERF_cue;
%     ERF_cue_high.avg    = squeeze(nanmean(temp_high,1));
%     ERF_cue_low         = ERF_cue;
%     ERF_cue_low.avg     = squeeze(nanmean(temp_low,1));
%     clear temp_low temp_high
% end
% 
% for irating = 1 : 7
%     rating_cnt_low(irating)       = size(find(data_epoch_MEG.trialinfo(:,2) == irating & rating_split == 1),1);
%     rating_cnt_high(irating)      = size(find(data_epoch_MEG.trialinfo(:,2) == irating & rating_split == 2),1);
%     rating_cnt_discarted(irating) = size(find(data_epoch_MEG.trialinfo(discarted,2) == irating),1);
% end
