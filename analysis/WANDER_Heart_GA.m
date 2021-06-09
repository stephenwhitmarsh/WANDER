
addpath('D:\fieldtrip\fieldtrip.git\trunk');
addpath('D:/analysis/WANDER/scripts/');
ft_defaults

restingstate    = 0;
rootpath        = 1;
force           = 0;
latency         = 'all';
timing          = 'cue';

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD
clear rating_ibi_avg
slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

for isubject = slist
    [HEF_avg{isubject},HEF_avg_low{isubject},HEF_avg_high{isubject},rating_ibi_avg(isubject,:)] = WANDER_Heart(isubject,0,timing,rootpath);
    %     F = ceil(2 * tiedrank(HEF_pks_trls.trialinfo(:,2)) / length(HEF_pks_trls.trialinfo(:,2)));
    HEF_avg_diff{isubject}      = HEF_avg_low{isubject};
    HEF_avg_diff{isubject}.avg  = (HEF_avg_low{isubject}.avg - HEF_avg_high{isubject}.avg) ./ (HEF_avg_low{isubject}.avg + HEF_avg_high{isubject}.avg);
%     
%     IBI_high{isubject}  = HEF_avg_high{isubject}.trialinfo;
%     IBI_low{isubject}   = HEF_avg_low{isubject}.trialinfo;
    IBI{isubject}       = HEF_avg{isubject}.trialinfo;
%     IBI_high_avg(isubject) = nanmean(HEF_avg_high{isubject}.trialinfo);
%     IBI_low{isubject} = HEF_avg_low{isubject}.trialinfo;  
close all;
end

% add trialcount discounting blocknr
for isubject = slist
    x = abs(diff(IBI{isubject}(:,18)));
    indx = find(x > 10);
    for i = 1 : length(IBI{isubject})
        m = find(i<= indx,1,'first');
        if isempty(m) 
            m = 4;
        end
        IBI{isubject}(i,19) = (m-1)*50 + IBI{isubject}(i,18);
    end

    for itrial = 1 : 200
        trial_indx = find(IBI{isubject}(:,19) == itrial);
        if isempty(trial_indx)
            fprintf('Could not find trial %d in subject %d \n',itrial,isubject);
        else
        IBI_trial_avg{isubject}(itrial)   = nanmean(IBI{isubject}(trial_indx,17));
        trialinfo_new{isubject}(itrial,:) = IBI{isubject}(trial_indx(1),:);
        end
    end
        
    F = ceil(2 * tiedrank(IBI{isubject}(:,2)) / length(IBI{isubject}(:,2)));
    
    IBI_low_avg(isubject)  = nanmean(IBI{isubject}(F == 1,17));   
    IBI_high_avg(isubject) = nanmean(IBI{isubject}(F == 2,17));
    IBI_low_std(isubject)  = nanstd(IBI{isubject}(F == 1,17));   
    IBI_high_std(isubject) = nanstd(IBI{isubject}(F == 2,17));    
    IBI_diff_avg(isubject) = (IBI_low_avg(isubject) - IBI_high_avg(isubject)) / (IBI_low_avg(isubject) + IBI_high_avg(isubject)) ;
    IBI_diff_std(isubject) = (IBI_low_std(isubject) - IBI_high_std(isubject)) / (IBI_low_std(isubject) + IBI_high_std(isubject));
    IBI_low_count(isubject) = size(find(F == 1),1);
    IBI_high_count(isubject) = size(find(F == 2),1);

end

IBI_low_avg_avg  = mean(IBI_low_avg(slist));
IBI_low_std_avg  = mean(IBI_low_std(slist));
IBI_high_avg_avg = mean(IBI_high_avg(slist));
IBI_high_std_avg = mean(IBI_high_std(slist));
IBI_diff_avg_avg = mean(IBI_diff_avg(slist));
IBI_diff_std_avg = mean(IBI_diff_std(slist));

[h,p,ci,stats] = ttest2(IBI_high_avg(slist),IBI_low_avg(slist))
figure; errorbar([IBI_high_avg_avg, IBI_low_avg_avg],[IBI_low_std_avg, IBI_high_std_avg]);
[h,p,ci,stats] = ttest2(IBI_high_std(slist),IBI_low_std(slist))
[h,p,ci,stats] = ttest(IBI_diff_avg(slist))


for isubject = slist
    F = ceil(2 * tiedrank(trialinfo_new{isubject}(:,2)) / length(trialinfo_new{isubject}(:,2)));
    IBI_trial_avg_high(isubject) = nanmean(IBI_trial_avg{isubject}(F == 2));
    IBI_trial_avg_low(isubject) = nanmean(IBI_trial_avg{isubject}(F == 1));
    IBI_trial_avg_rel(isubject) = (IBI_trial_avg_low(isubject) - IBI_trial_avg_high(isubject)) / (IBI_trial_avg_low(isubject) + IBI_trial_avg_high(isubject));
    trialinfo_new{isubject}(:,20) = F;
    trialinfo_new{isubject}(:,21) = isubject;

end

% append data over subjects
trialinfo_concatinated = [];
for isubject = slist
    trialinfo_concatinated = [trialinfo_concatinated; trialinfo_new{isubject} ];
end

header = {'stimcode','rating','response','stimduration','


IBI_trial_avg_high_avg = mean(IBI_trial_avg_high(slist));
IBI_trial_avg_high_std = std(IBI_trial_avg_high(slist));
IBI_trial_avg_low_avg  = mean(IBI_trial_avg_low(slist));
IBI_trial_avg_low_std  = std(IBI_trial_avg_low(slist));
IBI_trial_avg_rel_avg  = mean(IBI_trial_avg_rel(slist));
IBI_trial_avg_rel_std  = std(IBI_trial_avg_rel(slist));

figure; errorbar([IBI_trial_avg_high_avg, IBI_trial_avg_low_avg],[IBI_trial_avg_high_std, IBI_trial_avg_low_std]);
figure; errorbar([IBI_trial_avg_rel_avg],[IBI_trial_avg_rel_std]); 

avg_IBI_rating = nanmean(rating_ibi_avg(slist,:));
std_IBI_rating =  nanstd(rating_ibi_avg(slist,:));
figure; errorbar(avg_IBI_rating,std_IBI_rating)

avg_rel_IBI_rating = nanmean(rating_ibi_avg(slist,:)./nanmean(rating_ibi_avg(slist,:),2),1);
std_rel_IBI_rating = nanstd(rating_ibi_avg(slist,:)./nanmean(rating_ibi_avg(slist,:),2),1);
figure; errorbar(avg_rel_IBI_rating,std_rel_IBI_rating/5);  xlim([2 7]);
figure; plot(rating_ibi_avg'./nanmean(rating_ibi_avg,2)')

HEF_avg_GA = ft_timelockgrandaverage([],HEF_avg{slist});
HEF_avg_low_GA = ft_timelockgrandaverage([],HEF_avg_low{slist});
HEF_avg_high_GA = ft_timelockgrandaverage([],HEF_avg_high{slist});

HEF_avg_diff_GA = ft_timelockgrandaverage([],HEF_avg_diff{slist});

% Plot average
cfg         = [];
cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEG*1';
% 
% cfg.layout  = 'neuromag306cmb';
% cfg.channel = 'MEG*3';

%     fig = figure; ft_singleplotER(cfg,HEF_avg);
fig = figure; ft_singleplotER(cfg,HEF_avg_low_GA,HEF_avg_high_GA);

