function WANDER_behaviour_average

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

stimduration         = nan(26,200);
hits                 = nan(26,1);
correct_rejections   = nan(26,1);
misses               = nan(26,1);
FAs                  = nan(26,1);
rating               = nan(26,200);
stimduration_mean    = nan(26,8);
stimduration_correct = nan(26,200);

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

for isubject = slist
    behav                                                       = WANDER_behaviour(isubject,0,0);
    stimduration(isubject,1:length(behav.stimduration))   = behav.stimduration;
    hits(isubject)                                              = behav.hits;
    correct_rejections(isubject)                                = behav.correct_rejections;
    misses(isubject)                                            = behav.misses;
    FAs(isubject)                                               = behav.FAs;
    rating(isubject,1:length(behav.rating))              = behav.rating;
    stimduration_mean(isubject,:)                               = behav.stimduration_mean;
    trialinfo{isubject}                                         = behav.trialinfo;
end
trialinfo{17} = [trialinfo{17}; nan; nan]; % subject misses first last trials of first block
rating(17,3:end) = rating(17,1:end-2); rating(17,1:2) = [nan nan];

rating_count = nan(26,8);
for isubject = slist
    for irating = 1 : 8
        rating_count(isubject,irating) = sum(rating(isubject,:) == irating);
    end
end

% overview number and trialduration ratings
rating_count_avg = nanmean(rating_count,1);
rating_count_std = nanstd(rating_count,1);
stimduration_mean_avg = nanmean(stimduration_mean,1);
stimduration_mean_std = nanstd(stimduration_mean,1);
fig = figure;
subplot(3,1,1); errorbar(rating_count_avg,rating_count_std); hold; bar(rating_count_avg); ylim([0 75]);
subplot(3,1,2); errorbar(stimduration_mean_avg,stimduration_mean_std); hold; bar(stimduration_mean_avg); ylim([10 20]);

% overview rating over time (per ten trials)
clear rating_time_avg rating_time_std
i = 1;
stepsize = 1;
for itime = 1 : stepsize : 200
    rating_time_avg(i) = nanmean(nanmean(rating(slist,[itime : itime + stepsize-1]),2));
    rating_time_std(i) = nanstd(nanmean(rating(slist,[itime : itime + stepsize-1]),2));
    i = i + 1;
end
subplot(3,1,3); errorbar(rating_time_avg,rating_time_std); hold; bar(rating_time_avg,1); axis tight; ylim([4 7]);

print(fig,'-deps',['d:\analysis\WANDER\images\behaviour\overview_ratings.eps']);


% get triallengths for correct responses
for isubject = slist
    for irating = 1 : 8
        stimduration_mean(isubject,irating) = nanmean(stimduration(isubject,(rating(isubject,:) == irating)' & trialinfo{isubject} == 4));
    end
end

% some ANOVAs
[p,tbl,stats]               = anova1(stimduration_mean(slist,:));
figure;
[comparison,means,h,gnames] = multcompare(stats,'display','on');


[p,tbl,stats]               = anova1(stimduration(slist,:));
figure;
[comparison,means,h,gnames] = multcompare(stats,'display','on');

trialnr = [1:200];
mean_ratings = nanmean(rating(slist,:));

% [p,S,mu]  = polyfit(mean_ratings,trialnr,1)
% Y = polyconf(p,mean_ratings)
% [y,delta] = polyval(p,x,S,mu)

% linear fit over all trials
mdl = fitlm(mean_ratings,trialnr)
mdl = fitlm(trialnr,mean_ratings)

% linear fit over all trials over blocks
figure; plot(mean([mean_ratings(1:50); mean_ratings(51:100); mean_ratings(101:150); mean_ratings(151:200)]));
mdl = fitlm(mean([mean_ratings(1:50); mean_ratings(51:100); mean_ratings(101:150); mean_ratings(151:200)]),1:50);
% mdl = fitlm(1:50,mean([mean_ratings(1:50); mean_ratings(51:100); mean_ratings(101:150); mean_ratings(151:200)]))


stimduration_mean_avg = nanmean(stimduration(slist,:),2);
stimduration_mean_std = nanstd(stimduration(slist,:),1);

stimduration_mean_avg = nanmean(stimduration_mean(slist,:),1);
stimduration_mean_std = nanstd(stimduration_mean(slist,:),1);

figure; hold
boxplot(stimduration_mean);


plot([1:8],stimduration_mean_avg+stimduration_mean_std);
plot([1:8],stimduration_mean_avg-stimduration_mean_std);
plot([1:8],stimduration_mean_avg);
%
% stimduration_avg = nanmean(stimduration');
% stimduration_max = max(stimduration');
% stimduration_min = min(stimduration');

% figure; hold
% plot([1:26],stimduration_avg);
% plot([1:26],stimduration_max);
% plot([1:26],stimduration_min);

% plot performance of all subjects
fig = figure;
for isubject = 1 : 26
    subplot(5,6,isubject);
    ax  = gca;
    x   = [1 2 3 4];
    y   = [hits(isubject);correct_rejections(isubject);misses(isubject);FAs(isubject)];
    bar(y,'FaceColor',[0 0 0]);
    ax.XTickLabel = {'H','C','M','F'};
    for i1=1:4
        text(x(i1),y(i1),num2str(y(i1),'%d'),'HorizontalAlignment','center','VerticalAlignment','bottom')
    end
    xlim([0.5 4.5]);
    ylim([0 220]);
    title(isubject);
end

% plot mean over subjects
% fig = figure;
subplot(5,6,30);

ax  = gca;
hitrate             = (hits               ./ (hits + misses));
nanmean(hitrate(slist))
nanstd(hitrate(slist))

missrate            = (misses             ./ (hits + misses));
rejectionrate       = (correct_rejections ./ (correct_rejections + FAs));
nanmean(rejectionrate(slist))
nanstd(rejectionrate(slist))



FArate              = (FAs                ./ (correct_rejections + FAs));

x   = [1 2 3 4];
y   = [nanmean(hitrate);nanmean(rejectionrate);nanmean(missrate);nanmean(FArate)]*100;
sd  = [nanstd(hitrate);nanstd(rejectionrate);nanstd(missrate);nanstd(FArate)]*100;
% errorbar(y,sd); hold;
bar(y,'FaceColor',[0 0 0]);

ax.XTickLabel = {'H','C','M','F'};
ylabel('nr of trials');
for i1=1:4
    text(x(i1),y(i1),sprintf('%0.1f',y(i1)),'HorizontalAlignment','center','VerticalAlignment','bottom')
end
xlim([0.5 4.5]);
ylim([0 110]);
title('mean');

savefig(fig,'/shared/projects/project_wander/WANDER\images\poster\hitrates');


%
% set(fig,'PaperSize',[11*2 8.5*2 ]);
% set(fig,'PaperPosition',[0 0  11*2 8.5*2]);
% set(fig,'PaperOrientation','landscape');
% set(fig,'Position',[50 50 1200 800]);
% print(fig,'-dpdf',['d:\analysis\WANDER\images\behaviour\avg_responses.pdf']);
% print(fig,'-dpng',['d:\analysis\WANDER\images\behaviour\avg_responses.png']);

% stimdurating distribution
stimduration_resh = reshape(stimduration,1,200*26);
min(stimduration_resh)
max(stimduration_resh)



fig = figure;
hist(stimduration_resh,30);
title('Stimulation duration');
xlabel('time (s)');
ylabel('nr of trials');
set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0  11*2 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['/shared/projects/project_wander/WANDER\images\behaviour\avg_duration_dist.pdf']);
print(fig,'-dpng',['/shared/projects/project_wander/WANDER\images\behaviour\avg_duration_dist.png']);

% rating numbers for correct rejections
fig = figure;
iplot = 1;
for isubject = slist
    subplot(5,5,iplot);
    histogram(rating(isubject,trialinfo{isubject} == 4),0.5:8.5,'Normalization','count','EdgeColor',[0 0 0],'FaceColor',[0 0 0]); axis tight;
    title(isubject);
    ylim([0,120]);
    iplot = iplot + 1;
end
subplot(5,5,iplot);
histogram(rating,0.5:8.5,'Normalization','count','EdgeColor',[0 0 0],'FaceColor',[0.5 0.5 0.5]); axis tight;
title('mean');
ylabel('nr of trials');
xlabel('rating');
set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0  11*2 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['/shared/projects/project_wander/WANDER\images\behaviour\avg_rating_count.pdf']);
print(fig,'-dpng',['/shared/projects/project_wander/WANDER\images\behaviour\avg_rating_count.png']);

% plot nr of ratings with median split
fig = figure;
clear rating_split
iplot = 1;
for isubject = slist
    subplot(5,5,iplot);
    
    % make index according to median split based on ratings
    F = ceil(2 * tiedrank(rating(isubject,:)) / length(rating(isubject,:)));
    
    % before equalizing high/low bins
    rating_split(isubject,:)     = ones(size(F));
    rating_split(isubject,F==2)  = 2;
    
%     histogram(rating(isubject,rating_split(isubject,:) == 1 & trialinfo{isubject}' == 4),0.5:8.5,'Normalization','count','EdgeColor',[0 0 0],'FaceColor',[0 1 0]); hold;
%     histogram(rating(isubject,rating_split(isubject,:) == 2 & trialinfo{isubject}' == 4),0.5:8.5,'Normalization','count','EdgeColor',[0 0 0],'FaceColor',[1 0 0]); axis tight;
    histogram(rating(isubject,rating_split(isubject,:) == 1),0.5:8.5,'Normalization','count','EdgeColor',[0 0 0],'FaceColor',[0 1 0]); hold;
    histogram(rating(isubject,rating_split(isubject,:) == 2),0.5:8.5,'Normalization','count','EdgeColor',[0 0 0],'FaceColor',[1 0 0]); axis tight;    
    title(isubject);
    ylim([0,120]);
    iplot = iplot + 1;
    
end


title('mean');
ylabel('nr of trials');
xlabel('rating');
set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0  11*2 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['/shared/projects/project_wander/WANDER\images\behaviour\avg_rating_count_mediansplit.pdf']);
print(fig,'-dpng',['/shared/projects/project_wander/WANDER\images\behaviour\avg_rating_count_medianSplit.png']);

% duration per rating
fig = figure;
for isubject = 1 : 26
    subplot(5,6,isubject);
    bar(stimduration_mean(isubject,:),'FaceColor',[0 0 0]); axis tight;
    title(isubject);
    ylim([0,20]);
end

subplot(5,6,30);
bar(nanmean(stimduration_mean),'FaceColor',[0.5 0.5 0.5]); axis tight;
title('mean');
ylim([0,20]);
ylabel('trialduration');
xlabel('rating');

set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0  11*2 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['/shared/projects/project_wander/WANDER\images\behaviour\avg_rating_duration.pdf']);
print(fig,'-dpng',['/shared/projects/project_wander/WANDER\images\behaviour\avg_rating_duration.png']);

% ratings over time/block

rating_cr = rating;
for isubject = 1 : 26
    rating_cr(isubject,trialinfo{isubject} ~= 4) = NaN;
end
avg_rating_cr = nanmean(rating_cr,1);

rating_mistake = rating;
for isubject = 1 : 26
    rating_mistake(isubject,trialinfo{isubject} == 4) = NaN;
end
rating_mistake = nanmean(rating_mistake,1);


fig = figure;
subplot(4,4,1);
histogram(rating(:,1:50),'Normalization','count','BinMethod','auto','EdgeColor',[0 0 0],'FaceColor',[0 0 0]); axis tight;
title(sprintf('Block 1, M: %0.1f',mean(avg_rating_cr(1:50))));
ylim([0,300]);
ylabel('nr of trials');
xlabel('rating');

subplot(4,4,2);
histogram(rating(:,50:100),'Normalization','count','BinMethod','auto','EdgeColor',[0 0 0],'FaceColor',[0 0 0]); axis tight;
title(sprintf('Block 2, M: %0.1f',mean(avg_rating_cr(51:100))));
ylim([0,300]);
ylabel('nr of trials');
xlabel('rating');

subplot(4,4,3);
histogram(rating(:,100:150),'Normalization','count','BinMethod','auto','EdgeColor',[0 0 0],'FaceColor',[0 0 0]); axis tight;
title(sprintf('Block 3, M: %0.1f',mean(avg_rating_cr(101:150))));
ylim([0,300]);
ylabel('nr of trials');
xlabel('rating');

subplot(4,4,4);
histogram(rating(:,150:200),'Normalization','count','BinMethod','auto','EdgeColor',[0 0 0],'FaceColor',[0 0 0]); axis tight;
title(sprintf('Block 4, M: %0.1f',mean(avg_rating_cr(151:200))));
ylim([0,300]);
ylabel('nr of trials');
xlabel('rating');

subplot(4,4,5);
bar(nanmean(rating_cr(:,1:50)),1); axis tight;
ylabel('average rating');
xlabel('trial');
ylim([4,8]);

subplot(4,4,6);
bar(nanmean(rating_cr(:,51:100)),1); axis tight;
ylabel('average rating');
xlabel('trial');
ylim([4,8]);

subplot(4,4,7);
bar(nanmean(rating_cr(:,101:150)),1); axis tight;
ylabel('average rating');
xlabel('trial');
ylim([4,8]);

subplot(4,4,8);
bar(nanmean(rating_cr(:,151:200)),1); axis tight;
ylabel('average rating');
xlabel('trial');
ylim([4,8]);

set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0  11*2 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['/shared/projects/project_wander/WANDER\images\behaviour\avg_rating_block_distribution.pdf']);
print(fig,'-dpng',['/shared/projects/project_wander/WANDER\images\behaviour\avg_rating_block_distribution.png']);





