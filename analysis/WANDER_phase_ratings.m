function WANDER_phase_ratings

addpath d:/fieldtrip/
addpath('D:/analysis/WANDER/scripts/');
addpath('D:/analysis/WANDER/scripts/CircStat2012a/');
addpath('D:/analysis/WANDER/scripts/PhaseOppositionCode/');
ft_defaults

%%

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD
rootpath = 1;
timing = 'cue';

clear rating_* phase_end phase_begin phase_rand rating
clear IBI_trial_avg IBI_trial_onset IBI_trial_offset trialnr
[dataset_task, ~] = WANDER_subjectinfo;

for isubject = slist
    
    % load continuous data
    data_EGG = WANDER_filter_EGG(isubject,0,rootpath,0);
    
    for ipart = 1 : 4
        
        % create trial definition
        cfg = [];
        cfg.dataset                 = dataset_task{isubject,ipart};
        cfg.trialfun                = 'WANDER_trialfun_cuelocked';
        cfg.trialdef.stim_pre       = 0;            % time before onset stimulation - REMOVED FIRST SECOND DUE TO ONSET RESPONSE
        cfg.trialdef.stim_post      = 0;            % time after offset stimulation
        cfg.trialdef.onset          = [129 131];    % trigger for start stimulation
        cfg.trialdef.offset         = 3;            % trigger for end trial
        cfg.trialdef.stim           = [128 130];    % trigger for stimulations
        cfg.trialdef.rating         = 10:19;
        cfg.trialdef.performance    = 20:29;
        cfg                         = ft_definetrial(cfg);
        
%         epoch EGG data
        data_EGG_epoched{isubject}{ipart}     = ft_redefinetrial(cfg,data_EGG{ipart});
        data_EGG_epoched{isubject}{ipart}.trialinfo(:,7) = data_EGG_epoched{isubject}{ipart}.sampleinfo(:,1) - data_EGG_epoched{isubject}{ipart}.sampleinfo(1,1); % 0 at onset of each block
        data_EGG_epoched{isubject}{ipart}.trialinfo(:,8) = data_EGG_epoched{isubject}{ipart}.sampleinfo(:,2) - data_EGG_epoched{isubject}{ipart}.sampleinfo(1,1); 
        
    end
    
    data_EGG_comb{isubject} = ft_appenddata([],data_EGG_epoched{isubject}{:});
    
    clear data_EGG_epoched
    
    cfg = [];
    cfg.channel = {'phase','amplitude'};
    data_EGG_comb{isubject} = ft_selectdata(cfg,data_EGG_comb{isubject});
    
    % load heart beat data
    [HEF_avg{isubject},HEF_avg_low{isubject},HEF_avg_high{isubject},rating_ibi_avg(isubject,:)] = WANDER_Heart(isubject,0,timing,rootpath);
    IBI{isubject}       = HEF_avg{isubject}.trialinfo;
    
    % add trialcount discounting blocknr
    x = abs(diff(IBI{isubject}(:,18)));
    indx = find(x > 10);
    for i = 1 : length(IBI{isubject})
        m = find(i<= indx,1,'first');
        if isempty(m)
            m = 4;
        end
        IBI{isubject}(i,19) = (m-1)*50 + IBI{isubject}(i,18);
        IBI{isubject}(i,20) = m;
    end
    
    i = 1;
    for itrial = 1 : 200
        trial_indx = find(IBI{isubject}(:,19) == itrial);
        if isempty(trial_indx)
            fprintf('Could not find trial %d in subject %d \n',itrial,isubject);
        else
            IBI_trial_avg{isubject}(i)    = nanmean(IBI{isubject}(trial_indx,17));
            IBI_trial_onset{isubject}(i)  = IBI{isubject}(trial_indx(1),17);
            IBI_trial_offset{isubject}(i) = IBI{isubject}(trial_indx(end),17);
            trialnr{isubject}(i)          = IBI{isubject}(trial_indx(end),19);
            blocknr{isubject}(i)          = IBI{isubject}(trial_indx(end),20);
            
            i = i + 1;
        end
    end
    
    F = ceil(2 * tiedrank(IBI{isubject}(:,2)) / length(IBI{isubject}(:,2)));
    IBI_low_avg(isubject)       = nanmean(IBI{isubject}(F == 1,17));
    IBI_high_avg(isubject)      = nanmean(IBI{isubject}(F == 2,17));
    IBI_low_std(isubject)       = nanstd(IBI{isubject}(F == 1,17));
    IBI_high_std(isubject)      = nanstd(IBI{isubject}(F == 2,17));
    IBI_diff_avg(isubject)      = (IBI_low_avg(isubject) - IBI_high_avg(isubject)) / (IBI_low_avg(isubject) + IBI_high_avg(isubject)) ;
    IBI_diff_std(isubject)      = (IBI_low_std(isubject) - IBI_high_std(isubject)) / (IBI_low_std(isubject) + IBI_high_std(isubject));
    IBI_low_count(isubject)     = size(find(F == 1),1);
    IBI_high_count(isubject)    = size(find(F == 2),1);
end

%%

clear phase_end phase_begin rating correct
for isubject = slist
    for itrial = 1 : size(data_EGG_comb{isubject}.trial,2)
        phase_end{isubject}(itrial)     = data_EGG_comb{isubject}.trial{itrial}(1,end);
        phase_begin{isubject}(itrial)   = data_EGG_comb{isubject}.trial{itrial}(1,1);
        amp_end{isubject}(itrial)       = data_EGG_comb{isubject}.trial{itrial}(2,end);
        amp_begin{isubject}(itrial)     = data_EGG_comb{isubject}.trial{itrial}(2,1);
        amp_avg{isubject}(itrial)       = nanmean(data_EGG_comb{isubject}.trial{itrial}(2,:));
        rating{isubject}(itrial)        = data_EGG_comb{isubject}.trialinfo(itrial,2);
    end
end

for isubject = slist
    for itrial = 1 : size(data_EGG_comb{isubject}.trial,2)
        rating_split{isubject} = ceil(2 * tiedrank(rating{isubject}) / length(rating{isubject}));
    end
end

% combine trialinfo data
% 1)trialtype, 2)rating, 3)correct, 4)stimduration, 5)RT rating confirmation (onset =
% t0), 6)trialnr, 7)samplenr onset, 8)samplenr offset, 9)EGG phase onset,
% 10) EGG phase offset, 11) IBI first HB, 12) IBI last HB, 13) IBI avg over
% trial, 14) trialnr over blocks, 15) EGG amplitude onset, 16) EGG
% amplitude offset


for isubject = slist
    behav{isubject}       = data_EGG_comb{isubject}.trialinfo;
    behav{isubject}(:,9)  = phase_begin{isubject};
    behav{isubject}(:,10) = phase_end{isubject};
    behav{isubject}(:,11) = IBI_trial_onset{isubject};
    behav{isubject}(:,12) = IBI_trial_offset{isubject};
    behav{isubject}(:,13) = IBI_trial_avg{isubject};
    behav{isubject}(:,14) = trialnr{isubject};
    behav{isubject}(:,15) = amp_begin{isubject};
    behav{isubject}(:,16) = amp_end{isubject};
    behav{isubject}(:,17) = amp_avg{isubject};
    behav{isubject}(:,18) = isubject;
    behav{isubject}(:,19) = rating_split{isubject};
    behav{isubject}(:,20) = blocknr{isubject};
end


%% concatinate over subjects
behaviour = [];
for isubject = slist
    behaviour = [behaviour;  behav{isubject}];
end
save('D:\analysis\WANDER\data\behaviour\data_behaviour_all.mat','behav','behaviour')
dlmwrite('D:\analysis\WANDER\data\behaviour\data_behaviour_all.txt',behaviour,'precision','%.5f','roffset',1)

%% linear models

tbl = table(behaviour(:,2),behaviour(:,3),behaviour(:,7),behaviour(:,18),behaviour(:,19),'VariableNames',{'Rating','correct','Samplenr','Subjectnr','Blocknr'});

tbl = tbl(tbl.correct == 4,:);

lme = fitlme(tbl,'Rating~Samplenr+(1|Blocknr)+(Samplenr-1|Blocknr)');

lme = fitlme(tbl,'Rating~Samplenr+(Samplenr|Blocknr)');
lme = fitlme(tbl,'Rating~Samplenr+(Samplenr|Subjectnr)');

lme = fitlme(tbl,'Rating~Samplenr+(Samplenr|Blocknr)+(Samplenr|Subjectnr)');
lme = fitlme(tbl,'Rating~Samplenr+(Samplenr|Subjectnr)+(Samplenr|Subjectnr:Blocknr)');

lme = fitlme(tbl,'Rating~Samplenr+(Samplenr|Subjectnr)+(Samplenr|Subjectnr:Blocknr)');

stats = anova(lme)


%%

% plot relationship angle and rating
i = 1;
figure;
for isubject = slist
    subplot(5,5,i);
    scatter(behav{isubject}(:,10),behav{isubject}(:,2),'.');
    i = i + 1;
end

i = 1;
figure;
for isubject = slist
    subplot(5,5,i);
    correct_indx = behav{isubject}(:,3) == 4;
    only_correct = behav{isubject}(correct_indx,:);
    
    [n,bin] = histc(only_correct(:,10),-pi:2*pi/18:pi);
    bin_rating = nan(1,18);
    for ibin = 1 : 18
        bin_rating(ibin) = mean(only_correct(bin==ibin,2));
    end
    bar([-pi:2*pi/18:pi-pi/18]+pi/18/2,bin_rating);
    axis tight;
    i = i + 1;
end

i = 1;
figure;
for isubject = slist
    phase_nrbins = 9;
    subplot(5,5,i);
    correct_indx = behav{isubject}(:,3) == 4;
    only_correct = behav{isubject}(correct_indx,:);
    
    trial_count_low(isubject)  = sum(only_correct(:,19) == 1);
    trial_count_high(isubject) = sum(only_correct(:,19) == 2);
    
    indx_high  = find(only_correct(:,19) == 2);
    indx_low   = find(only_correct(:,19) == 1);
    
    mincount(isubject) = min(trial_count_low(isubject),trial_count_high(isubject));

    r = indx_high(randperm(size(indx_high,1)));
    high = only_correct(r(1:mincount(isubject)),:);
    
    r = indx_low(randperm(size(indx_low,1)));
    low = only_correct(r(1:mincount(isubject)),:);
    
    [n,bin] = histc(low(:,10),-pi:2*pi/phase_nrbins:pi);
    bin_rating_low  = nan(1,phase_nrbins);
    for ibin = 1 : phase_nrbins
        bin_rating_low(ibin)    = mean(low(bin==ibin,2));
    end
    
    [n,bin] = histc(high(:,10),-pi:2*pi/phase_nrbins:pi);
    bin_rating_high  = nan(1,phase_nrbins);
    for ibin = 1 : phase_nrbins
        bin_rating_high(ibin)    = mean(high(bin==ibin,2));
    end

    bar([-pi:2*pi/phase_nrbins:pi-pi/phase_nrbins]+pi/phase_nrbins/2,[bin_rating_low; bin_rating_high]',1.2);
    
    axis tight;
    i = i + 1;
    MI_high(isubject) = (log(phase_nrbins)-(-sum((bin_rating_high / sum(bin_rating_high))  .* log((bin_rating_high /sum(bin_rating_high))))))  /log(phase_nrbins);
    MI_low(isubject)  = (log(phase_nrbins)-(-sum((bin_rating_low  / sum(bin_rating_low ))  .* log((bin_rating_low  /sum(bin_rating_low ))))))   /log(phase_nrbins);

    title(sprintf('H: %0.4f L: %0.4f',MI_high(isubject),MI_low(isubject)));
end

figure;
subplot(1,2,1);
boxplot(mincount(slist));
ylabel('Count');

subplot(1,2,2);
boxplot([MI_high(slist);MI_low(slist)]');
set(gca,'XTickLabel',{'Low','High'});
ylabel('MI');

[h,p,ci,stats] = ttest(MI_high(slist),MI_low(slist));








header = {'trialtype', 'rating', 'correct', 'stimduration', 'RT', 'trialnr', 'samplenr_onset', 'samplenr_offset', 'EGG_phase_onset','EGG_phase_offset', 'IBI_first_HB', 'IBI_last_HB', 'IBI_avg_over','trialnr_over_blocks', 'EGG_amplitude_onset', 'EGG_amplitude_offset','EGG_amplitude_avg','subjectnr','rating_split'};

% slist = [slist nan nan nan];

% plot original with trend
for ilist = 0 : 4
    fig = figure;
    try
        i = 1;
        for isubject = slist((1:5)+ilist*5)
            for iblock = 1 : 4
                trials = find(blocknr{isubject} == iblock);
                subplot(5,4,iblock+(i-1)*4); hold on
                time = (behav{isubject}(trials,7) - behav{isubject}(trials(1),7)) / 2000;
                stairs(time,behav{isubject}(trials,2),'color','k','LineWidth',0.5,'Marker','o','MarkerFaceColor','w','MarkerSize',2);
                axis tight
                
                [p,~,mu] = polyfit(time, behav{isubject}(trials,2), 8);
                f = polyval(p,time,[],mu);
                plot(time,f,'color','r')
                
                ylim([0 9]);
                title(sprintf('Subject %d Block %d',isubject,iblock));
            end
            i = i + 1;
        end
    catch
    end
    set(fig,'PaperSize',[8.5*2 11*2 ]);
    set(fig,'PaperPosition',[0 0 11*2 8.5*2]);
    set(fig,'PaperOrientation','landscape');
    set(fig,'Position',[50 50 1200 800]);
    print(fig,'-dpdf',['d:\analysis\WANDER\images\trended_ratings_' num2str(ilist) '.pdf']);
end


% plot detrended
for ilist = 0 : 4
    fig = figure;
    try
        i = 1;
        for isubject = slist((1:5)+ilist*5)
            for iblock = 1 : 4
                trials = find(blocknr{isubject} == iblock);
                subplot(5,4,iblock+(i-1)*4); hold on
                time = (behav{isubject}(trials,7) - behav{isubject}(trials(1),7)) / 2000;
                
                [p,~,mu] = polyfit(time, behav{isubject}(trials,2), 4);
                f = polyval(p,time,[],mu);
                stairs(time,behav{isubject}(trials,2)-f,'color','k','LineWidth',0.5,'Marker','o','MarkerFaceColor','w','MarkerSize',2);
                axis tight
                
                ylim([-5 5]);
                title(sprintf('Subject %d Block %d',isubject,iblock));
            end
            i = i + 1;
        end
    catch
    end
    set(fig,'PaperSize',[8.5*2 11*2 ]);
    set(fig,'PaperPosition',[0 0 11*2 8.5*2]);
    set(fig,'PaperOrientation','landscape');
    set(fig,'Position',[50 50 1200 800]);
    print(fig,'-dpdf',['d:\analysis\WANDER\images\detrended_ratings_' num2str(ilist) '.pdf']);
end


%
%
% save('D:\analysis\WANDER\data\behaviour\data_behavior_WANDER.mat','behaviour')
% save('D:\analysis\WANDER\data\behaviour\header_behavior_WANDER.mat','header')
%
%
% dlmwrite('D:\analysis\WANDER\data\behaviour\data_behaviour.txt',behaviour)
%
% dlmwrite('D:\analysis\WANDER\data\behaviour\data_behaviour.txt','-append',behaviour)
% xlswrite('D:\analysis\WANDER\data\behaviour\data_behaviour.txt',[header; num2str(behaviour)])
%%

for isubject = slist
    EGG_amp_high_begin(isubject)    = nanmean(amp_begin{isubject}(rating_split{isubject} == 2));
    EGG_amp_low_begin(isubject)     = nanmean(amp_begin{isubject}(rating_split{isubject} == 1));
    EGG_amp_diff_begin(isubject)    = (EGG_amp_low_begin(isubject) - EGG_amp_high_begin(isubject)) / (EGG_amp_low_begin(isubject) + EGG_amp_high_begin(isubject) );
    
    EGG_amp_high_end(isubject)      = nanmean(amp_end{isubject}(rating_split{isubject} == 2));
    EGG_amp_low_end(isubject)       = nanmean(amp_end{isubject}(rating_split{isubject} == 1));
    EGG_amp_diff_end(isubject)    = (EGG_amp_low_end(isubject) - EGG_amp_high_end(isubject)) / (EGG_amp_low_end(isubject) + EGG_amp_high_end(isubject) );
    
    EGG_amp_high_avg(isubject)      = nanmean(amp_avg{isubject}(rating_split{isubject} == 2));
    EGG_amp_low_avg(isubject)       = nanmean(amp_avg{isubject}(rating_split{isubject} == 1));
    EGG_amp_diff_avg(isubject)    = (EGG_amp_low_avg(isubject) - EGG_amp_high_avg(isubject)) / (EGG_amp_low_avg(isubject) + EGG_amp_high_avg(isubject) );
    
    EGG_amp_high_norm_avg(isubject)      = EGG_amp_high_avg(isubject) / (EGG_amp_low_avg(isubject) + EGG_amp_high_avg(isubject) );
    EGG_amp_low_norm_avg(isubject)       = EGG_amp_low_avg(isubject) / (EGG_amp_low_avg(isubject) + EGG_amp_high_avg(isubject) );
    
    EGG_amp_high_norm_begin(isubject)      = EGG_amp_high_begin(isubject) / (EGG_amp_low_begin(isubject) + EGG_amp_high_begin(isubject) );
    EGG_amp_low_norm_begin(isubject)       = EGG_amp_low_begin(isubject) / (EGG_amp_low_begin(isubject) + EGG_amp_high_begin(isubject) );
    
    EGG_amp_high_norm_end(isubject)      = EGG_amp_high_end(isubject) / (EGG_amp_low_end(isubject) + EGG_amp_high_end(isubject) );
    EGG_amp_low_norm_end(isubject)       = EGG_amp_low_end(isubject) / (EGG_amp_low_end(isubject) + EGG_amp_high_end(isubject) );
    
end



% Bayesian tests
addpath('D:/analysis/WANDER/scripts/Mariana/');

% Prior: corresponds to an effect differing from 0 with a p-value of 0.05
nobs = 22; % number of observations
xref = +0.443; % reference effect size (should correspond to a significant effect)
disp('reference effect significance') % check that pref < 0.05
tref = xref*sqrt(nobs);
pref = 2*tcdf(-abs(tref),nobs-1);

[~, bf_begin, res_Bayes_begin]              = run_ttest_bayes(EGG_amp_high_begin(slist),        EGG_amp_low_begin(slist), xref);
[~, bf_end, res_Bayes_end]                  = run_ttest_bayes(EGG_amp_high_end(slist),          EGG_amp_low_end(slist), xref);
[~, bf_avg, res_Bayes_avg]                  = run_ttest_bayes(EGG_amp_high_avg(slist),          EGG_amp_low_avg(slist), xref);

[h,p,ci,stats_begin]                        = ttest(EGG_amp_high_begin(slist),                  EGG_amp_low_begin(slist));
[h,p,ci,stats_end]                          = ttest(EGG_amp_high_end(slist),                    EGG_amp_low_end(slist));
[h,p,ci,stats_avg]                          = ttest(EGG_amp_high_avg(slist),                    EGG_amp_low_avg(slist));

[~, bf_norm_begin, res_Bayes_norm_begin]    = run_ttest_bayes(EGG_amp_high_norm_begin(slist),   EGG_amp_low_norm_begin(slist), xref);
[~, bf_norm_end, res_Bayes_norm_end]        = run_ttest_bayes(EGG_amp_high_norm_end(slist),     EGG_amp_low_norm_end(slist), xref);
[~, bf_norm_avg, res_Bayes_norm_avg]        = run_ttest_bayes(EGG_amp_high_norm_avg(slist),     EGG_amp_low_norm_avg(slist), xref);

[h,p,ci,stats_norm_begin]                   = ttest(EGG_amp_high_norm_begin(slist),             EGG_amp_low_norm_begin(slist));
[h,p,ci,stats_norm_end]                     = ttest(EGG_amp_high_norm_end(slist),               EGG_amp_low_norm_end(slist));
[h,p,ci,stats_norm_avg]                     = ttest(EGG_amp_high_norm_avg(slist),               EGG_amp_low_norm_avg(slist));

fig = figure;
subplot(2,3,1);
bar([EGG_amp_high_begin(slist); EGG_amp_low_begin(slist)]');
title(sprintf('Begin - t: %0.3f, bf: %0.3f',stats_begin.tstat,bf_begin));
axis tight;

subplot(2,3,2);
bar([EGG_amp_high_end(slist); EGG_amp_low_end(slist)]');
title(sprintf('End - t: %0.3f, bf: %0.3f',stats_end.tstat,bf_end));
axis tight;

subplot(2,3,3);
bar([EGG_amp_high_avg(slist); EGG_amp_low_avg(slist)]');
title(sprintf('Avg - t: %0.3f, bf: %0.3f',stats_avg.tstat,bf_avg));
axis tight;

subplot(2,3,4);
bar([EGG_amp_diff_begin(slist)]');
title(sprintf('Begin - t: %0.3f, bf: %0.3f',stats_norm_begin.tstat,bf_norm_begin));
axis tight;

subplot(2,3,5);
bar([EGG_amp_diff_end(slist)]');
title(sprintf('End - t: %0.3f, bf: %0.3f',stats_norm_end.tstat,bf_norm_end));
axis tight;

subplot(2,3,6);
bar([EGG_amp_diff_avg(slist)]');
title(sprintf('Avg - t: %0.3f, bf: %0.3f',stats_norm_avg.tstat,bf_norm_avg));
axis tight;
%
%
% figure;
% boxplot([EGG_amp_high_norm_begin(slist); EGG_amp_low_norm_begin(slist)]')

set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0 11*4 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['d:\analysis\WANDER\images\EGG_power_diff.pdf']);
print(fig,'-dpng',['d:\analysis\WANDER\images\EGG_power_diff.png']);

%%

% IBI_trial_avg{isubject}(itrial)         = nanmean(IBI{isubject}(trial_indx,17));
% IBI_trial_onset{isubject}(itrial,:)     = IBI{isubject}(trial_indx(1),:);
% IBI_trial_offset{isubject}(itrial,:)    = IBI{isubject}(trial_indx(end),:);

% investigate circular distribution

for isubject = slist
    disp(num2str(isubject));
    for irand = 1 : 1000
        
        F = ceil(2 * tiedrank(rating{isubject}) / length(rating{isubject}));
        if irand > 1
            r = randperm(size(F,2));
            F = F(r);
        end
        
%         
%         circ_stat_begin{isubject}(irand)        = circ_stats(phase_begin{isubject});
%         circ_stat_begin_high{isubject}(irand)	= circ_stats(phase_begin{isubject}(F==1));
%         circ_stat_begin_low{isubject}(irand)	= circ_stats(phase_begin{isubject}(F==2));
%         
%         circ_r_begin(isubject,irand)            = circ_r(phase_begin{isubject}');
%         circ_r_begin_high(isubject,irand)       = circ_r(phase_begin{isubject}(F==1)');
%         circ_r_begin_low(isubject,irand)        = circ_r(phase_begin{isubject}(F==2)');
%         
%         circ_stat_end{isubject}(irand)          = circ_stats(phase_end{isubject});
%         circ_stat_end_high{isubject}(irand)     = circ_stats(phase_end{isubject}(F==1));
%         circ_stat_end_low{isubject}(irand)      = circ_stats(phase_end{isubject}(F==2));
%         
%         circ_r_end(isubject,irand)              = circ_r(phase_end{isubject}');
%         circ_r_end_high(isubject,irand)         = circ_r(phase_end{isubject}(F==1)');
%         circ_r_end_low(isubject,irand)          = circ_r(phase_end{isubject}(F==2)');
    end
end

fig = figure; hold;
plot(circ_r_end(slist,2:end));

plot(median(circ_r_end(slist,2:end),2));
plot(circ_r_begin(slist,1));


fig = figure;
%
% [h,p,ci,stats] = ttest(circ_r_begin(slist))
%
% subplot(1,4,1);
% boxplot([circ_r_begin(slist)']);
% title(sprintf('Begin - t: %0.3f, p: %0.3f',stats.tstat,p));
%
% [h,p,ci,stats] = ttest(circ_r_end(slist))
%
% subplot(1,4,2);
% boxplot([circ_r_end(slist)']);
% title(sprintf('End - t: %0.3f, p: %0.3f',stats.tstat,p));

[~, bf, res_Bayes]  = run_ttest_bayes(circ_r_begin_high(slist), circ_r_begin_low(slist), xref);
[h,p,ci,stats] = ttest(circ_r_begin_high(slist),circ_r_begin_low(slist))

subplot(1,2,1);
boxplot([circ_r_begin_high(slist)' circ_r_begin_low(slist)']);
title(sprintf('Begin - t: %0.3f, p: %0.3f, bf: %0.3f',stats.tstat,p,bf));

[~, bf, res_Bayes]  = run_ttest_bayes(circ_r_end_high(slist), circ_r_end_low(slist), xref);
[h,p,ci,stats] = ttest(circ_r_end_high(slist),circ_r_end_low(slist))

subplot(1,2,2);
boxplot([circ_r_end_high(slist)' circ_r_begin_low(slist)']);
title(sprintf('End - t: %0.3f, p: %0.3f, bf: %0.3f',stats.tstat,p,bf));

set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0 11*4 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['d:\analysis\WANDER\images\EGG_phase_corr.pdf']);
print(fig,'-dpng',['d:\analysis\WANDER\images\EGG_phase_corr.png']);


% Use VanRullen 2016 - Non-significant difference
for isubject = slist
    [p_circWW_begin(isubject), p_POS_begin(isubject), p_zPOS_begin(isubject)]  = PhaseOpposition(phase_begin{isubject}(rating_split{isubject}==1),phase_begin{isubject}(rating_split{isubject}==2),5000);
    [p_circWW_end(isubject),   p_POS_end(isubject),   p_zPOS_end(isubject)]    = PhaseOpposition(phase_end{isubject}(rating_split{isubject}==1),phase_end{isubject}(rating_split{isubject}==2),5000);
end
result_begin = combine_pvalues(p_POS_begin(slist),2,1);
result_end   = combine_pvalues(p_POS_end(slist),2,1);

% checking consistency
figure;
plot([p_POS_begin(slist); p_zPOS_begin(slist)]');

fig = figure;
subplot(1,2,1);
bar(p_POS_begin(slist));
title(sprintf('Begin: p = %0.3f',result_begin));
axis tight;
subplot(1,2,2);
bar(p_POS_end(slist));
title(sprintf('End: p = %0.3f',result_end));
axis tight;


fig = figure;
i = 1;
for isubject = slist
    F = ceil(2 * tiedrank(rating{isubject}) / length(rating{isubject}));

    % Calculate ITC+POS for the selected subject
    [ ITC_all(isubject), ITC_A(isubject), ITC_B(isubject), POS(isubject) ] = GUTSEE_calc_POS(phase_end{isubject}(F==1),phase_end{isubject}(F==2));
    
    % calculate some surrogate POS distribution
    [ surrPOS_distr(isubject,:) ] = GUTSEE_calc_surrogate_POS(phase_end{isubject}(F==1),phase_end{isubject}(F==2), 5000 );
    
    % single subject stat
    p_POS = length(find(surrPOS_distr(isubject,:)>POS(isubject)))/length(surrPOS_distr(isubject,:));

    subplot(5,5,i);
    hist(surrPOS_distr(isubject,:), 50);
    vline2(POS(isubject), 'r');
    title(sprintf('p: %0.3f',p_POS));
    i = i + 1;
end
POS_avg = mean(POS(slist));
% second level stat

for iperm = 1 : 100000
    
    ilist = 1;
    clear r
    for isubject = slist
        r(ilist) = surrPOS_distr(isubject,randi([1, size(surrPOS_distr,2)]));
        ilist = ilist + 1;
    end
    perm_dist(iperm) = mean(r);
end

figure; hold;
hist(perm_dist,100);
vline2(POS_avg, 'r');








% Compute the phase opposition statistics only for the subject selected
[p_circWW, p_POS, p_zPOS] = PhaseOpposition(phases_hits, phases_misses, 1000)

for isubject = slist
    [p_circWW(isubject), p_POS(isubject), p_zPOS(isubject)] = PhaseOpposition(phase_end{isubject}(rating_split{~isubject}),phase_end{isubject}(rating_split{isubject}),50000);
end
result = combine_pvalues(p_POS(slist),2,1)

% Use circ toolbox
%
%   [rho pval ts] = circ_corrcl(alpha, x)
%   Correlation coefficient between one circular and one linear random
%   variable.
%
%   Input:
%     alpha   sample of angles in radians
%     x       sample of linear random variable
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value

% for isubject = slist
%     [rho(isubject) pval(isubject)] = circ_corrcl(phase_begin{isubject}, rating{isubject});
% end
%
% figure;
% i = 1;
% for isubject = slist
%     subplot(5,5,i);
% %     polarhistogram('BinEdges',phase_bins,'BinCounts',rating_bincount_end{isubject},'FaceAlpha',0.5);
%     polarhistogram(phase_end{isubject}(rating_split{isubject}),'normalization','count','FaceAlpha',0.5,'FaceColor',[0 1 0]); hold;
%     polarhistogram(phase_end{isubject}(~rating_split{isubject}),'normalization','count','FaceAlpha',0.5,'FaceColor',[1 0 0]); axis tight;
%     %     polarhistogram(phase_end{isubject}(rating_split{isubject}),'normalization','count','DisplayStyle','stairs'); hold; axis tight
%     %     polarhistogram(phase_end{isubject}(~rating_split{isubject}),'normalization','count','DisplayStyle','stairs'); axis tight
%     %     circ_plot(phase_end{isubject}(rating_split{isubject})','pretty','ro',true,'linewidth',2,'color','r');
% %     circ_plot(phase_begin{isubject}(rating_split{isubject})','hist',[],30,true,true,'linewidth',2,'color','r');
%     i = i + 1;
% end
%
% figure;
% i = 1;
% for isubject = slist
%     subplot(5,5,i);
%     circ_plot(phase_end{isubject}(rating_split{~isubject})','hist',[],30,true,true,'linewidth',2,'color','r'); hold;
%     i = i + 1;
% end

window = 4;
figure;
i = 1;
for isubject = slist
    
    behav{isubject}(:,20) = nan(size(size(behav{isubject},1),1));
    for iblock = 1 : 4
        trials = find(blocknr{isubject} == iblock);
        trials = trials(window+1:end-window-1);
        for itrial = trials
            behav{isubject}(itrial,20) = nanmean(behav{isubject}(itrial-window:itrial+window,2));
        end
    end
    behav{isubject}(:,21) = behav{isubject}(:,2) - behav{isubject}(:,20);
    subplot(5,5,i); hold;
    plot(behav{isubject}(:,2),'k');
    %     plot(behav{isubject}(:,20),'r','linewidth',2);
    ylim([0 9]);
    title(num2str(isubject));
    
    i = i + 1;
end

figure;
i = 1;
for isubject = slist
    subplot(5,5,i); hold;
    plot(behav{isubject}(:,21),'k');
    %     plot(behav{isubject}(:,20),'r','linewidth',2);
    title(num2str(isubject));
    i = i + 1;
end

