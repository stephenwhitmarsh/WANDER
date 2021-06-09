function [behav] = WANDER_behaviour(isubject,force,rootpath)

if rootpath == 1
    fname_behav = ['d:\analysis\WANDER\data\behaviour\s' num2str(isubject) '_behav.mat'];
else
    fname_behav = ['/shared/projects/project_wander/WANDER/data/behaviour/s' num2str(isubject) '_behav.mat'];
end

if exist(fname_behav,'file') && force ~= 1
    fprintf('Returning behavioural analysis \n');
    load(fname_behav);
else
    fprintf('Behavioural analysis not found, creating it now! \n');
    [dataset_task, dataset_rs] = WANDER_subjectinfo(rootpath);
    
    for ipart = 1:4
        % define trials
        cfg = [];
        cfg.dataset                 = dataset_task{isubject,ipart};
        cfg.trialfun                = 'WANDER_trialfun_probelocked';
        cfg.trialdef.stim_pre       = 0;            % time before onset stimulation
        cfg.trialdef.stim_post      = 0;            % time after offset stimulation
        cfg.trialdef.onset          = [129 131];    % trigger for start stimulation
        cfg.trialdef.offset         = 3;            % trigger for end trial
        cfg.trialdef.stim           = [128 130];    % trigger for stimulations
        cfg.trialdef.rating         = 10:19;
        cfg.trialdef.performance    = 20:29;
        trl{ipart}                  = ft_definetrial(cfg);
    end
    
    trialinfo                   = [ trl{1}.trl; trl{2}.trl; trl{3}.trl; trl{4}.trl; ];
    behav.stimduration          = trialinfo(:,7);
    behav.hits                  = size(find(trialinfo(:,6) == 1),1);
    behav.correct_rejections    = size(find(trialinfo(:,6) == 4),1);
    behav.misses                = size(find(trialinfo(:,6) == 2),1);
    behav.FAs                   = size(find(trialinfo(:,6) == 3),1);
    behav.rating                = 9-trialinfo(:,5);
    behav.trialinfo             = trialinfo(:,6);
    
    for irating = 1 : 8
        behav.stimduration_mean(irating)    = nanmean(behav.stimduration(behav.rating == irating));
        behav.stimduration_std(irating)     = nanstd(behav.stimduration(behav.rating == irating));
    end
    save(fname_behav,'behav');
end

if rootpath == 1
    
    fig = figure;
    
    subplot(2,2,1);
    hist(behav.rating',1:1:8);
    axis tight;
    xlabel('attention rating');
    ylabel('nr of trials');
    
    subplot(2,2,2);
    ax = gca;
    x = [1 2 3 4];
    y = [behav.hits;behav.correct_rejections;behav.misses;behav.FAs];
    bar(y);
    ax.XTickLabel = {'Hits','CR','Misses','FAs'};
    ylabel('nr of trials');
    for i1=1:4
        text(x(i1),y(i1),num2str(y(i1),'%d'),...
            'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
    xlim([0.5 4.5]);
    ylim([0 200]);
    
    subplot(2,2,3);
    hist(behav.stimduration,8:1:30);
    xlabel(sprintf('stimduration (M: %.2f, SD: %.2f)', mean(behav.stimduration),std(behav.stimduration)));
    ylabel('nr of trials');
    axis tight;
    
    subplot(2,2,4);
    bar(behav.stimduration_mean);
    axis tight;
    xlabel('attention rating');
    ylabel('stimulation duration');
    
%     print(fig,'-dpdf',['d:\analysis\WANDER\images\behaviour\s' num2str(isubject) '.pdf']);
    
end

