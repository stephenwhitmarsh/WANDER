% combine data
clear all

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

force = 0;
rootpath = 1;

if rootpath == 1
    datapath = 'z:';
    addpath('Z:/WANDER/scripts/');
else
    datapath = '/shared/projects/project_wander';
end

% get ROI (stat_FFT)
% from WANDER_source_individual_onlyalpha_GRAD_GA.m
load(fullfile(datapath,'/WANDER/data/stat/stat_FFT'),'stat_FFT');
load(fullfile(datapath,'/WANDER/data/stat/int_stat_FFT'),'int_stat_FFT');

for isubject = slist
    
    % need to add trialinfo of trials to source
    source = WANDER_source_individual_onlyalpha_GRAD_keeptrials(isubject,force,rootpath);
    [EYE,EYE_avg,EYE_avg_low, EYE_avg_high,EYE_fix,EYE_avg_fix,EYE_avg_fix_low,EYE_avg_fix_high] = WANDER_EYE(isubject,force,rootpath);
    load(fullfile(datapath,'WANDER','data','ERF',['s' num2str(isubject) '_ratings_singletrial.mat'])); %dip
        
    % select those trialnr as specificied in trialinfo that are retained in all measurements
    [c_all{isubject}] = intersect(intersect(source.trialinfo(:,7), dip.trialinfo(:,7)), EYE.trialinfo(:,7));
    
    % get the index
    [~, ~, c_eye]       = intersect(c_all{isubject}, EYE.trialinfo(:,7));
    [~, ~, c_dip]       = intersect(c_all{isubject}, dip.trialinfo(:,7));
    [~, ~, c_source]    = intersect(c_all{isubject}, source.trialinfo(:,7));

    cfg = [];
    cfg.trials = c_eye;
    EYE_sel{isubject} = ft_selectdata(cfg,EYE);
        
    cfg = [];
    cfg.trials = c_dip;
    dip_sel{isubject} = ft_selectdata(cfg,dip);
    
    cfg = [];
    cfg.trials = c_source;
    source_sel{isubject} = ft_selectdata(cfg,source);

    % take the rating from the source data (could be any after selecting trials)
    ratings{isubject} = source_sel{isubject}.trialinfo(:,2);
    trialinfo{isubject} = source_sel{isubject}.trialinfo;
    
    roi = find(stat_FFT.negclusterslabelmat == 1);
    
    source_roi{isubject} = nanmean(source_sel{isubject}.pow(roi,:),1);
    
    % EYE average over last 10 seconds
    for itrial = 1:size(EYE_sel{isubject}.trial,2)
        eye_roi{isubject}(itrial) = nanmean(EYE_sel{isubject}.trial{itrial}(1:10000));
    end
    
    % dip average power over last 10 seconds
    ifreq = find(dip_sel{isubject}.freq > 16,1,'first');
    for itrial = 1:size(dip_sel{isubject}.trial,2)
        dip_roi{isubject}(itrial) = nanmean(dip_sel{isubject}.pow{itrial}(1,ifreq-1:ifreq+1));
    end
    
    % block regressor
    blocknr{isubject} = ones(size(dip_sel{isubject}.trial,2),1);
    for itrial = 2:size(dip_sel{isubject}.trial,2)
        if source_sel{isubject}.trialinfo(itrial,6) < source_sel{isubject}.trialinfo(itrial-1,6)
            blocknr{isubject}(itrial:end) = blocknr{isubject}(itrial:end) + 1;
        end
    end
end



% nr of trials per subject after artefact rejection
for isubject = slist
    
    % need to add trialinfo of trials to source
    source = WANDER_source_individual_onlyalpha_GRAD_keeptrials(isubject,force,rootpath);
    [EYE,EYE_avg,EYE_avg_low, EYE_avg_high,EYE_fix,EYE_avg_fix,EYE_avg_fix_low,EYE_avg_fix_high] = WANDER_EYE(isubject,force,rootpath);
    load(fullfile(datapath,'WANDER','data','ERF',['s' num2str(isubject) '_ratings_singletrial.mat'])); %dip
        

    
    % select those trialnr as specificied in trialinfo that are retained in all measurements
    [c_all{isubject}] = intersect(intersect(source.trialinfo(:,7), dip.trialinfo(:,7)), EYE.trialinfo(:,7));
    
    nrtrials.alpha(isubject)    = size(source.trialinfo, 1); 
    nrtrials.dip(isubject)      = size(dip.trialinfo, 1);
    nrtrials.eye(isubject)      = size(EYE.trialinfo, 1);
    nrtrials.intcpt(isubject)   = size(c_all{isubject}, 1);
    
end

save('Z:/WANDER/revision_EJN/nrtrials', 'nrtrials');
mean(nrtrials.alpha(slist))
std(nrtrials.alpha(slist))
mean(nrtrials.dip(slist))
std(nrtrials.dip(slist))
mean(nrtrials.eye(slist))
std(nrtrials.eye(slist))
mean(nrtrials.intcpt(slist))
std(nrtrials.intcpt(slist))


% 
% 
% 
% 
% q = [];
% pow = [];
% ii = 1;
% for isubject = slist
% % 
% %     nrbins = 2;
% %     [ratings_binned, ratings_old, ratings_stat] = bin_data(dip_sel{isubject}.trialinfo(:,2),nrbins);
%         % median split
%         
%     F                                   = ceil(2 * tiedrank(dip_sel{isubject}.trialinfo(:,2)) / length(dip_sel{isubject}.trialinfo(:,2)));
%     rating_split                        = ones(size(F));
%     rating_split(F==2)                  = 2;
%     
%     cfg = [];
%     cfg.avgoverrpt = 'yes';
%     cfg.channel = 1;
%     cfg.trials = find(rating_split == 1);
%     pow_low{isubject} = ft_selectdata(cfg,dip_sel{isubject});
% 
%     cfg.trials = find(rating_split == 2);
%     pow_high{isubject} = ft_selectdata(cfg,dip_sel{isubject});    
% end
% 
% i = 1;
% clear dat_high dat_low
% for isubject = slist
%     dat_low(i,:) = pow_low{isubject}.pow{1};
%     dat_high(i,:) = pow_high{isubject}.pow{1};   
%     i = i + 1;
% end
% pow_low_avg = nanmean(dat_low);
% pow_high_avg = nanmean(dat_high);
% 
% figure; hold;
% plot(pow_low{1}.freq,pow_low_avg);
% plot(pow_low{1}.freq,pow_high_avg);
% 
% figure; 
% plot(pow_low{1}.freq,(pow_low_avg - pow_high_avg) ./
% end
% pow_low_avg = nanmean(dat_low);
% pow_high_avg = nanmean(dat_high);
% 
% figure; hold;
% plot(pow_low{1}.freq,pow_low_avg);
% plot(pow_low{1}.freq,pow_high_avg);
% 
% figure; 
% plot(pow_low{1}.freq,(pow_low_avg - pow_high_avg) ./ (pow_high_avg + pow_low_avg));


clear b_* VIF cormat predictors
for isubject = slist
%     predictors{isubject} = [eye_roi{isubject}', source_roi{isubject}', dip_roi{isubject}', blocknr{isubject}, trialinfo{isubject}(:,6), trialinfo{isubject}(:,8)];  
%     predictors{isubject} = [eye_roi{isubject}', source_roi{isubject}', dip_roi{isubject}', blocknr{isubject}, trialinfo{isubject}(:,6), trialinfo{isubject}(:,8) ];  
%         predictors{isubject} = [eye_roi{isubject}', source_roi{isubject}', dip_roi{isubject}', blocknr{isubject}, trialinfo{isubject}(:,6) ]; 

    predictors{isubject} = [eye_roi{isubject}', source_roi{isubject}', dip_roi{isubject}', blocknr{isubject}, trialinfo{isubject}(:,8), trialinfo{isubject}(:,4)];  

    % cut off first trial
    predictors{isubject} = predictors{isubject}(2:end,:);
 
    % add previous rating as predictor
    predictors{isubject}(:,7) = dip_sel{isubject}.trialinfo(1:end-1,2);

%     predictors{isubject} = [eye_roi{isubject}', source_roi{isubject}', dip_roi{isubject}'];
%   predictors{isubject} = [ dip_roi{isubject}'];
    predictors{isubject} = zscore(predictors{isubject});

    ratings{isubject} = zscore(dip_sel{isubject}.trialinfo(2:end,2));
    
%    
%     F                                   = ceil(2 * tiedrank(ratings{isubject}) / length(ratings{isubject}));
%     rating_split                        = ones(size(F));
%     rating_split(F==2)                  = 2;
%     
%     ratings{isubject} = zscore(rating_split);

%     figure; imagesc(predictors{isubject});
    
    b_all(isubject,:) = glmfit(predictors{isubject},ratings{isubject},'normal');
    
%     b_eye(isubject,:) = glmfit(zscore(eye_roi{isubject}'),ratings{isubject},'normal');
%     b_source(isubject,:) = glmfit(zscore(source_roi{isubject}'),ratings{isubject},'normal');
%     b_dip(isubject,:) = glmfit(zscore(dip_roi{isubject}'),ratings{isubject},'normal');

    matX = predictors{isubject};
    VIF(isubject, :) = diag(inv(corrcoef(matX)))';
    
    cormat(isubject,:,:) = corr(predictors{isubject});
% end
% pow_low_avg = nanmean(dat_low);
% pow_high_avg = nanmean(dat_high);
% 
% figure; hold;
% plot(pow_low{1}.freq,pow_low_avg);
% plot(pow_low{1}.freq,pow_high_avg);
% 
% figure; 
% plot(pow_low{1}.freq,(pow_low_avg - pow_high_avg) ./ (pow_high_avg + pow_low_avg));


clear b_* VIF cormat predictors
for isubject = slist
%     predictors{isubject} = [eye_roi{isubject}', source_roi{isubject}', dip_roi{isubject}', blocknr{isubject}, trialinfo{isubject}(:,6), trialinfo{isubject}(:,8)];  
%     predictors{isubject} = [eye_roi{isubject}', source_roi{isubject}', dip_roi{isubject}', blocknr{isubject}, trialinfo{isubject}(:,6), trialinfo{isubject}(:,8) ];  
%         predictors{isubject} = [eye_roi{isubject}', source_roi{isubject}', dip_roi{isubject}', blocknr{isubject}, trialinfo{isubject}(:,6) ]; 

    predictors{isubject} = [eye_roi{isubject}', source_roi{isubject}', dip_roi{isubject}', blocknr{isubject}, trialinfo{isubject}(:,8), trialinfo{isubject}(:,4)];  

    % cut off first trial
    predictors{isubject} = predictors{isubject}(2:end,:);
 
    % add previous rating as predictor
    predictors{isubject}(:,7) = dip_sel{isubject}.trialinfo(1:end-1,2);

%     predictors{isubject} = [eye_roi{isubject}', source_roi{isubject}', dip_roi{isubject}'];
%   predictors{isubject} = [ dip_roi{isubject}'];
    predictors{isubject} = zscore(predictors{isubject});

    ratings{isubject} = zscore(dip_sel{isubject}.trialinfo(2:end,2));
    
%    
%     F                                   = ceil(2 * tiedrank(ratings{isubject}) / length(ratings{isubject}));
%     rating_split                        = ones(size(F));
%     rating_split(F==2)                  = 2;
%     
%     ratings{isubject} = zscore(rating_split);

%     figure; imagesc(predictors{isubject});
    
    b_all(isubject,:) = glmfit(predictors{isubject},ratings{isubject},'normal');
    
%     b_eye(isubject,:) = glmfit(zscore(eye_roi{isubject}'),ratings{isubject},'normal');
%     b_source(isubject,:) = glmfit(zscore(source_roi{isubject}'),ratings{isubject},'normal');
%     b_dip(isubject,:) = glmfit(zscore(dip_roi{isubject}'),ratings{isubject},'normal');

    matX = predictors{isubject};
    VIF(isubject, :) = diag(inv(corrcoef(matX)))';
    
    cormat(isubject,:,:) = corr(predictors{isubject});
end

% create xml file for jerome
% SubjectID,Pupil,Alpha,SteadyState,Blocknr,Trialnr,Timestamp,Rating,Length, 

clear trialdata
load(fullfile(datapath,'WANDER','data','SSEF','rms'));

trialdata = [];
for isubject = slist
    for itrial = 1 : size(trialinfo{isubject},1)
        trialdata = [trialdata; isubject, eye_roi{isubject}(itrial), source_roi{isubject}(itrial), dip_roi{isubject}(itrial), blocknr{isubject}(itrial), trialinfo{isubject}(itrial,6), trialinfo{isubject}(itrial,8), trialinfo{isubject}(itrial,2), trialinfo{isubject}(itrial,4) r{isubject}(itrial) ];
    end
end

dat = table;
dat.SubjectID   = trialdata(:,1);
dat.Pupil       = trialdata(:,2);
dat.Alpha       = trialdata(:,3);
dat.SteadyState = trialdata(:,4);
dat.Blocknr     = trialdata(:,5);
dat.Trialnr     = trialdata(:,6);
dat.Timestamp   = trialdata(:,7);

% xlswrite(fullfile('C:','Users','swhitmarsh','Dropbox','WANDER_article','WANDER','R','trialdata_triallength_2.xls'),trialdata)
xlswrite(fullfile(datapath,'WANDER','data','trialdata_triallength_2.xls'),trialdata);






%% ANALYSE FURTHER WITH R






load(fullfile(datapath,'/WANDER/data/stat/stat_FFT'),'stat_FFT');

figure; plot(xcorr(ratings{isubject}))

% 
% cormat_avg = squeeze(nanmean(cormat(slist,:,:),1));
% [h,p,ci,stats] = ttest(cormat(slist,:,:));
% figure; imagesc(squeeze(stats.tstat));


figure; boxplot(b_all(slist,:))
figure; boxplot(VIF(slist,:))

[h,p,ci,stats] = ttest(b_all(slist,:))

[h,p,ci,stats] = ttest(b_eye(slist,:))
[h,p,ci,stats] = ttest(b_source(slist,:))
[h,p,ci,stats] = ttest(b_dip(slist,:))

 cfg =[];
 cfg.avgoverrpt = 'yes';
 t = ft_selectdata(cfg,dip);
 
 b1 = glmfit(predictors,ratings,'normal')
 [b2,bint,r,rint,stats] = regress(ratings,predictors,0.05)
 
    

 i = 1;
for isubject = slist
xc(isubject) = xcorr(
i = i + 1;
end
 
 
 