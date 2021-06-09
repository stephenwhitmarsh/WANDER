function [HEF_avg,HEF_avg_low,HEF_avg_high,rating_ibi_avg] = WANDER_Heart(isubject,force,timing,rootpath)

addpath('d:/analysis/WANDER/scripts/heart_functions-master/');

if ~(strcmp(timing,'cue') || strcmp(timing,'probe'))
    fprintf('Use "cue" or "probe" as third argument.\n');
    return
end

if rootpath == 1
    fname_heart = ['i:\analysis\WANDER\data\heart\s' num2str(isubject) '.mat'];
else
    fname_heart = ['/home/swhitmarsh/WANDER/data/heart/s' num2str(isubject) '.mat'];
end

if exist(fname_heart,'file') && force ~= 1
    fprintf('Returning ERF\n');
    load(fname_heart);
else
    fprintf('Heart-locked ERF not found, creating it now! \n');
    data_epoch_MEG = WANDER_ICA(isubject,0,rootpath,0);
    artdef = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
    
    [dataset_task, dataset_rs] = WANDER_subjectinfo;
    
    % start with cuelocked
    for ipart = 1 : 4
        
        cfg = [];
        cfg.fsample = 1000;
        cfg.dataset = dataset_task{isubject,ipart};
        cfg.channel = 'BIO003';
        [HeartBeats, R_time] = heart_peak_detect(cfg);
        
        clear R RR_sort
        for i = 1 : size(HeartBeats,2)
            R(i) = HeartBeats(i).R_sample;
        end
        figure; hist( diff(R),100)
        [Y RR_sort] = sort(diff(R),'descend');
        R_high = [HeartBeats(RR_sort(1:5)).R_sample];
        R_low  = [HeartBeats(RR_sort(end-5:end)).R_sample];
        R_diff = [NaN diff(R)];
        
        cfg = [];
        cfg.dataset = dataset_task{isubject,ipart};
        cfg.channel = 'BIO003';
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 1;
        data = ft_preprocessing(cfg);
%         
%         cfg = [];
%         cfg.viewmode = 'vertical';
%         cfg.blocksize = 15;
%         cfg.artfctdef.R.artifact      = [R'-50 R'+50];
%         cfg.artfctdef.R_high.artifact = [R_high'-50 R_high'+50];
%         cfg.artfctdef.R_low.artifact  = [R_low'-50 R_low'+50];
%         ft_databrowser(cfg,data);
        
%         str = input('Let me know when you are ready with an ENTER','s');
        
%         apply trial segmentation to MEG data
        cfg = [];
        cfg.trl = [R'-300 R'+600 ones(size(R,2),1)*-300];
        cfg.trl = cfg.trl(2:end-1,:);
        cfg.dataset = dataset_task{isubject,ipart}; %temp = rmfield(data_epoch_MEG{ipart},'trialinfo');
        cfg.hpfilter = 'yes';
        cfg.hpfreq = 1;
        Hlocked{ipart} = ft_preprocessing(cfg);
        
%         % artefact rejection
%         cfg = [];
%         cfg.reject                      = 'nan';
%         cfg.artfctdef.minaccepttim      = 0;
%         cfg.artfctdef                   = artdef{ipart};
%         cfg.channel                     = 'MEG';
%         temp_MEG                        = ft_rejectartifact(cfg,temp_MEG);        
        
        % create index of trialnr for each sample in data
        trialindx = [];
        for itrial = unique(data_epoch_MEG{ipart}.trialinfo(:,6))'
            trialindx(data_epoch_MEG{ipart}.sampleinfo(itrial,1) : data_epoch_MEG{ipart}.sampleinfo(itrial,2)) = itrial;
        end
        
        % add trial number and trialinfo to each heart beat segment
        
        for ipeak = 1 : size(Hlocked{ipart}.sampleinfo,1)
            Hlocked{ipart}.trialinfo(ipeak,:) = nan(1,18);
            try
                Hlocked{ipart}.trialinfo(ipeak,1:6) = data_epoch_MEG{ipart}.trialinfo(trialindx(HeartBeats(ipeak).R_sample),:);
                Hlocked{ipart}.trialinfo(ipeak,18) = trialindx(HeartBeats(ipeak).R_sample);   
            catch
            end
            Hlocked{ipart}.trialinfo(ipeak,7)  = HeartBeats(ipeak).P_sample;
            Hlocked{ipart}.trialinfo(ipeak,8)  = HeartBeats(ipeak).P_time;
            Hlocked{ipart}.trialinfo(ipeak,9)  = HeartBeats(ipeak).Q_sample;
            Hlocked{ipart}.trialinfo(ipeak,10) = HeartBeats(ipeak).Q_time;
            Hlocked{ipart}.trialinfo(ipeak,11) = HeartBeats(ipeak).R_sample;
            Hlocked{ipart}.trialinfo(ipeak,12) = HeartBeats(ipeak).R_time;
            Hlocked{ipart}.trialinfo(ipeak,13) = HeartBeats(ipeak).S_sample;
            Hlocked{ipart}.trialinfo(ipeak,14) = HeartBeats(ipeak).S_time;
            Hlocked{ipart}.trialinfo(ipeak,15) = HeartBeats(ipeak).T_sample;
            Hlocked{ipart}.trialinfo(ipeak,16) = HeartBeats(ipeak).T_time;
            Hlocked{ipart}.trialinfo(ipeak,17) = R_diff(ipeak);
        end
        
        % remove beats outside of trial time
        cfg = [];
        cfg.trials = find(~isnan(Hlocked{ipart}.trialinfo(:,1)));
        Hlocked{ipart} = ft_selectdata(cfg,Hlocked{ipart});
        
        % remove beats of which the next doesn't show an R peak, in one subject
        if isubject == 7
            cfg = [];
            cfg.trials = find(Hlocked{ipart}.trialinfo(:,17) < 1300 & Hlocked{ipart}.trialinfo(:,17) > 700) ;        
            Hlocked{ipart} = ft_selectdata(cfg,Hlocked{ipart});
        end
%         
%         % only correct rejections
%         cfg         = [];
%         cfg.trials  = Hlocked{ipart}.trialinfo(:,3) == 4;
%         Hlocked{ipart} = ft_selectdata(cfg,Hlocked{ipart});
        
    end
    
    % concatinate trials
    Hlocked                             = ft_appenddata([],Hlocked{:});
    
    % reject all but correct rejections
    cfg                                 = [];
    cfg.trials                          = find(Hlocked.trialinfo(:,3) == 4);
    
    % plot IBI vs rating
    for irating = 1 : 8
        indx = Hlocked.trialinfo(:,2) == irating;
        rating_ibi_avg(irating) = nanmean(Hlocked.trialinfo(indx,17));
        rating_ibi_std(irating) = nanstd(Hlocked.trialinfo(indx,17));
    end
%     figure; errorbar(rating_ibi_avg,rating_ibi_std);
    
    % put data in matrix for nanmean averaging
    dat = NaN(size(Hlocked.trial,2),size(Hlocked.label,1),size(Hlocked.trial{1},2));
    for itrial = 1 : size(Hlocked.trial,2)
        fprintf('Adding trial %d or %d to data matrix \n',itrial,size(Hlocked.trial,2));
        dat(itrial,:,:) = Hlocked.trial{itrial};
    end
    
    % manually nanmean over all trials
    HEF_avg = [];
    HEF_avg.avg         = squeeze(nanmean(dat,1));
    HEF_avg.var         = squeeze(nanstd(dat,1));
    HEF_avg.time        = Hlocked.time{1};
    HEF_avg.dimord      = 'chan_time';
    HEF_avg.label       = Hlocked.label;
    HEF_avg.trialinfo   = Hlocked.trialinfo; % changed
    
    % make index according to median split based on ratings
    F = ceil(2 * tiedrank(Hlocked.trialinfo(:,2)) / length(Hlocked.trialinfo(:,2)));
    
    HEF_avg_low             = HEF_avg;
    HEF_avg_low.avg         = squeeze(nanmean(dat(F == 1,:,:),1));
    HEF_avg_low.var         = squeeze(nanstd(dat(F == 1,:,:),1));
    HEF_avg_low.trialinfo   = Hlocked.trialinfo(F == 1,:);
    
    HEF_avg_high            = HEF_avg;
    HEF_avg_high.avg        = squeeze(nanmean(dat(F == 2,:,:),1));
    HEF_avg_high.var        = squeeze(nanstd(dat(F == 2,:,:),1));
    HEF_avg_high.trialinfo  = Hlocked.trialinfo(F == 2,:);
    
    % plot IBI for high/low ratings
    figure;
    subplot(1,2,1); boxplot(HEF_avg_high.trialinfo(:,17)); ylim([500 1700]);
    subplot(1,2,2); boxplot(HEF_avg_low.trialinfo(:,17)); ylim([500 1700]);
    
    % Plot topo HEF
    cfg         = [];
    cfg.layout  = 'neuromag306mag';
    cfg.channel = 'MEG*1';
    
    %     cfg.layout  = 'neuromag306cmb';
    %     cfg.channel = 'MEG*3';
    
    cfg.channel = 'MISC007';
    
%     fig = figure; ft_singleplotER(cfg,HEF_avg);
    
%     fig = figure; ft_singleplotER(cfg,HEF_avg_low,HEF_avg_high);
    disp('saving data');
    save(fname_heart,'HEF_avg','HEF_avg_low','HEF_avg_high','rating_ibi_avg','rating_ibi_std');
end





