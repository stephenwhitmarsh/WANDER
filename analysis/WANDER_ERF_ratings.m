function ERF_probe_ratings = WANDER_ERF_ratings(isubject,force,rootpath)

if rootpath == 1
    fname_ERF_ratings = ['W:\WANDER\data\ERF\s' num2str(isubject) '_ratings.mat'];
else
    fname_ERF_ratings = ['/shared/projects/project_wander/WANDER/data/ERF/s' num2str(isubject) '_ratings.mat'];
end

if exist(fname_ERF_ratings,'file') && force ~= 1
    fprintf('Returning ERF\n');
    load(fname_ERF_ratings);
else
    fprintf('ERF not found, creating it now! \n');
    
    data_epoch_MEG          = WANDER_ICA(isubject,0,rootpath,0);
    artdef                  = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
    
    % now lock to end of trial
    for ipart = 1 : 4
        cfg = [];
        cfg.reject                      = 'nan'; % make bugreport, does not work with cfg.artfctdef.reject
        cfg.artfctdef.minaccepttim      = 0;
        cfg.artfctdef                   = artdef{ipart};
        data_epoch_MEG{ipart}           = ft_rejectartifact(cfg,data_epoch_MEG{ipart});
        for itrial = 1 : size(data_epoch_MEG{ipart}.trial,2)
            fprintf('Cutting off beginning of trial %d block %d \n',itrial,ipart);
            data_epoch_MEG{ipart}.trial{itrial} = data_epoch_MEG{ipart}.trial{itrial}(:,1501:end);
            data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial}(1501:end);
            fprintf('Shifting timecourse of trial %d block %d \n',itrial,ipart);
            data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial} - data_epoch_MEG{ipart}.time{itrial}(end) + 1.000;
        end
    end
    
    % concatinate trials
    data_epoch_MEG                      = ft_appenddata([],data_epoch_MEG{:});
    
    % reject all but correct rejections, only grad
    cfg                                 = [];
    cfg.channel                         = 'MEG*1';
    cfg.trials                          = find(data_epoch_MEG.trialinfo(:,3) == 4);
    data_epoch_MEG                      = ft_selectdata(cfg,data_epoch_MEG);
    
    % maximum length trial
    maxlength = 0;
    for itrial = 1 : size(data_epoch_MEG.trial,2)
        if size(data_epoch_MEG.trial{itrial},2) > maxlength
            maxlength = size(data_epoch_MEG.trial{itrial},2);
        end
    end
    
    % put data in matrix
    dat = NaN(size(data_epoch_MEG.trial,2),size(data_epoch_MEG.label,1),maxlength);
    for itrial = 1 : size(data_epoch_MEG.trial,2)
        fprintf('Adding trial %d to data matrix \n',itrial);
        triallength = size(data_epoch_MEG.trial{itrial},2);
        dat(itrial,:,maxlength-triallength+1:maxlength) = data_epoch_MEG.trial{itrial};
    end
    
    nrbins = 4;
    [ratings_binned, ratings_old, ratings_stat] = bin_data(data_epoch_MEG.trialinfo(:,2),nrbins);
    
    cfg = [];
    cfg.vartrllength    = 2;
    ERF_probe           = ft_timelockanalysis(cfg,data_epoch_MEG);
    ERF_probe           = rmfield(ERF_probe,'var');
    ERF_probe           = rmfield(ERF_probe,'cfg');

    for irating = 1:nrbins
        ERF_probe_ratings{irating} = rmfield(ERF_probe,'dof');
        ERF_probe_ratings{irating}.avg = squeeze(nanmean(dat(ratings_binned == irating,:,:),1));
        ERF_probe_ratings{irating}.nrtrials = length(find(ratings_binned == irating));
        ERF_probe_ratings{irating}.perctrials = length(find(ratings_binned == irating)) / length(ratings_binned);
        
    end
    
    fprintf('Saving data, this will take a while... \n');
    save(fname_ERF_ratings,'ERF_probe_ratings','-v7.3');
    
end