function [EYE,EYE_avg,EYE_avg_low, EYE_avg_high,EYE_fix,EYE_avg_fix,EYE_avg_fix_low,EYE_avg_fix_high] = WANDER_EYE(isubject,force,rootpath)

if rootpath == 1
    fname_source_eye = ['z:\WANDER\data\EYE\s' num2str(isubject) '.mat'];
else
    fname_source_eye = ['/shared/projects/project_wander/WANDER/data/EYE/s' num2str(isubject) '.mat'];
end

% addpath('D:\fieldtrip\fieldtrip.git\trunk');
% addpath('D:/analysis/WANDER/scripts/');
% addpath D:\analysis\WANDER\scripts\Inpaint_nans

ft_defaults

if exist(fname_source_eye,'file') && force ~= 1
    fprintf('Returning Eyelink analysis \n');
    load(fname_source_eye);
else
    fprintf('Eyelink analysis not found, creating it now! \n');
    [dataset_task, dataset_rs] = WANDER_subjectinfo;
    
    %         artdef = WANDER_blink_detection(isubject,0,rootpath,0);
    
    for ipart = 1 : 4
        
        % load data as continuous for high-pass filter
        cfg = [];
        cfg.hpfilter                    = 'yes';
        cfg.hpfilter                    = 2;
        cfg.detrend                     = 'no';
        cfg.continuous                  = 'yes';
        cfg.demean                      = 'no';
        cfg.dftfilter                   = 'yes';
        cfg.lpfilter                    = 'yes';
        cfg.lpfreq                      = 150;
        cfg.hpfilter                    = 'no';
        cfg.dataset                     = dataset_task{isubject,ipart};
        cfg.channel                     = {'MISC009'};
        data_EYE{ipart}                 = ft_preprocessing(cfg);
        
        % define trials
        cfg                             = [];
        cfg.dataset                     = dataset_task{isubject,ipart};
        cfg.trialfun                    = 'WANDER_trialfun_cuelocked';
        cfg.trialdef.stim_pre           = 1;            % time before onset stimulation
        cfg.trialdef.stim_post          = 3;            % time after offset stimulation
        cfg.trialdef.onset              = [129 131];    % trigger for start stimulation
        cfg.trialdef.offset             = 3;            % trigger for end trial
        cfg.trialdef.stim               = [128 130];    % trigger for stimulations
        cfg.trialdef.rating             = 10:19;
        cfg.trialdef.performance        = 20:29;
        cfg                             = ft_definetrial(cfg);
        data_EYE{ipart}                 = ft_redefinetrial(cfg, data_EYE{ipart});
        data_EYE{ipart}.trl             = cfg.trl;
        
        % redefine timeaxis to probelocked
        for itrial = 1 : size(data_EYE{ipart}.trial,2)
            data_EYE{ipart}.time{itrial} = data_EYE{ipart}.time{itrial} - data_EYE{ipart}.time{itrial}(end) + 3;
        end
        
        cfg = [];
        cfg.latency = [-10 3];
        data_EYE{ipart} = ft_selectdata(cfg,data_EYE{ipart});
        
        
        %         % semi-automatic blink detection based on vertical EOG
        %         cfg = [];
        %         cfg.artfctdef.zvalue.channel    = {'BIO002','BIO012'};
        %         cfg.artfctdef.zvalue.cutoff     = 0.20;
        %         cfg.artfctdef.zvalue.bpfilter   = 'yes';
        %         cfg.artfctdef.zvalue.bpfilter   = 'yes';
        %         cfg.artfctdef.zvalue.bpfilttype = 'but';
        %         cfg.artfctdef.zvalue.bpfreq     = [4 15];
        %         cfg.artfctdef.zvalue.bpfiltord  = 4;
        %         cfg.artfctdef.zvalue.hilbert    = 'yes';
        %         cfg.artfctdef.zvalue.interactive= 'yes';
        %         cfg.artfctdef.zvalue.artpadding = 0.1;
        %
        %         [cfg, artifact_eye] = ft_artifact_zvalue(cfg, data_EYE{ipart});
        %
        %         % replace artifacts with nans
        %         cfg = [];
        %         cfg.reject                      = 'nan';
        %         cfg.artfctdef.zvalue.artifact   = artifact_eye;
        %         data_EYE{ipart}                 = ft_rejectartifact(cfg,data_EYE{ipart});
        %         data_EYE_fix{ipart}             = data_EYE{ipart};
        %         eyechan                         = find(strcmp('MISC009', data_EYE{ipart}.label));
        
        % add trialnum
        for itrial = 1 : size(data_EYE{ipart}.trialinfo,1)
            data_EYE{ipart}.trialinfo(itrial,7) = itrial + (ipart-1)*50;
        end
        
        % semi-automatic blink detection on EYELINK
        cfg = [];
        cfg.trl = data_EYE{ipart}.trl;
        cfg.artfctdef.jump.channel    = 'MISC009';
        cfg.artfctdef.jump.cutoff = 2;
        
        %         cfg.artfctdef.zvalue.cutoff     = 1;
        cfg.artfctdef.jump.bpfilter   = 'no';
        %         cfg.artfctdef.zvalue.bpfilttype = 'but';
        %         cfg.artfctdef.zvalue.bpfreq     = [30 200];
        %         cfg.artfctdef.zvalue.bpfiltord  = 4;
        %         cfg.artfctdef.zvalue.hilbert    = 'yes';
        cfg.artfctdef.jump.interactive= 'no';
        cfg.artfctdef.jump.artpadding = 0.15;
        
        [cfg, artifact_jump] = ft_artifact_jump(cfg, data_EYE{ipart});
        
        % replace artifacts with nans
        cfg = [];
        cfg.reject                      = 'nan';
        cfg.artfctdef.zvalue.artifact   = [artifact_jump];
        data_EYE_fix{ipart}    = ft_rejectartifact(cfg,data_EYE{ipart});
        
        % remove data that is beyond threshold
        eyechan                         = find(strcmp('MISC009', data_EYE_fix{ipart}.label));
        for itrial = 1 : size(data_EYE_fix{ipart}.trial,2)
            remove = find(abs(data_EYE_fix{ipart}.trial{itrial}(eyechan,:)) > 8);
            data_EYE_fix{ipart}.trial{itrial}(eyechan,remove) = nan;
        end
        
        % remove trials that have to many nans
        remove = [];
        for itrial = 1 : size(data_EYE_fix{ipart}.trial,2)
            if mean(isnan(data_EYE_fix{ipart}.trial{itrial}(eyechan,:))) > 0.5
                remove = [remove itrial];
            end
            mean(isnan(data_EYE_fix{ipart}.trial{itrial}(eyechan,:)))
        end
        
        remove_indx = true(size(data_EYE_fix{ipart}.trialinfo,1),1);
        remove_indx(remove) = false;
        
        cfg                     = [];
        cfg.trials              = find(remove_indx);
        data_EYE_fix{ipart}  = ft_selectdata(cfg,data_EYE_fix{ipart});
        
        % interpolate NaNs
        for itrial = 1 : size(data_EYE_fix{ipart}.trial,2)
            if any(~isnan(data_EYE_fix{ipart}.trial{itrial}(eyechan,:))) > 0
                data_EYE_fix{ipart}.trial{itrial}(eyechan,:) = fixgaps(data_EYE_fix{ipart}.trial{itrial}(eyechan,:));
                %                 data_EYE_fix{isubject,ipart}.trial{itrial}(eyechan,:) = inpaint_nans(data_EYE_fix{isubject,ipart}.trial{itrial}(eyechan,:),4);
            end
        end
        
        % remove trials that still have (trialing) nans
        remove = [];
        for itrial = 1 : size(data_EYE_fix{ipart}.trial,2)
            if any(isnan(data_EYE_fix{ipart}.trial{itrial}(eyechan,:)))
                remove = [remove itrial];
            end
            mean(isnan(data_EYE_fix{ipart}.trial{itrial}(eyechan,:)))
        end
        
        remove_indx = true(size(data_EYE_fix{ipart}.trialinfo,1),1);
        remove_indx(remove) = false;
        
        cfg                     = [];
        cfg.trials              = find(remove_indx);
        data_EYE_fix{ipart}  = ft_selectdata(cfg,data_EYE_fix{ipart});
        
        % %
        %         cfg = [];
        %         cfg.viewmode = 'vertical';
        %         cfg.detrend = 'yes';
        %         cfg.channel = 'MISC009';
        %         ft_databrowser(cfg, data_EYE_fix{isubject,ipart}  );
    end
    
    EYE     = ft_appenddata([],data_EYE{:});
    EYE_fix = ft_appenddata([],data_EYE_fix{:});
    
    clear data*
    [dataset_task, dataset_rs] = WANDER_subjectinfo;
    
    % reject all but correct rejections
    cfg                 = [];
    cfg.trials          = find(EYE.trialinfo(:,3) == 4);
    EYE                 = ft_selectdata(cfg,EYE);
    cfg.trials          = find(EYE_fix.trialinfo(:,3) == 4);
    EYE_fix             = ft_selectdata(cfg,EYE_fix);
    
    % put data in matrix
    dat     = NaN(size(EYE.trial,2),10001+3000);
    dat_fix = NaN(size(EYE_fix.trial,2),10001+3000);
    
    %     dat_fix = NaN(size(EYE_fix{isubject}.trial,2),size(EYE_fix{isubject}.trial{1},1),size(EYE_fix{isubject}.trial{1},2));
    
    for itrial = 1 : size(EYE.trial,2)
        fprintf('Adding trial %d to data matrix \n',itrial);
        temp            = EYE.trial{itrial};
        dat(itrial,end-size(temp,2)+1:end)   = temp;
    end
    for itrial = 1 : size(EYE_fix.trial,2)
        fprintf('Adding trial %d to data matrix \n',itrial);
        temp            = EYE_fix.trial{itrial};
        dat_fix(itrial,end-size(temp,2)+1:end)   = temp;
    end
    
    %       % put data in matrix
    %     dat     = NaN(size(EYE{isubject}.trial,2),size(EYE{isubject}.trial{1},1),size(EYE{isubject}.trial{1},2));
    %     dat_fix = NaN(size(EYE_fix{isubject}.trial,2),size(EYE_fix{isubject}.trial{1},1),size(EYE_fix{isubject}.trial{1},2));
    %
    %     for itrial = 1 : size(EYE{isubject}.trial,2)
    %         fprintf('Adding trial %d to data matrix \n',itrial);
    %         dat(itrial,:,:)     = EYE{isubject}.trial{itrial};
    %     end
    %     for itrial = 1 : size(EYE_fix{isubject}.trial,2)
    %         fprintf('Adding trial %d to data matrix \n',itrial);
    %         dat_fix(itrial,:,:)     = EYE_fix{isubject}.trial{itrial};
    %     end
    time = linspace(-10,3,13*1000+1);
    %
    %         figure;
    %         plot(time,nanmean(dat_fix));
    %         figure;
    %         imagesc(dat_fix)
    %
    % manually average over all trials, excluding nans
    cfg = [];
    cfg.vartrllength            = 2;
    EYE_avg           = ft_timelockanalysis(cfg,EYE);
    EYE_avg.avg       = nanmean(dat,1);
    EYE_avg           = rmfield(EYE_avg,'var');
    
    cfg = [];
    cfg.vartrllength            = 2;
    EYE_avg_fix       = ft_timelockanalysis(cfg,EYE_fix);
    EYE_avg_fix.avg   = nanmean(dat_fix,1);
    EYE_avg_fix       = rmfield(EYE_avg_fix,'var');
    
    % make index according to median split based on ratings
    F = ceil(2 * tiedrank(EYE.trialinfo(:,2)) / length(EYE.trialinfo(:,2)));
    rating_split        = ones(size(F));
    rating_split(F==2)  = 2;
    
    EYE_avg_high          = EYE_avg;
    EYE_avg_high.avg      = nanmean(dat(rating_split==2,:,:),1);
    EYE_avg_low           = EYE_avg;
    EYE_avg_low.avg       = nanmean(dat(rating_split==1,:,:),1);
    
    F = ceil(2 * tiedrank(EYE_fix.trialinfo(:,2)) / length(EYE_fix.trialinfo(:,2)));
    rating_split        = ones(size(F));
    rating_split(F==2)  = 2;
    
    EYE_avg_fix_high      = EYE_avg_fix;
    EYE_avg_fix_high.avg  = nanmean(dat_fix(rating_split==2,:,:),1);
    EYE_avg_fix_low       = EYE_avg_fix;
    EYE_avg_fix_low.avg   = nanmean(dat_fix(rating_split==1,:,:),1);
    
    %         figure; hold; title(num2str(isubject));
    for irating = 1:8
        if any(EYE_fix.trialinfo(:,2)==irating)
            rating_present(irating) = true;
            EYE_avg_fix_rating{irating} = EYE_avg_fix;
            EYE_avg_fix_rating{irating}.avg = nanmean(dat_fix(EYE_fix.trialinfo(:,2) == irating,:),1);
            %                 plot(time,EYE_avg_fix_rating{isubject}{irating}.avg);
        else
            rating_present(irating) = false;
        end
    end
    
    EYE = rmfield(EYE,'cfg');
    EYE_avg = rmfield(EYE_avg,'cfg');
    EYE_avg_low = rmfield(EYE_avg_low,'cfg');
    EYE_avg_high = rmfield(EYE_avg_high,'cfg');
    EYE_fix = rmfield(EYE_fix,'cfg');
    EYE_avg_fix = rmfield(EYE_avg_fix,'cfg');
    EYE_avg_fix_low = rmfield(EYE_avg_fix_low,'cfg');
    EYE_avg_fix_high = rmfield(EYE_avg_fix_high,'cfg');
    
    % save data
    fprintf('Saving data in %s. This can take a while.\n',fname_source_eye);
    save(fname_source_eye,'EYE','EYE_avg','EYE_avg_low','EYE_avg_high','EYE_fix','EYE_avg_fix','EYE_avg_fix_low','EYE_avg_fix_high','EYE_avg_fix_rating','rating_present','-v7.3');
end
