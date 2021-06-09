
function [data_MEG] = WANDER_epoch_MEG(isubject,force,rootpath,restingstate)

if restingstate == 1
    if rootpath == 1
        fname_epoched_MEG = ['w:\WANDER\data\trial\s' num2str(isubject) '_epoched_MEG_restingstate.mat'];
    else
        fname_epoched_MEG = ['/shared/projects/project_wander/WANDER/data/trial/s' num2str(isubject) '_epoched_MEG_restingstate.mat'];
    end
else
    if rootpath == 1
        fname_epoched_MEG = ['w:\WANDER\data\trial\s' num2str(isubject) '_epoched_MEG_cuelocked.mat'];
    else
        fname_epoched_MEG = ['/shared/projects/project_wander/WANDER/data/trial/s' num2str(isubject) '_epoched_MEG_cuelocked.mat'];
    end
end

if exist(fname_epoched_MEG,'file') && force ~= 1
    fprintf('Returning epoched MEG data\n');
    load(fname_epoched_MEG);
else
    fprintf('Epoched MEG data not found, creating it now! \n');
    addpath('D:/analysis/WANDER/scripts/');
    [dataset_task, dataset_rs] = WANDER_subjectinfo(rootpath);
    
    if restingstate == 1
        % define trial
        cfg                         = [];
        cfg.dataset                 = dataset_rs{isubject};
        cfg.trialfun                = 'WANDER_trialfun_restingstate';
        cfg                         = ft_definetrial(cfg);
        
        % load data as continuous segment
        cfg.continuous              = 'yes';
        cfg.hpfilter                = 'no';
        cfg.detrend                 = 'no';
        cfg.continuous              = 'yes';
        cfg.demean                  = 'yes';
        cfg.dftfilter               = 'yes';
        cfg.dftfreq                 = [50 100 150 200];
        cfg.lpfilter                = 'yes';
        cfg.lpfreq                  = 150;
        cfg.channel                 = 'all';
        cfg.hpfilter                = 'yes';
        cfg.hpfreq                  = 0.5;
        cfg.hpfilttype              = 'fir';
        data_MEG                    = ft_preprocessing(cfg);
    else
        for ipart = 1:4
            
            % load data as continuous segment
            cfg                         = [];
            cfg.continuous              = 'yes';
            cfg.hpfilter                = 'no';
            cfg.detrend                 = 'no';
            cfg.continuous              = 'yes';
            cfg.demean                  = 'yes';
            cfg.dftfilter               = 'yes';
            cfg.dftfreq                 = [50 100 150 200];
            cfg.lpfilter                = 'yes';
            cfg.lpfreq                  = 150;
            cfg.dataset                 = dataset_task{isubject,ipart};
            cfg.channel                 = 'all';
            
            cfg.hpfilter                = 'yes';
            cfg.hpfreq                  = 0.5;
            cfg.hpfilttype              = 'fir';
            data_MEG{ipart}             = ft_preprocessing(cfg);
            
            % define trials
            cfg                         = [];
            cfg.dataset                 = dataset_task{isubject,ipart};
            cfg.trialfun                = 'WANDER_trialfun_cuelocked';
            cfg.trialdef.stim_pre       = 1.5;          % time before onset stimulation
            cfg.trialdef.stim_post      = 1;            % time after offset stimulation
            cfg.trialdef.onset          = [129 131];    % trigger for start stimulation
            cfg.trialdef.offset         = 3;            % trigger for end trial
            cfg.trialdef.stim           = [128 130];    % trigger for stimulations
            cfg.trialdef.rating         = 10:19;
            cfg.trialdef.performance    = 20:29;
            cfg                         = ft_definetrial(cfg);
            
            % epoch data
            data_MEG{ipart}             = ft_redefinetrial(cfg,data_MEG{ipart});
            
            % select relevant data with workaround since preprocessing doesnt take labels and lists at same time
            cfg                         = [];
            cfg.channel                 = 'MEG';
            temp_MEG                    = ft_selectdata(cfg,data_MEG{ipart});
            cfg.channel                 = {'BIO001','BIO002','BIO003','BIO012'};
            temp_BIO                    = ft_selectdata(cfg,data_MEG{ipart});
            cfg.channel                 = {'MISC*'};
            temp_MISC                   = ft_selectdata(cfg,data_MEG{ipart});
            data_MEG{ipart}             = ft_appenddata([],temp_MEG,temp_BIO,temp_MISC);
        end
    end
    save(fname_epoched_MEG,'data_MEG','-v7.3');
end

