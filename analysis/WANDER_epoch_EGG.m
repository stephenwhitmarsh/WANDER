
function [data_EGG] = WANDER_epoch_EGG(isubject,rootpath,force)

timing = 'cue';

if rootpath == 1
    fname_epoched_EGG = ['w:\WANDER\data\trial\s' num2str(isubject) '_epoched_EGG_' timing '.mat'];
else
    fname_epoched_EGG = ['/shared/projects/project_wander/WANDER/data/trial/s' num2str(isubject) '_epoched_EGG_' timing '.mat'];
end


if exist(fname_epoched_EGG,'file') && force~=1
    fprintf('Returning epoched EGG data\n');
    load(fname_epoched_EGG);
else
    fprintf('Epoched EGG data not found, creating it now! \n');
    WANDER_subjectinfo;
    
    fprintf('Loading filtered EGG data \n');
    data_EGG = WANDER_filter_EGG(isubject,0,rootpath,0);

    for ipart = 1:4
        
        % define trials
        cfg = [];
        cfg.dataset                 = dataset{isubject,ipart};        
        cfg.trialfun                = 'WANDER_trialfun_cuelocked';
        cfg.trialdef.stim_pre       = 1.5;          % time before onset stimulation
        cfg.trialdef.stim_post      = 1;            % time after offset stimulation
        cfg.trialdef.onset          = [129 131];    % trigger for start stimulation
        cfg.trialdef.offset         = 3;            % trigger for end trial
        cfg.trialdef.stim           = [128 130];    % trigger for stimulations
        cfg.trialdef.rating         = 10:19;
        cfg.trialdef.performance    = 20:29;
        cfg                         = ft_definetrial(cfg);
        
        % load and preprocess data
        data_EGG{ipart}             = ft_redefinetrial(cfg,data_EGG{ipart});
        
    end
    save(fname_epoched_EGG,'data_EGG','-v7.3');
end

