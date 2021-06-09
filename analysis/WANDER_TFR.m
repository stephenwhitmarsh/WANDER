function [TFR] = WANDER_TFR(isubject,force,timing,rootpath,restingstate)

if ~(strcmp(timing,'cue') || strcmp(timing,'probe'))
    fprintf('Use "cue" or "probe" as third argument.\n');
    return
end

if restingstate == 1
    if rootpath == 1
        fname_TFR = ['w:\WANDER\data\TFR\s' num2str(isubject) '_TFR_rs.mat'];
    else
        fname_TFR = ['/shared/projects/project_wander/WANDER/data/TFR/s' num2str(isubject) '_TFR_rs.mat'];
    end
else
    if rootpath == 1
        fname_TFR = ['w:\WANDER\data\TFR\s' num2str(isubject) '_TFR_' timing '.mat'];
    else
        fname_TFR = ['/shared/projects/project_wander/WANDER/data/TFR/s' num2str(isubject) '_TFR_' timing '.mat'];
    end
end

if exist(fname_TFR,'file') && force ~= 1
    fprintf('Returning TFR\n');
    load(fname_TFR);
else
    
    if restingstate == 1
        data_MEG = WANDER_ICA(isubject,0,rootpath,restingstate);

        fprintf('TFR restingstate not found, creating it now! \n');
        
        cfg              = [];
        cfg.pad          = 'nextpow2';
        cfg.channel      = 'all';
        cfg.method       = 'mtmconvol';
        cfg.keeptrials   = 'no';
        cfg.channel      = {'MEG'};
        cfg.foi          = 1:1:30;
        cfg.taper        = 'hanning';
        cfg.toi          = 0:0.05:12.5*60;
        cfg.t_ftimwin    = ones(size(cfg.foi));
        TFR              = ft_freqanalysis(cfg, data_MEG);
        
    else
        fprintf('TFR not found, creating it now! \n');
        
        % load data
        data_epoch_MEG = WANDER_ICA(isubject,0,rootpath,restingstate);
        
        if strcmp(timing,'cue')
            for ipart  = 1 : 4
                for itrial = 1 : size(data_epoch_MEG{ipart}.trial,2)
                    fprintf('Cutting off beginning of trial %d block %d \n',itrial,ipart);
                    data_epoch_MEG{ipart}.trial{itrial} = data_epoch_MEG{ipart}.trial{itrial}(:,1501:end);
                    data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial}(1501:end);                    
                    fprintf('Cutting off end of trial %d block %d \n',itrial,ipart);
                    data_epoch_MEG{ipart}.trial{itrial} = data_epoch_MEG{ipart}.trial{itrial}(:,1:end-1000);
                    data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial}(1:end-1000);
                end
            end
        end
        
        if strcmp(timing,'probe')
            for ipart  = 1 : 4
                for itrial = 1 : size(data_epoch_MEG{ipart}.trial,2)
                    fprintf('Cutting off beginning of trial %d block %d \n',itrial,ipart);
                    data_epoch_MEG{ipart}.trial{itrial} = data_epoch_MEG{ipart}.trial{itrial}(:,1501:end);
                    data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial}(1501:end);
                    fprintf('Cutting off end of trial %d block %d \n',itrial,ipart);
                    data_epoch_MEG{ipart}.trial{itrial} = data_epoch_MEG{ipart}.trial{itrial}(:,1:end-1000);
                    data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial}(1:end-1000);                    
                    fprintf('Shifting timecourse of trial %d block %d \n',itrial,ipart);
                    data_epoch_MEG{ipart}.time{itrial}  = data_epoch_MEG{ipart}.time{itrial} - data_epoch_MEG{ipart}.time{itrial}(end);
                end
            end
        end
        
        for ipart = 1 : 4
            cfg              = [];
            cfg.pad          = 'nextpow2';
            cfg.channel      = 'all';
            cfg.method       = 'mtmconvol';
            cfg.keeptrials   = 'yes';
            cfg.channel      = {'MEG'};
            cfg.foi          = 1:1:30;
            cfg.taper        = 'hanning';
            if strcmp(timing,'cue')
                cfg.toi      = 0:0.05:30;
            end
            if strcmp(timing,'probe')
                cfg.toi      = -30:0.05:0;
            end
            cfg.t_ftimwin    = ones(size(cfg.foi));
            TFR{ipart}       = ft_freqanalysis(cfg, data_epoch_MEG{ipart});
        end
    end
    
    save(fname_TFR,'TFR','-v7.3');
end

