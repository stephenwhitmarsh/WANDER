
function [artdef] = WANDER_blink_detection(isubject,force,rootpath,restingstate)

if restingstate == 1
    if rootpath == 1
        fname_artdef = ['d:\analysis\WANDER\data\artefacts\s' num2str(isubject) '_artdef_EOG_cuelocked_rs.mat'];
    else
        fname_artdef = ['/home/swhitmarsh/WANDER/data/artefacts/s' num2str(isubject) '_artdef_EOG_cuelocked_rs.mat'];
    end
else
    if rootpath == 1
        fname_artdef = ['d:\analysis\WANDER\data\artefacts\s' num2str(isubject) '_artdef_EOG_cuelocked.mat'];
    else
        fname_artdef = ['/home/swhitmarsh/WANDER/data/artefacts/s' num2str(isubject) '_artdef_EOG_cuelocked.mat'];
    end
end

if exist(fname_artdef,'file') && force ~= 1
    fprintf('Returning blink detection artefact definitions \n');
    load(fname_artdef);
else
    fprintf('Blink detection artefact definition not found, creating it now! \n');
    addpath('D:/analysis/WANDER/scripts/');
    WANDER_subjectinfo;
    
    fprintf('Loading epoched data \n');
    data_epoched = WANDER_epoch_MEG(isubject,0,rootpath,restingstate);
    clear temp
    
    if restingstate == 1
        cfg = [];
        
        % remove events from visualization
        data_epoched.cfg.event       = [];
        
        % semi-automatic eyeblink detection
        cfg.artfctdef.zvalue.channel        = {'BIO002','BIO012'};
        cfg.artfctdef.zvalue.cutoff         = 2;
        cfg.artfctdef.zvalue.trlpadding     = 0;
        cfg.artfctdef.zvalue.artpadding     = 0.1;
        cfg.artfctdef.zvalue.fltpadding     = 0;
        cfg.artfctdef.zvalue.bpfilter       = 'yes';
        cfg.artfctdef.zvalue.bpfilttype     = 'but';
        cfg.artfctdef.zvalue.bpfreq         = [1 15];
        cfg.artfctdef.zvalue.bpfiltord      = 4;
        cfg.artfctdef.zvalue.hilbert        = 'yes';
        cfg.artfctdef.zvalue.interactive    = 'no';
        cfg                                 = ft_artifact_zvalue(cfg,data_epoched);
        
        artdef = cfg.artfctdef;
    else
        for ipart = 1:4
            
            cfg = [];
            
            % remove events from visualization
            data_epoched{ipart}.cfg.event       = [];
            
            % semi-automatic eyeblink detection
            cfg.artfctdef.zvalue.channel        = {'BIO002','BIO012'};
            cfg.artfctdef.zvalue.cutoff         = 2;
            cfg.artfctdef.zvalue.trlpadding     = 0;
            cfg.artfctdef.zvalue.artpadding     = 0.1;
            cfg.artfctdef.zvalue.fltpadding     = 0;
            cfg.artfctdef.zvalue.bpfilter       = 'yes';
            cfg.artfctdef.zvalue.bpfilttype     = 'but';
            cfg.artfctdef.zvalue.bpfreq         = [1 15];
            cfg.artfctdef.zvalue.bpfiltord      = 4;
            cfg.artfctdef.zvalue.hilbert        = 'yes';
            cfg.artfctdef.zvalue.interactive    = 'no';
            cfg                                 = ft_artifact_zvalue(cfg,data_epoched{ipart});
            
            artdef{ipart} = cfg.artfctdef;
        end
    end
    save(fname_artdef,'artdef');
    
end
