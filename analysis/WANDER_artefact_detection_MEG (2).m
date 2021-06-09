function [artdef] = WANDER_artefact_detection_MEG(isubject,force,rootpath,restingstate)

if restingstate == 1
    if rootpath == 1
        fname_artdef = ['d:\analysis\WANDER\data\artefacts\s' num2str(isubject) '_artdef_MEG_cuelocked_rs.mat'];
    else
        fname_artdef = ['/home/swhitmarsh/WANDER/data/artefacts/s' num2str(isubject) '_artdef_MEG_cuelocked_rs.mat'];
    end
else
    if rootpath == 1
        fname_artdef = ['d:\analysis\WANDER\data\artefacts\s' num2str(isubject) '_artdef_MEG_cuelocked.mat'];
    else
        fname_artdef = ['/home/swhitmarsh/WANDER/data/artefacts/s' num2str(isubject) '_artdef_MEG_cuelocked.mat'];
    end
end

if exist(fname_artdef,'file') && force ~= 1
    fprintf('Returning artefact definitions \n');
    load(fname_artdef);
else
    fprintf('Artefact definition not found, creating it now! \n');
    addpath('D:/analysis/WANDER/scripts/');
    WANDER_subjectinfo;
    
    fprintf('Loading epoched data \n');
    data_MEG = WANDER_epoch_MEG(isubject,0,rootpath,restingstate);
    
    if restingstate == 1
       
        % semi-automatic muscle detection
        cfg                                 = [];
        cfg.artfctdef.zvalue.channel        = 'MEG';
        cfg.artfctdef.zvalue.cutoff         = 20;
        cfg.artfctdef.zvalue.trlpadding     = 0;
        cfg.artfctdef.zvalue.fltpadding     = 0;
        cfg.artfctdef.zvalue.artpadding     = 0.1;
        cfg.artfctdef.zvalue.bpfilter       = 'yes';
        cfg.artfctdef.zvalue.bpfreq         = [110 140];
        cfg.artfctdef.zvalue.bpfiltord      = 9;
        cfg.artfctdef.zvalue.bpfilttype     = 'but';
        cfg.artfctdef.zvalue.hilbert        = 'yes';
        cfg.artfctdef.zvalue.boxcar         = 0.2;
        cfg.artfctdef.zvalue.interactive    = 'no';
        [cfg, artifact_muscle]              = ft_artifact_zvalue(cfg,data_MEG);
        
        % visual MEG artefact rejection
%         cfg                                 = [];
%         cfg.viewmode                        = 'vertical';
%         cfg.ploteventlabels                 = 'colorvalue';
%         cfg.preproc.demean                  = 'yes';
%         cfg.ylim                            = 'maxmin';
%         cfg.channel                         = {'MEGMAG'};
%         cfg.blocksize = 60;
%         cfg.artfctdef.MEG.artifact          = [];
%         cfg.artfctdef.muscle.artifact       = artifact_muscle;
%         %         cfg.event                           = 1; % hack to clear events for display
%         cfg.selectfeature                   = 'MEG';
%         cfg                                 = ft_databrowser(cfg,data_MEG);
%         
        artdef = cfg.artfctdef;
    else
        for ipart = 1:4
            
            % remove events from visualization
            data_epoched{ipart}.cfg.event       = [];
            
            % semi-automatic muscle detection
            cfg                                 = [];
            cfg.artfctdef.zvalue.channel        = 'MEG';
            cfg.artfctdef.zvalue.cutoff         = 20;
            cfg.artfctdef.zvalue.trlpadding     = 0;
            cfg.artfctdef.zvalue.fltpadding     = 0;
            cfg.artfctdef.zvalue.artpadding     = 0.1;
            cfg.artfctdef.zvalue.bpfilter       = 'yes';
            cfg.artfctdef.zvalue.bpfreq         = [110 140];
            cfg.artfctdef.zvalue.bpfiltord      = 9;
            cfg.artfctdef.zvalue.bpfilttype     = 'but';
            cfg.artfctdef.zvalue.hilbert        = 'yes';
            cfg.artfctdef.zvalue.boxcar         = 0.2;
            cfg.artfctdef.zvalue.interactive    = 'no';
            [cfg, artifact_muscle]              = ft_artifact_zvalue(cfg,data_MEG{ipart});
            
            % visual MEG artefact rejection
            cfg                                 = [];
            cfg.viewmode                        = 'vertical';
            cfg.ploteventlabels                 = 'colorvalue';
            cfg.preproc.demean                  = 'yes';
            cfg.ylim                            = 'maxmin';
            cfg.channel                         = {'MEGMAG'};
            cfg.artfctdef.MEG.artifact          = [];
            cfg.artfctdef.muscle.artifact       = artifact_muscle;
            cfg.selectfeature                   = 'MEG';
            cfg                                 = ft_databrowser(cfg,data_MEG{ipart});
%             
            artdef{ipart} = cfg.artfctdef;
        end
    end
    save(fname_artdef,'artdef');
    
end

