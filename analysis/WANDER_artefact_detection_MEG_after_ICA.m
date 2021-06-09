function [artdef] = WANDER_artefact_detection_MEG_after_ICA(isubject,force,rootpath,restingstate)

if restingstate == 1
    if rootpath == 1
        fname_artdef = ['d:\analysis\WANDER\data\artefacts\s' num2str(isubject) '_artdef_MEG_cuelocked_rs.mat'];
    else
        fname_artdef = ['/shared/projects/project_wander/WANDER/data/artefacts/s' num2str(isubject) '_artdef_MEG_cuelocked_rs.mat'];
    end
else
    if rootpath == 1
        fname_artdef = ['d:\analysis\WANDER\data\artefacts\s' num2str(isubject) '_artdef_MEG_cuelocked.mat'];
    else
        fname_artdef = ['/shared/projects/project_wander/WANDER/data/artefacts/s' num2str(isubject) '_artdef_MEG_cuelocked.mat'];
    end
end

if ~exist(fname_artdef,'file')
    fprintf('Artefact definition not found! \n');
else
    load(fname_artdef);
    
    if restingstate == 1
        
        data_MEG = WANDER_ICA(isubject,0,rootpath,restingstate);
        
        % visual MEG artefact rejection
        cfg                                 = [];
        cfg.viewmode                        = 'vertical';
        cfg.ploteventlabels                 = 'colorvalue';
        cfg.preproc.demean                  = 'yes';
        cfg.ylim                            = 'maxmin';
        cfg.channel                         = {'MEGMAG'};
        cfg.blocksize = 60;
        %         cfg.artfctdef.MEG.artifact          = artdef.MEG.artifact;
        try
            cfg.artfctdef.muscle.artifact    = artdef.zvalue.artifact;
        catch
        end
        try
            cfg.artfctdef.MEG.artifact       = artdef.MEG.artifact;
            cfg.artfctdef.muscle.artifact    = artdef.muscle.artifact;
        catch
        end
        %         cfg.event                           = 1; % hack to clear events for display
        cfg.selectfeature                   = 'MEG';
        cfg                                 = ft_databrowser(cfg,data_MEG);
        
        artdef = cfg.artfctdef;
    else
        data_MEG = WANDER_ICA(isubject,0,rootpath,restingstate);
        
        for ipart = 1:4
            
            % visual MEG artefact rejection
            cfg                                 = [];
            cfg.viewmode                        = 'vertical';
            cfg.ploteventlabels                 = 'colorvalue';
            cfg.preproc.demean                  = 'yes';
            cfg.ylim                            = 'maxmin';
            cfg.channel                         = {'MEGMAG'};
            cfg.blocksize = 60;
            
            try
                cfg.artfctdef.MEG.artifact      = artdef{ipart}.MEG.artifact;
                cfg.artfctdef.muscle.artifact   = artdef{ipart}.muscle.artifact;
            catch
            end
            
            %         cfg.event                           = 1; % hack to clear events for display
            cfg.selectfeature                   = 'MEG';
            cfg                                 = ft_databrowser(cfg,data_MEG{ipart});
            artdef{ipart} = cfg.artfctdef;
        end
    end
    save(fname_artdef,'artdef');
    
end

