function data_ICA = WANDER_ICA(isubject,force,rootpath,restingstate)

if restingstate == 1
    if rootpath == 1
        fname_ICA  = ['d:\analysis\WANDER\data\ICA\s' num2str(isubject) '_epoched_MEG_ICA_rs.mat'];
    else
        fname_ICA  = ['/home/swhitmarsh/WANDER/data/ICA\s' num2str(isubject) '_epoched_MEG_ICA_rs.mat'];
    end
else
    if rootpath == 1
        fname_ICA  = ['d:\analysis\WANDER\data\ICA\s' num2str(isubject) '_epoched_MEG_ICA.mat'];
    else
        fname_ICA  = ['/home/swhitmarsh/WANDER/data/ICA\s' num2str(isubject) '_epoched_MEG_ICA.mat'];
    end
end

if exist(fname_ICA,'file') && force ~= 1
    fprintf('Returning ICA cleaned data\n');
    load(fname_ICA);
else
    fprintf('ICA cleaned data not found, creating it now! \n');
    WANDER_subjectinfo;
    
    % load epoched data
    data_epoched = WANDER_epoch_MEG(isubject,0,1,restingstate);
    
    % load artefact definition
    artdef_MEG = WANDER_artefact_detection_MEG(isubject,0,1,restingstate);
    
    if restingstate == 1
        
        % partial rejection
        cfg = [];
        cfg.artfctdef           = artdef_MEG;
        cfg.artfctdef.reject    = 'partial';
        data_artefact_rejected  = ft_rejectartifact(cfg,data_epoched);
        
        % determine rank
        cfg = [];
        cfg.trials              = 1;
        cfg.channel             = 'MEG';
        temp                    = ft_selectdata(cfg,data_epoched);
        data_rank               = rank(temp.trial{1}*temp.trial{1}');
        clear temp
        
        % decompose
        cfg                     = [];
        cfg.method              = 'runica';
        cfg.runica.pca          = data_rank;
        cfg.channel             = 'MEG';
        cfg.randomseed          = WANDER_plantseeds(isubject,0);
        comp                    = ft_componentanalysis(cfg, data_artefact_rejected);
        
        EOGvi = find(ismember(data_artefact_rejected.label,'BIO002'));
        if isempty(EOGvi)
            EOGvi = find(ismember(data_artefact_rejected.label,'BIO012'));
        end
        ECGi  = find(ismember(data_artefact_rejected.label,'BIO003'));
        
        % correlate components with EOG & ECG
        for itrial = 1 : size(data_artefact_rejected.trial,2)
            for icomp = 1 : size(comp.label,1)
                fprintf('Correlating trial %d of %d\n',itrial,size(data_artefact_rejected.trial,2));
                %                 EOGh(itrial,icomp) = corr(data_artefact_rejected.trial{itrial}(EOGhi,:)',comp.trial{itrial}(icomp,:)');
                EOGv(itrial,icomp) = corr(data_artefact_rejected.trial{itrial}(EOGvi,:)',comp.trial{itrial}(icomp,:)');
                ECG(itrial,icomp)  = corr(data_artefact_rejected.trial{itrial}(ECGi,:)',comp.trial{itrial}(icomp,:)');
            end
        end
        EOGv_avg  = mean(EOGv,1);
        rejected_EOG = false(1,data_rank);
        
        % determine to-be-rejected components based on thresold, iteratively
        for iteration = 1 : 5
            rejected_EOG = abs(EOGv_avg) > std(EOGv_avg(~rejected_EOG))*3;
        end
        
        % plot correlations
        fig = figure;
        ah = axes('parent', fig);
        hold(ah, 'on');
        for ibar = 1 : numel(EOGv_avg)
            bar(ibar, EOGv_avg(ibar), 'parent', ah, 'facecolor', [1 0 0] .* rejected_EOG(ibar));
        end
        title('average correlation component with vertical EOG');
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_corr_comp_EOGv_rs.pdf']);
        
        % plot removed EOGv components
        cfg = [];
        cfg.component = find(rejected_EOG);
        cfg.layout = 'neuromag306mag.lay';
        cfg.marker = 'off';
        cfg.comment = 'no';
        fig = figure; ft_topoplotIC(cfg,comp);
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_topo_comp_EOGv_rs.pdf']);
        
        % determine to-be-rejected components based on threshold, iteratively
        ECG_avg  = mean(ECG,1);
        rejected_ECG = false(1,data_rank);
        for iteration = 1 : 5
            rejected_ECG = abs(ECG_avg) > std(ECG_avg(~rejected_ECG))*3;
        end
        
        % plot correlations
        fig = figure;
        ah = axes('parent', fig);
        hold(ah, 'on');
        for ibar = 1 : numel(ECG_avg)
            bar(ibar, ECG_avg(ibar), 'parent', ah, 'facecolor', [1 0 0] .* rejected_ECG(ibar));
        end
        title('average correlation component with ECG');
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_corr_comp_ECG_rs.pdf']);
        
        % plot removed ECG components
        cfg = [];
        cfg.component = find(rejected_ECG);
        cfg.layout = 'neuromag306mag.lay';
        cfg.marker = 'off';
        cfg.comment = 'no';
        fig = figure; ft_topoplotIC(cfg,comp);
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_topo_comp_ECG_rs.pdf']);
        
        % recombine remaining components
        cfg           = [];
        cfg.component = [find(rejected_EOG) find(rejected_ECG)];
        data_ICA      = ft_rejectcomponent(cfg,comp,data_epoched);

    else
        
        for ipart = 1 : 4
            
            % partial rejection
            cfg = [];
            cfg.artfctdef           = artdef_MEG{ipart};
            cfg.artfctdef.reject    = 'partial';
            data_artefact_rejected  = ft_rejectartifact(cfg,data_epoched{ipart});
            
            % determine rank
            cfg = [];
            cfg.trials              = 1;
            cfg.channel             = 'MEG';
            temp                    = ft_selectdata(cfg,data_epoched{ipart});
            data_rank               = rank(temp.trial{1}*temp.trial{1}');
            clear temp
            
            % decompose
            cfg                     = [];
            cfg.method              = 'runica';
            cfg.runica.pca          = data_rank;
            cfg.channel             = 'MEG';
            cfg.randomseed          = WANDER_plantseeds(isubject,0);
            comp                    = ft_componentanalysis(cfg, data_artefact_rejected);
            
            EOGvi = find(ismember(data_artefact_rejected.label,'BIO002'));
            if isempty(EOGvi)
                EOGvi = find(ismember(data_artefact_rejected.label,'BIO012'));
            end
            ECGi  = find(ismember(data_artefact_rejected.label,'BIO003'));
            
            % correlate components with EOG & ECG
            for itrial = 1 : size(data_artefact_rejected.trial,2)
                for icomp = 1 : size(comp.label,1)
                    fprintf('Correlating trial %d of %d\n',itrial,size(data_artefact_rejected.trial,2));
                    %                 EOGh(itrial,icomp) = corr(data_artefact_rejected.trial{itrial}(EOGhi,:)',comp.trial{itrial}(icomp,:)');
                    EOGv(itrial,icomp) = corr(data_artefact_rejected.trial{itrial}(EOGvi,:)',comp.trial{itrial}(icomp,:)');
                    ECG(itrial,icomp)  = corr(data_artefact_rejected.trial{itrial}(ECGi,:)',comp.trial{itrial}(icomp,:)');
                end
            end
            EOGv_avg  = mean(EOGv,1);
            rejected_EOG = false(1,data_rank);
            
            % determine to-be-rejected components based on thresold, iteratively
            for iteration = 1 : 5
                rejected_EOG = abs(EOGv_avg) > std(EOGv_avg(~rejected_EOG))*3;
            end
            
            % plot correlations
            fig = figure;
            ah = axes('parent', fig);
            hold(ah, 'on');
            for ibar = 1 : numel(EOGv_avg)
                bar(ibar, EOGv_avg(ibar), 'parent', ah, 'facecolor', [1 0 0] .* rejected_EOG(ibar));
            end
            title('average correlation component with vertical EOG');
            print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_corr_comp_EOGv_part ' num2str(ipart) '.pdf']);
            
            % plot removed EOGv components
            cfg = [];
            cfg.component = find(rejected_EOG);
            cfg.layout = 'neuromag306mag.lay';
            cfg.marker = 'off';
            cfg.comment = 'no';
            fig = figure; ft_topoplotIC(cfg,comp);
            print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_topo_comp_EOGv_part ' num2str(ipart) '.pdf']);
            
            % determine to-be-rejected components based on threshold, iteratively
            ECG_avg  = mean(ECG,1);
            rejected_ECG = false(1,data_rank);
            for iteration = 1 : 5
                rejected_ECG = abs(ECG_avg) > std(ECG_avg(~rejected_ECG))*3;
            end
            
            % plot correlations
            fig = figure;
            ah = axes('parent', fig);
            hold(ah, 'on');
            for ibar = 1 : numel(ECG_avg)
                bar(ibar, ECG_avg(ibar), 'parent', ah, 'facecolor', [1 0 0] .* rejected_ECG(ibar));
            end
            title('average correlation component with ECG');
            print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_corr_comp_ECG_part ' num2str(ipart) '.pdf']);
            
            % plot removed ECG components
            cfg = [];
            cfg.component = find(rejected_ECG);
            cfg.layout = 'neuromag306mag.lay';
            cfg.marker = 'off';
            cfg.comment = 'no';
            fig = figure; ft_topoplotIC(cfg,comp);
            print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_topo_comp_ECG_part ' num2str(ipart) '.pdf']);
            
            % recombine remaining components
            cfg           = [];
            cfg.component = [find(rejected_EOG) find(rejected_ECG)];
            data_ICA      = ft_rejectcomponent(cfg,comp,data_epoched);
            data_ICA{ipart} = ft_rejectcomponent(cfg,comp,data_epoched{ipart});
        
            close all
        end
    end
    
    % save data
    save(fname_ICA,'data_ICA','-v7.3');
end


