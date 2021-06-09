function data_EGG_epoched = WANDER_ICA_probelocked(isubject,force,timing)

fname_ICA_probelocked    = ['d:\analysis\WANDER\data\ICA\s' num2str(isubject) '_epoched_MEG_ICA_' timing '.mat'];

if exist(fname_ICA_probelocked,'file') && force ~= 1
    fprintf('Returning ICA cleaned data\n');
    data_EGG_epoched = load(fname_ICA_probelocked);
else
    fprintf('ICA cleaned data not found, creating it now! \n');
    addpath('D:/analysis/WANDER/scripts/');
    WANDER_subjectinfo;
    
    % load epoched data
    temp = WANDER_epoch_MEG(isubject,0,timing);
    data_epoched = temp.data_MEG;
    clear temp
    
    % load artefact definition
    artdef_MEG = WANDER_artefact_detection_MEG(isubject,0,timing);
        
    for ipart = 1 : 4
        
        % partial rejection
        cfg = [];
        cfg.artfctdef           = artdef_MEG.artdef{ipart};
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
        comp                    = ft_componentanalysis(cfg, data_artefact_rejected);
        
        EOGhi = find(ismember(data_artefact_rejected.label,'BIO001'));
        ECGi  = find(ismember(data_artefact_rejected.label,'BIO003'));
        EOGvi = find(ismember(data_artefact_rejected.label,'BIO002'));
        if isempty(EOGvi)
            EOGvi = find(ismember(data_artefact_rejected.label,'BIO012'));
        end

        % correlate components with EOG 
        for itrial = 1 : size(data_artefact_rejected.trial,2)
            for icomp = 1 : size(comp.label,1)
                fprintf('Correlating trial %d\n',itrial);
                EOGh(itrial,icomp) = corr(data_artefact_rejected.trial{itrial}(EOGhi,:)',comp.trial{itrial}(icomp,:)');
                EOGv(itrial,icomp) = corr(data_artefact_rejected.trial{itrial}(EOGvi,:)',comp.trial{itrial}(icomp,:)');
                ECG(itrial,icomp)  = corr(data_artefact_rejected.trial{itrial}(ECGi,:)',comp.trial{itrial}(icomp,:)');
            end
        end

        % determine EOGh components based on thresold 
        close all
        EOGh = mean(EOGh,1);
        EOGhb = abs(EOGh) > std(EOGh)*2; find(EOGhb)
        fig = figure;
        ah = axes('parent', fig);
        hold(ah, 'on');
        for ibar = 1 : numel(EOG1)
            bar(ibar, EOGh(ibar), 'parent', ah, 'facecolor', [1 0 0] .* EOGhb(ibar));
        end
        title('average correlation component with horizontal EOG');
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_part' num2str(ipart) '_corr_comp_EOGh_' timing '.pdf']);
        
        % plot removed EOGh components
        cfg = [];
        cfg.component = find(EOGhb);
        cfg.layout = 'neuromag306mag.lay';
        cfg.marker = 'off';
        cfg.comment = 'no';
        fig = figure; ft_topoplotIC(cfg,comp);
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_part' num2str(ipart) '_topo_comp_EOGh_' timing '.pdf']);
        
        % determine EOGv components based on thresold
        EOGv  = mean(EOGv,1);
        EOGvb = abs(EOGv) > std(EOGv)*2; find(EOGvb)
        fig = figure;
        ah = axes('parent', fig);
        hold(ah, 'on');
        for ibar = 1 : numel(EOGv)
            bar(ibar, EOGv(ibar), 'parent', ah, 'facecolor', [1 0 0] .* EOGvb(ibar));
        end
        title('average correlation component with vertical EOG');
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_part' num2str(ipart) '_corr_comp_EOGv_' timing '.pdf']);
        
        % plot removed EOGv components
        cfg = [];
        cfg.component = find(EOGvb);
        cfg.layout = 'neuromag306mag.lay';
        cfg.marker = 'off';
        cfg.comment = 'no';
        fig = figure; ft_topoplotIC(cfg,comp);
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_part' num2str(ipart) '_topo_comp_EOGv_' timing '.pdf']);
        
        % determine ECG components based on thresold
        ECG  = mean(ECG,1);
        ECGb = abs(ECG) > std(ECG)*2; find(ECGb)
        fig = figure;
        ah = axes('parent', fig);
        hold(ah, 'on');
        for ibar = 1 : numel(ECG)
            bar(ibar, ECG(ibar), 'parent', ah, 'facecolor', [1 0 0] .* ECGb(ibar));
        end
        title('average correlation component with ECG');
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_part' num2str(ipart) '_corr_comp_ECG_' timing '.pdf']);
        
        % plot removed components
        cfg = [];
        cfg.component = find(ECGb);
        cfg.layout = 'neuromag306mag.lay';
        cfg.marker = 'off';
        cfg.comment = 'no';
        figure; ft_topoplotIC(cfg,comp);
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_part' num2str(ipart) '_topo_comp_ECG_' timing '.pdf']);
        
        % recombine remaining components
        cfg             = [];
        cfg.component   = [find(EOGhb) find(EOGvb) find(ECGb)];
        data_ICA{ipart} = ft_rejectcomponent(cfg,comp,data_epoched{ipart});
         
%         % browse components
%         cfg = [];
%         cfg.viewmode    = 'component';
%         cfg.component   = 1:10;
%         cfg.layout      = 'neuromag306mag.lay';
%         ft_databrowser(cfg,comp);
%         
%         % select components to remove
%         badcomp         = input('What MEG components do you want removed? e.g.: [1 2:5 9]');
%         
%         % recombine remaining components
%         cfg             = [];
%         cfg.component   = badcomp;
%         cfg.component   = [find(EOG1b) find(EOG2b) find(ECGb)];
%         data_ICA{ipart} = ft_rejectcomponent(cfg,comp,data_epoched{ipart});


    end
    
    % save data
    save(fname_ICA_probelocked,'data_ICA','-v7.3');
end


