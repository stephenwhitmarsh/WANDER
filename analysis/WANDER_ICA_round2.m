function data_EGG_epoched = WANDER_ICA_round2(isubject,force,rootpath)

if isempty(rootpath)
    rootpath = 'd:\analysis\WANDER\data\';
end
if rootpath == 1
    fname_ICA_round2  = ['d:\analysis\WANDER\data\ICA\s' num2str(isubject) '_epoched_MEG_ICA_round2.mat'];
else
    fname_ICA_round2  = ['/home/swhitmarsh/WANDER/data/ICA/s' num2str(isubject) '_epoched_MEG_ICA_round2.mat'];
end

if exist(fname_ICA_round2,'file') && force ~= 1
    fprintf('Returning ICA cleaned data, round 2\n');
    data_EGG_epoched = load(fname_ICA_round2);
else
    fprintf('ICA cleaned data, round 2, not found, creating it now! \n');
    WANDER_subjectinfo;
    
    % load epoched data
    data_epoched = WANDER_ICA(isubject,0,rootpath);
    
    % load artefact definition
    
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
        ECGi  = find(ismember(data_artefact_rejected.label,'BIO003'));
        if isempty(EOGvi)
            EOGvi = find(ismember(data_artefact_rejected.label,'BIO012'));
        end
        
        % correlate components with EOG
        for itrial = 1 : size(data_artefact_rejected.trial,2)
            for icomp = 1 : size(comp.label,1)
                fprintf('Correlating trial %d\n',itrial);
                EOGv(itrial,icomp) = corr(data_artefact_rejected.trial{itrial}(EOGvi,:)',comp.trial{itrial}(icomp,:)');
                ECG(itrial,icomp)  = corr(data_artefact_rejected.trial{itrial}(ECGi,:)',comp.trial{itrial}(icomp,:)');
            end
        end
        
        % determine EOGv components based on thresold
        EOGv  = mean(EOGv,1);
        EOGvb = abs(EOGv) > std(EOGv)*3; find(EOGvb)
        fig = figure;
        ah = axes('parent', fig);
        hold(ah, 'on');
        for ibar = 1 : numel(EOGv)
            bar(ibar, EOGv(ibar), 'parent', ah, 'facecolor', [1 0 0] .* EOGvb(ibar));
        end
        title('average correlation component with vertical EOG');
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_corr_comp_EOGv_part' num2str(ipart) '_round2.pdf']);
        
        % plot removed EOGv components
        cfg = [];
        cfg.component = find(EOGvb);
        cfg.layout = 'neuromag306mag.lay';
        cfg.marker = 'off';
        cfg.comment = 'no';
        fig = figure; ft_topoplotIC(cfg,comp);
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_topo_comp_EOGv_part' num2str(ipart) '_round2.pdf']);
        
        % determine ECG components based on thresold
        ECG  = mean(ECG,1);
        ECGb = abs(ECG) > std(ECG)*3; find(ECGb)
        fig = figure;
        ah = axes('parent', fig);
        hold(ah, 'on');
        for ibar = 1 : numel(ECG)
            bar(ibar, ECG(ibar), 'parent', ah, 'facecolor', [1 0 0] .* ECGb(ibar));
        end
        title('average correlation component with ECG');
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_corr_comp_ECG_part' num2str(ipart) '_round2.pdf']);
        
        % plot removed components
        cfg = [];
        cfg.component = find(ECGb);
        cfg.layout = 'neuromag306mag.lay';
        cfg.marker = 'off';
        cfg.comment = 'no';
        fig = figure; ft_topoplotIC(cfg,comp);
        print(fig,'-painters','-dpdf','-r100',['d:\analysis\WANDER\images\ICA\s' num2str(isubject) '_topo_comp_ECG_part' num2str(ipart) '_round2.pdf']);
        
        % recombine remaining components
        cfg             = [];
        cfg.component   = [find(EOGvb) find(ECGb)];
        data_ICA{ipart} = ft_rejectcomponent(cfg,comp,data_epoched{ipart});
        
        close all
    end
    
    % save data
    save(fname_ICA_round2,'data_ICA','-v7.3');
end


