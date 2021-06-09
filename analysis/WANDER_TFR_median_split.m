function [TFR_corr] = WANDER_TFR_median_split(isubject,force,timing)

if ~(strcmp(timing,'cue') || strcmp(timing,'probe'))
    fprintf('Use "cue" or "probe" as third argument.\n');
    return
end

WANDER_subjectinfo;

fname_TFR_corr = ['h:\analysis\WANDER\data\TFR\s' num2str(isubject) '_TFR_corr_' timing '.mat'];

if exist(fname_TFR_corr,'file') && force ~= 1
    fprintf('Returning TFR correlation \n');
    load(fname_TFR_corr);
else
    fprintf('TFR correlation not found, creating it now! \n');   
    TFR = WANDER_TFR(isubject,0,timing);
    
    % only combined gradiometers - will eventually run with magnetometers as well
    % only correct rejections
    % only selected timesegment
    for ipart = 1 : 4
        TFR{ipart}  = ft_combineplanar([],TFR{ipart});
        cfg         = [];
        cfg.channel = 'MEG*3';
        cfg.trials  = TFR{ipart}.trialinfo(:,3) == 4;
        cfg.latency = [-30 0];
        TFR{ipart} = ft_selectdata(cfg,TFR{ipart});
    end
    
    % combine blocks.
    % NOTE: from this point on gradiometer definition is invalid!!!
    TFR_comb             = TFR{1};
    TFR_comb.powspctrm   = [TFR{1}.powspctrm; TFR{2}.powspctrm; TFR{3}.powspctrm; TFR{4}.powspctrm];
    TFR_comb.trialinfo   = [TFR{1}.trialinfo; TFR{2}.trialinfo; TFR{3}.trialinfo; TFR{4}.trialinfo];
    TFR_comb             = rmfield(TFR_comb,'cumtapcnt');
    TFR_comb             = rmfield(TFR_comb,'grad');
%     TFR_corr.dimord      = 'chan_freq_time';
%     TFR_corr             = rmfield(TFR_corr,'powspctrm');
    clear TFR  
    
    % load artifact definition and reshape to match TFR
    fname_TFR_phaselocking = ['h:\analysis\WANDER\data\TFR\s' num2str(isubject) '_TFR_phaselocked_' timing '.mat'];
    load(fname_TFR_phaselocking,'artefact_mask')
    artefact_mask = reshape(artefact_mask,size(TFR_comb.powspctrm,1), 1, 1, size(TFR_comb.powspctrm,4));
    
    % remove measures during artefact periods
    TFR_comb.powspctrm(repmat(artefact_mask,1,size(TFR_comb.powspctrm,2),size(TFR_comb.powspctrm,3),1)) = NaN;    
    
    fprintf('Calculating correlation subject %d. This will take a while! \n',isubject);
    TFR_corr = TFR_comb;
    TFR_dimord = 'chan_freq_time';
    for ichan = 1 : size(TFR_comb.powspctrm,2)
        fprintf('Channel %d of %d \n',ichan,size(TFR_comb.powspctrm,2));
        for ifreq = 1 : size(TFR_comb.powspctrm,3)
            for itime = 1 : size(TFR_comb.powspctrm,4)
                [TFR_corr.rho(ichan,ifreq,itime), TFR_corr.pval(ichan,ifreq,itime)] = corr(TFR_comb.powspctrm(:,ichan,ifreq,itime),TFR_comb.trialinfo(:,2),'type','spearman','rows','pairwise');
            end
        end
    end
    
%     TFR_corr.rho    = reshape(rho,nrchan,nrfreq,nrtime);
%     TFR_corr.pval   = reshape(pval,nrchan,nrfreq,nrtime);
    
    for itime = 1:size(TFR_corr.rho,3)
        TFR_corr.nr(itime) = sum(~isnan(TFR_comb.powspctrm(:,1,1,itime)));
    end
    save(fname_TFR_corr,'TFR_corr','-v7.3');
end
% 
% % TFR multiplot average
% fig             = figure; 
% cfg             = [];
% cfg.parameter   = 'pval';
% cfg.layout      = 'neuromag306cmb';
% cfg.zlim        = 'maxabs';
% ft_singleplotTFR(cfg,TFR_corr);
% hold;
% ax = axis;
% perc = squeeze(TFR_corr.nr) ./ max(TFR_corr.nr) * ax(4);
% plot(TFR_corr.time,perc,'r','linewidth',2);
% print(fig,'-dpdf',['d:\analysis\WANDER\images\TFR\s' num2str(isubject) '_' timing '_TFR_corr_cue.pdf']);
% print(fig,'-dpng',['d:\analysis\WANDER\images\TFR\s' num2str(isubject) '_' timing '_TFR_corr_cue.png']);


