function [TFR] = WANDER_TFR_probelocked(isubject,force,timing)

% if ~strcmp(timing,'cue') || ~strcmp(timing,'probe')
%     fprintf('Use "cue" or "probe" as third argument.\n');
%     return
% end

WANDER_subjectinfo;

fname_TFR_probelocked = ['d:\analysis\WANDER\data\TFR\s' num2str(isubject) '_TFR_probelocked.mat'];

if exist(fname_TFR_probelocked,'file') && force ~= 1
    fprintf('Returning TFR\n');
    load(fname_TFR_probelocked);
else
    fprintf('TFR not found, creating it now! \n');
    %     data_EGG_epoched = WANDER_ICA_probelocked(isubject);
    data_epoch_MEG = WANDER_epoch_MEG(isubject,0,timing);
    data_epoch_MEG = data_epoch_MEG.data_MEG;
    
    for ipart = 1 : 4
        % TFR
        cfg              = [];
        cfg.output       = 'pow';
        cfg.channel      = 'all';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.keeptrials   = 'yes';
        cfg.channel      = {'MEG'};
        cfg.foi          = 1:1:30;
        if strcmp(timing,'cue');
            cfg.toi      = -2:0.05:10;
        end
        if strcmp(timing,'probe');
            cfg.toi      = -30:0.05:1;
        end
        cfg.t_ftimwin    = ones(size(cfg.foi));
        TFR{ipart}       = ft_freqanalysis(cfg, data_epoch_MEG{ipart});
    end
    
    save(fname_TFR_probelocked,'TFR','-v7.3');
end