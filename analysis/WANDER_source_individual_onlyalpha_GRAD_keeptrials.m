function [source] = WANDER_source_individual_onlyalpha_GRAD_keeptrials(isubject,force,rootpath)

disp(['Working on Subjectnr. ' num2str(isubject)]);

if rootpath == 1
    fname_source_alpha = ['z:\WANDER\data\beamformer\s' num2str(isubject) '_source_individual_alpha_GRAD_keeptrials.mat'];
else
    fname_source_alpha = ['/shared/projects/project_wander/WANDER/data/beamformer/s' num2str(isubject) '_source_individual_alpha_GRAD_keeptrials.mat'];
end

if exist(fname_source_alpha,'file') && force ~= 1
    fprintf('Returning TFR-EGG phaselocking \n');
    load(fname_source_alpha);
else
    fprintf('TFR-EGG phaselocking not found, creating it now! \n');
    
    % load data
    data_MEG            = WANDER_ICA(isubject,0,rootpath,0);
    artdef              = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
    
    for iblock = 1 : 4
        
        % select GRAD channels
        cfg                             = [];
        cfg.channel                     = {'MEG*2','MEG*3'};
        data_MEG{iblock}                = ft_selectdata(cfg,data_MEG{iblock});
        data_MEG{iblock}                = rmfield(data_MEG{iblock},'cfg');
        
        % add trialnum and sampleinfo
        for itrial = 1 : size(data_MEG{iblock}.trialinfo,1)
            data_MEG{iblock}.trialinfo(itrial,7) = itrial + (iblock-1)*50;
            data_MEG{iblock}.trialinfo(itrial,8) = data_MEG{iblock}.sampleinfo(itrial,1) - data_MEG{iblock}.sampleinfo(1,1);
            data_MEG{iblock}.trialinfo(itrial,9) = data_MEG{iblock}.sampleinfo(itrial,2) - data_MEG{iblock}.sampleinfo(1,1);
        end
        
        %         % end trial at last stimulation
        %         for itrial = 1 : size(data_MEG{iblock}.trial,2)
        %             fprintf('Cutting off end of trial %d block %d \n',itrial,iblock);
        %             triallength = size(data_MEG{iblock}.trial{itrial},2);
        %             data_MEG{iblock}.trial{itrial} = data_MEG{iblock}.trial{itrial}(:,1:triallength-1000);
        %             data_MEG{iblock}.time{itrial}  = data_MEG{iblock}.time{itrial}(1:triallength-1000);
        %         end
        %
        %         % start trial at first stimulation
        %         for itrial = 1 : size(data_MEG{iblock}.trial,2)
        %             fprintf('Cutting off beginning of trial %d block %d \n',itrial,iblock);
        %             triallength = size(data_MEG{iblock}.trial{itrial},2);
        %             data_MEG{iblock}.trial{itrial} = data_MEG{iblock}.trial{itrial}(:,triallength-10000:end);
        %             data_MEG{iblock}.time{itrial}  = data_MEG{iblock}.time{itrial}(triallength-10000:end);
        %         end
        %
        %         % remove sampleinfo which else creates unexpected behaviour in
        %         % rejectartefact
        %         data_MEG{iblock}      (isubject,force,rootpath,restingstate)          = rmfield(data_MEG{iblock},'sampleinfo');
        %
        %         % remove artefacts
        %         cfg                             = [];
        %         cfg.artfctdef                   = artdef{iblock};
        %         cfg.artfctdef.reject            = 'complete';
        %         data_MEG_clean{iblock}          = ft_rejectartifact(cfg,data_MEG{iblock});
        %
        
        % start trial at first stimulation
        for itrial = 1 : size(data_MEG{iblock}.trial,2)
            fprintf('Cutting off beginning of trial %d block %d \n',itrial,iblock);
            data_MEG{iblock}.trial{itrial}          = data_MEG{iblock}.trial{itrial}(:,1500:end);
            data_MEG{iblock}.time{itrial}           = data_MEG{iblock}.time{itrial}(1500:end);
            data_MEG{iblock}.sampleinfo(itrial,1)   = data_MEG{iblock}.sampleinfo(itrial,1) + 1499;
        end
        
        % lock data to end of trial
        for itrial = 1 : size(data_MEG{iblock}.trial,2)
            fprintf('Shifting timecourse of trial %d block %d \n',itrial,iblock);
            data_MEG{iblock}.time{itrial}  = data_MEG{iblock}.time{itrial} - data_MEG{iblock}.time{itrial}(end) + 1.000;
        end
        
        % select last 10 seconds of trial + second after offset
        cfg = [];
        cfg.latency = [-10, 0];
        data_MEG{iblock} = ft_selectdata(cfg,data_MEG{iblock});
               
        % remove artefactual segments
        cfg                             = [];
        cfg.reject                      = 'complete';
        cfg.artfctdef                   = artdef{iblock};
        cfg.artfctdef.minaccepttim      = 3;
        data_MEG_clean{iblock}          = ft_rejectartifact(cfg,data_MEG{iblock});
        
    end
    
    clear data_MEG
    data_MEG_append = ft_appenddata([],data_MEG_clean{:});
    clear data_MEG_clean
    
    % load headmodel
    [~, headmodel, ~, ~, ~] = WANDER_grid(isubject,rootpath,0);
    
    % load sensor information
    [dataset_task, ~] = WANDER_subjectinfo(rootpath);
    hdr = ft_read_header(dataset_task{isubject,1});
    
    % load common filter
    [common_filter, leadfield] = WANDER_common_filter_DICS_individual_onlyalpha(isubject,rootpath,0);
    
    % get individual alpha
%     if isubject ~= 13
        [FFT_MEG] = WANDER_FFT_MEG_alpha(isubject,0,'cue',rootpath,0);
        foi = FFT_MEG.max_freq;
%     else
%         foi = 11;
%     end
    
    cfg = [];
    cfg.foi             = foi;
    cfg.output          = 'powandcsd';
    cfg.method          = 'mtmfft';
    cfg.channel         = {'MEG*2', 'MEG*3'};
    cfg.keeptrials      = 'yes';
    cfg.tapsmofrq       = 0.5;
    cfg.pad             = 'nextpow2';
    FFT                 = ft_freqanalysis(cfg, data_MEG_append);
    FFT                 = rmfield(FFT,'cfg');
    
    cfg                 = [];
    cfg.channel         = {'MEG*2', 'MEG*3'};
    cfg.grad            = hdr.grad;
    cfg.headmodel       = headmodel;
    cfg.keeptrials      = 'yes';
    cfg.rawtrial        = 'yes';
    cfg.method          = 'dics';
    cfg.grid            = leadfield;
    cfg.grid.filter     = common_filter;
    source              = ft_sourceanalysis(cfg, FFT);
    
    % save data
    fprintf('Saving data in %s. This can take a while.\n',fname_source_alpha);
    save(fname_source_alpha,'source','-v7.3');
end