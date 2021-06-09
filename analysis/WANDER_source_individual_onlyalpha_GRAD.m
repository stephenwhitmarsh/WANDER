function [FFT_avg, FFT_avg_high, FFT_avg_low, source, source_low, source_high] = WANDER_source_individual_onlyalpha_GRAD(isubject,force,rootpath)

disp(['Working on Subjectnr. ' num2str(isubject)]);

if rootpath == 1
    fname_source_alpha = ['D:\WANDER\data\beamformer\s' num2str(isubject) '_source_individual_alpha_GRAD.mat'];
else
    fname_source_alpha = ['/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/WANDER/data/beamformer/s' num2str(isubject) '_source_individual_alpha_GRAD.mat'];
end

if exist(fname_source_alpha,'file') && force ~= 1
    fprintf('Returning TFR-EGG phaselocking \n');
    load(fname_source_alpha);
else
    fprintf('TFR-EGG phaselocking not found, creating it now! \n');
    
    % load data
    data_MEG            = WANDER_ICA(isubject,0,rootpath,0);

    for iblock = 1 : 4
        
        % select all channels, only correct rejections
        cfg                 = [];
        cfg.channel         = {'MEG*2','MEG*3'};
        cfg.trials          = find(data_MEG{iblock}.trialinfo(:,3) == 4);
        data_MEG{iblock}    = ft_selectdata(cfg,data_MEG{iblock});
        data_MEG{iblock}    = rmfield(data_MEG{iblock},'cfg');
      
        %% ADDED THE FOLLOWING TO MAKE TRIAL LENGTH CONSISTENT OVER MEASURES - 25-01-2018
        
        % start trial at first stimulation
        for itrial = 1 : size(data_MEG{iblock}.trial,2)
            fprintf('Cutting off beginning of trial %d block %d \n',itrial,iblock);
            data_MEG{iblock}.trial{itrial}          = data_MEG{iblock}.trial{itrial}(:,1500:end);
            data_MEG{iblock}.time{itrial}           = data_MEG{iblock}.time{itrial}(1500:end);
            data_MEG{iblock}.sampleinfo(itrial,1)   = data_MEG{iblock}.sampleinfo(itrial,1) + 1499;   
        end
%          
%         size(data_MEG{iblock}.trial{1})
%         data_MEG{iblock}.sampleinfo(1,1) - data_MEG{iblock}.sampleinfo(1,2)
%         
        % lock data to end of trial
        for itrial = 1 : size(data_MEG{iblock}.trial,2)
            fprintf('Shifting timecourse of trial %d block %d \n',itrial,iblock);
            data_MEG{iblock}.time{itrial}  = data_MEG{iblock}.time{itrial} - data_MEG{iblock}.time{itrial}(end) + 1.000;           
        end
        
        % select last 10 seconds of trial
        cfg = [];
        cfg.latency = [-10, 0];
        data_MEG{iblock} = ft_selectdata(cfg,data_MEG{iblock});    

        % SWITCHED TO ABOVE ON 29-1-2018
        
%         % cut off end of trial (end at last stim)
%         for itrial = 1 : size(data_MEG{iblock}.trial,2)
%             fprintf('Cutting off end of trial %d block %d \n',itrial,iblock);
%             triallength = size(data_MEG{iblock}.trial{itrial} ,2);
%             data_MEG{iblock}.trial{itrial}          = data_MEG{iblock}.trial{itrial}(:,1:triallength-1000);
%             data_MEG{iblock}.time{itrial}           = data_MEG{iblock}.time{itrial}(1:triallength-1000);
%             data_MEG{iblock}.sampleinfo(itrial,2)   = data_MEG{iblock}.sampleinfo(itrial,2) - 1000;   
%         end
%         
%         size(data_MEG{iblock}.trial{1})
%         data_MEG{iblock}.sampleinfo(1,1) - data_MEG{iblock}.sampleinfo(1,2)   
          % ADDED THE PRECEDING TO MAKE TRIAL LENGTH CONSISTENT OVER MEASURES - 25-01-2018
       
        % remove artefactual segments
        artdef                          = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);        
        cfg                             = [];
        cfg.reject                      = 'partial'; 
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
    
    % load alpha freq
    FFT_MEG = WANDER_FFT_MEG_alpha(isubject,0,'cue',rootpath,0);
    
    % freq analysis
%     cfg                 = [];
%     cfg.foi             = FFT_MEG.max_freq;
%     cfg.output          = 'powandcsd';
%     cfg.method          = 'mtmfft';
%     cfg.channel         = {'MEG*2', 'MEG*3'};
%     cfg.keeptrials      = 'yes'; 
%     cfg.tapsmofrq       = 0.5;
   
% will now just try with same freq range for each subject
    cfg                 = [];
    cfg.foi             = 11;
    cfg.output          = 'powandcsd';
    cfg.method          = 'mtmfft';
    cfg.channel         = {'MEG*2', 'MEG*3'};
    cfg.keeptrials      = 'yes'; 
    cfg.tapsmofrq       = 3;    
    
    
    
    % tried at 14-11-2018
%     cfg.tapsmofrq       = 1;   
    cfg.pad             = 'nextpow2';
    FFT                 = ft_freqanalysis(cfg, data_MEG_append);
    
    % median split based on ratings of 'distraction', so low = high
    % attention and vice-versa... yeah, stupid, but didn't want to change
    % half-way.
    F                   = ceil(2 * tiedrank(data_MEG_append.trialinfo(:,2)) / length(data_MEG_append.trialinfo(:,2)));  
    rating_split        = ones(size(F)); 
    rating_split(F==2)  = 2;
    
    cfg                 = [];
    cfg.trials          = find(rating_split == 1);
    FFT_low             = ft_selectdata(cfg,FFT);
    cfg.trials          = find(rating_split == 2);
    FFT_high            = ft_selectdata(cfg,FFT);
    
    % source analysis
    cfg                 = [];
    cfg.channel         = {'MEG*2', 'MEG*3'};
    cfg.grad            = hdr.grad;
    cfg.headmodel       = headmodel;
    cfg.keeptrials      = 'no';
    cfg.rawtrial        = 'no';
    cfg.method          = 'dics';
    cfg.grid            = leadfield;
    cfg.grid.filter     = common_filter; 
    source              = ft_sourceanalysis(cfg, FFT);
    source_low          = ft_sourceanalysis(cfg, FFT_low);
    source_high         = ft_sourceanalysis(cfg, FFT_high);

    source              = rmfield(source,'cfg');
    source_low          = rmfield(source_low,'cfg');
    source_high         = rmfield(source_high,'cfg');
    
    % average over trials
    cfg = [];
    cfg.avgoverrpt      = 'yes';
    FFT_avg             = ft_selectdata(cfg,FFT);
    FFT_avg_high        = ft_selectdata(cfg,FFT_high);
    FFT_avg_low         = ft_selectdata(cfg,FFT_low);
    
    % save data
    fprintf('Saving data in %s. This can take a while.\n',fname_source_alpha);
    save(fname_source_alpha,'FFT_avg','FFT_avg_high','FFT_avg_low','source','source_low','source_high','-v7.3');
end