function WANDER_SSEF_GA_topo_and_dipoles(rootpath, force)

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 3 SD
rootpath = 1;

for isubject = slist
    disp(['subject ', num2str(isubject)]);
    
    if rootpath == 1
        fname = ['d:\WANDER\data\SSEF\SSEF_phase_' num2str(isubject) '.mat'];
    else
        fname = ['/shared/projects/project_wander/WANDER/data/SSEF/SSEF_phase_' num2str(isubject) '.mat'];
    end

    % load directly from file
    temp = load(fname,'data_MEG_stim_avg');
    data_MEG_stim_avg{isubject}             = temp.data_MEG_stim_avg;
    temp = load(fname,'data_MEG_stim_avg_high');
    data_MEG_stim_avg_high{isubject}        = temp.data_MEG_stim_avg_high;
    temp = load(fname,'data_MEG_stim_avg_low');
    data_MEG_stim_avg_low{isubject}         = temp.data_MEG_stim_avg_low;
    %     temp = load(fname_SSEF_phase,'data_MEG_stim_MI');
    %     data_MEG_stim_MI{isubject}              = temp.data_MEG_stim_MI;
    %     temp = load(fname_SSEF_phase,'data_MEG_stim_MI_high');
    %     data_MEG_stim_MI_high{isubject}         = temp.data_MEG_stim_MI_high;
    %     temp = load(fname_SSEF_phase,'data_MEG_stim_MI_low');
    %     data_MEG_stim_MI_low{isubject}          = temp.data_MEG_stim_MI_low;
    %
    %     data_MEG_stim_MI_low{isubject}.fsample  = 1000;
    %     data_MEG_stim_MI_high{isubject}.fsample = 1000;
    %     data_MEG_stim_MI{isubject}.fsample      = 1000;
    %
    % for some reason it reverts to single trial rather than average data
    % structure - this is a fix
    %     data_MEG_stim_MI_low{isubject}          = ft_timelockanalysis([],data_MEG_stim_MI_low{isubject});
    %     data_MEG_stim_MI_high{isubject}         = ft_timelockanalysis([],data_MEG_stim_MI_high{isubject});
    %     data_MEG_stim_MI{isubject}              = ft_timelockanalysis([],data_MEG_stim_MI{isubject});
    
    %     % normalize data as percentage of total
    %     data_MEG_MI_diff{isubject}              = data_MEG_stim_MI{isubject};
    %     data_MEG_MI_diff{isubject}.avg          = (data_MEG_stim_MI_high{isubject}.avg - data_MEG_stim_MI_low{isubject}.avg) ./ (data_MEG_stim_MI_high{isubject}.avg + data_MEG_stim_MI_low{isubject}.avg);
    %     %     data_MEG_MI_diff{isubject}.avg          = (data_MEG_stim_MI_high{isubject}.avg - data_MEG_stim_MI_low{isubject}.avg);
    %
    %     data_MEG_stim_MI_low_norm{isubject}     = data_MEG_stim_MI_low{isubject};
    %     data_MEG_stim_MI_high_norm{isubject}    = data_MEG_stim_MI_high{isubject};
    %     data_MEG_stim_MI_low_norm{isubject}.avg      = data_MEG_stim_MI_low{isubject}.avg  ./ (data_MEG_stim_MI_high{isubject}.avg + data_MEG_stim_MI_low{isubject}.avg);
    %     data_MEG_stim_MI_high_norm{isubject}.avg     = data_MEG_stim_MI_high{isubject}.avg ./ (data_MEG_stim_MI_high{isubject}.avg + data_MEG_stim_MI_low{isubject}.avg);
    %
    data_MEG_stim_diff{isubject}            = data_MEG_stim_avg{isubject};
    data_MEG_stim_diff{isubject}.avg        = (data_MEG_stim_avg_low{isubject}.avg - data_MEG_stim_avg_high{isubject}.avg);
    
end

% calculate rms and time of maximum response
for isubject = slist
    
    % square ERF
    data_MEG_stim_rms{isubject}          = data_MEG_stim_avg{isubject};
    data_MEG_stim_rms{isubject}.avg      = data_MEG_stim_rms{isubject}.avg .^2;
    data_MEG_stim_rms_high{isubject}     = data_MEG_stim_avg_high{isubject};
    data_MEG_stim_rms_high{isubject}.avg = data_MEG_stim_rms_high{isubject}.avg.^2;
    data_MEG_stim_rms_low{isubject}      = data_MEG_stim_avg_low{isubject};
    data_MEG_stim_rms_low{isubject}.avg  = data_MEG_stim_rms_low{isubject}.avg.^2;
    
    % root mean ERF
    cfg = [];
    cfg.avgovertime = 'no';
    data_MEG_stim_rms{isubject}          = ft_selectdata(cfg,data_MEG_stim_rms{isubject});
    data_MEG_stim_rms{isubject}.avg      = sqrt(data_MEG_stim_rms{isubject}.avg);
    data_MEG_stim_rms_high{isubject}     = ft_selectdata(cfg,data_MEG_stim_rms_high{isubject});
    data_MEG_stim_rms_high{isubject}.avg = sqrt(data_MEG_stim_rms_high{isubject}.avg);
    data_MEG_stim_rms_low{isubject}      = ft_selectdata(cfg,data_MEG_stim_rms_low{isubject});
    data_MEG_stim_rms_low{isubject}.avg  = sqrt(data_MEG_stim_rms_low{isubject}.avg);
    
    % diff
    data_MEG_stim_rms_diff{isubject}     = data_MEG_stim_rms_high{isubject};
    %     data_MEG_stim_rms_diff{isubject}.avg = (data_MEG_stim_rms_low{isubject}.avg - data_MEG_stim_rms_high{isubject}.avg) ./ (data_MEG_stim_rms_low{isubject}.avg + data_MEG_stim_rms_high{isubject}.avg) ;
    data_MEG_stim_rms_diff{isubject}.avg = (data_MEG_stim_rms_low{isubject}.avg - data_MEG_stim_rms_high{isubject}.avg);
    
end

i = 1;
figure;
for isubject = slist
    
    % find time of maximum response on RMS
    cfg = [];
    cfg.channel         = 'MEG*1';
    cfg.channel         = {'MEG0721', 'MEG0731', 'MEG0911', 'MEG0921', 'MEG0931', 'MEG0941', 'MEG1021', 'MEG1031', 'MEG1041', 'MEG1111', 'MEG1121', 'MEG1131', 'MEG1141', 'MEG1211', 'MEG1221', 'MEG1231', 'MEG1241', 'MEG1311', 'MEG1321', 'MEG1331', 'MEG1341', 'MEG1411', 'MEG1421', 'MEG1431', 'MEG1441', 'MEG2021', 'MEG2031', 'MEG2131', 'MEG2211', 'MEG2221', 'MEG2231', 'MEG2241', 'MEG2311', 'MEG2321', 'MEG2331', 'MEG2341', 'MEG2411', 'MEG2421', 'MEG2431', 'MEG2441', 'MEG2511', 'MEG2521', 'MEG2531', 'MEG2541', 'MEG2611', 'MEG2621', 'MEG2631', 'MEG2641'};
    cfg.avgoverchan     = 'yes';
    temp                = ft_selectdata(cfg,data_MEG_stim_avg{isubject});
    
    [pks_pos,locs_pos]  = findpeaks(temp.avg);
    [pks_neg,locs_neg]  = findpeaks(-temp.avg);
    pks                 = [pks_pos, pks_neg];
    locs                = [locs_pos, locs_neg];
    
    [amp,ipks]          = sort(abs(pks),'descend');
    
    max_time(isubject)  = temp.time(locs(ipks(1)));
    max_time(3)         = 0.037;
    max_time(8)         = 0.03;
    
    max_timecourse_rms(isubject,:) = temp.avg;
    
    % find maximum sensor on non-RMS data for max timepoint
    cfg = [];
    cfg.channel                         = 'MEG*1';
    cfg.avgoverchan                     = 'no';
    temp                                = ft_selectdata(cfg,data_MEG_stim_avg{isubject});
    temp_high                           = ft_selectdata(cfg,data_MEG_stim_avg_high{isubject});
    temp_low                            = ft_selectdata(cfg,data_MEG_stim_avg_low{isubject});
    
    [~, max_sensor(isubject)]           = max(temp.avg(:,locs(ipks(1))));
    [~, min_sensor(isubject)]           = min(temp.avg(:,locs(ipks(1))));
    
    max_label{isubject}                 = temp.label{max_sensor(isubject)};
    min_label{isubject}                 = temp.label{min_sensor(isubject)};
    
    % timecourse of maximum sensor
    max_timecourse_sel(isubject,:)      = temp.avg(max_sensor(isubject),:);
    min_timecourse_sel(isubject,:)      = temp.avg(min_sensor(isubject),:);
    max_timecourse_sel_high(isubject,:) = temp_high.avg(max_sensor(isubject),:);
    min_timecourse_sel_high(isubject,:) = temp_high.avg(min_sensor(isubject),:);
    max_timecourse_sel_low(isubject,:)  = temp_low.avg(max_sensor(isubject),:);
    min_timecourse_sel_low(isubject,:)  = temp_low.avg(min_sensor(isubject),:);
    
    
    subplot(25,1,i); hold;
    %     plot(max_timecourse(isubject,:)); axis tight;
    plot(max_timecourse_rms(isubject,:)); axis tight;
    
    i = i + 1;
end

%% flip orientations

flip_SS = round(rand(1,26));
flip_SS(flip_SS==0) = -1;
temp_SS = mean(max_timecourse_sel_high(slist,:) ,1);
for iteration = 1 : 100
    disp(['orienting iteration',num2str(iteration)]);
    for isubject = slist
        corr(max_timecourse_sel(isubject,:)',temp_SS')
        if corr(max_timecourse_sel(isubject,:)',temp_SS') < 0
            flip_SS(isubject)                   = -flip_SS(isubject);
            max_timecourse_sel(isubject,:)      = -max_timecourse_sel(isubject,:);
            max_timecourse_sel_high(isubject,:) = -max_timecourse_sel_high(isubject,:);
            max_timecourse_sel_low(isubject,:)  = -max_timecourse_sel_low(isubject,:);
            
        end
    end
    temp_SS = mean(max_timecourse_sel_high(slist,:) ,1);
    
end

%% plot single steady-state response

figure;
time = linspace(1,64,64);
subplot(3,1,1);
hist(max_time(slist)*1000,[0:1:63]); axis tight
title('max time');
xlabel('ms');
subplot(3,1,2); hold;

patch([time,time(end:-1:1)], [mean(max_timecourse_sel_high(slist,:)) + std(max_timecourse_sel_high(slist,:)), mean(max_timecourse_sel_high(slist,end:-1:1)) - std(max_timecourse_sel_high(slist,end:-1:1))],[1 0 0],'LineStyle','None','FaceAlpha',0.5);
patch([time,time(end:-1:1)], [mean(max_timecourse_sel_low(slist,:)) + std(max_timecourse_sel_low(slist,:)), mean(max_timecourse_sel_low(slist,end:-1:1)) - std(max_timecourse_sel_low(slist,end:-1:1))],[0 0 1],'LineStyle','None','FaceAlpha',0.5);
plot(mean(max_timecourse_sel_high(slist,:)),'r');
plot(mean(max_timecourse_sel_low(slist,:)),'b');
axis tight

subplot(3,1,3);
plot(max_timecourse_sel(slist,:)'); axis tight;
title('timecourse per subject');
xlabel('ms');

%% plot topo per subject
fig = figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    
    cfg = [];
    cfg.zlim                = 'maxabs';
    cfg.layout              = 'neuromag306mag';
    cfg.channel             = {'MEG*1'};
    cfg.highlight           = 'on';
    cfg.comment             = 'no';
    cfg.highlightchannel    = max_label{isubject};
    cfg.xlim                = [max_time(isubject) max_time(isubject)];
%     cfg.xlim                = [0.032 0.063];
    cfg.marker = 'off';
    ft_topoplotER(cfg,data_MEG_stim_avg{isubject});
    
    title(num2str(isubject));
    i = i + 1;
end
print(fig,'d:\WANDER\images\article\ERF_dipole_all_subjects.png','-r300','-dpng');

%% dipole fitting
for isubject = slist
    
    % load headmodel
    [mri_segmented{isubject}, headmodel{isubject}, subject_grid{isubject}, template_grid{isubject}, mri_realigned{isubject}] = WANDER_grid(isubject,rootpath,0);

    % load sensor information
    [dataset_task, ~] = WANDER_subjectinfo(rootpath);
    hdr{isubject} = ft_read_header(dataset_task{isubject,1});
    
    cfg = [];
    cfg.avgovertime = 'yes';
    cfg.latency = [max_time(isubject) max_time(isubject)];
    data_MEG_stim_maxtime{isubject} = ft_selectdata(cfg,data_MEG_stim_avg{isubject});
    
    % dipole fit
    cfg = [];
    cfg.numdipoles          = 1;
    cfg.grid                = subject_grid{isubject};    %   cfg.grid.pos    = template_grid.pos;
    cfg.vol                 = headmodel{isubject};
    cfg.grad                = hdr{isubject}.grad;
    cfg.channel             = 'MEG*1';
    %     cfg.channel     = {'MEG*2' 'MEG*3'};
    
    cfg.nonlinear           = 'no';
    cfg.reducerank          = 2; % to get out single orientation
    source_mag{isubject}    = ft_dipolefitting(cfg, data_MEG_stim_maxtime{isubject});
    
end

save('dipoles','source_mag*','-v7.3');

load('d:\WANDER\data\dipoles','source_mag*');


%% match positions with template

template_mri            = ft_read_mri('D:\WANDER\scripts\colin27_t1_tal_lin.nii');
% template_mri            = ft_read_mri('Z:\WANDER\fieldtrip21022019\template\anatomy\single_subj_T1_1mm.nii');
% template_mri = mri_realigned{isubject}
for isubject = slist
    
    pos_subject{isubject} = source_mag{isubject}.dip.pos;
    ori_subject{isubject} = source_mag{isubject}.dip.mom;
    
    index = find(subject_grid{isubject}.pos(:,1) == source_mag{isubject}.dip.pos(:,1) & subject_grid{isubject}.pos(:,2) == source_mag{isubject}.dip.pos(:,2) & subject_grid{isubject}.pos(:,3) == source_mag{isubject}.dip.pos(:,3));
    pos_template{isubject} = template_grid{1}.pos(index,:);
    ori_template{isubject} = source_mag{isubject}.dip.mom;
    
end

%% plot positions of dipoles

i = 1;
fig = figure;
for isubject = slist
    subplot(5,5,i); hold;
    ft_plot_slice(template_mri.anatomy, 'transform', template_mri.transform, 'location', pos_template{isubject}, 'plotmarker' ,pos_template{isubject},'markersize',30,'markercolor',[1 1 1])
% ft_plot_slice(template_mri.anatomy, 'location', pos_template{isubject}, 'plotmarker' ,pos_template{isubject},'markersize',30,'markercolor',[1 1 1])

    view(0,90);
    i = i + 1;
    title(num2str(isubject));
end

%% Export to BrainNet

% https://www.nitrc.org/docman/view.php/504/1190/BrainNet%20Viewer%20Manual%201.41.pdf
% The node file is defined as an ASCII text file with the suffix ‘node’. In the node file, there
% are 6 columns: columns 1-3 represent node coordinates, column 4 represents node
% colors, column 5 represents node sizes, and the last column represents node labels.
% Please note, a symbol ‘-‘(no ‘’) in column 6 means no labels. The user may put the
% modular information of the nodes into column 4, like ‘1, 2, 3…’ or other information to
% be shown by color. Column 5 could be set as nodal degree, centrality, T-value, etc. to
% emphasize nodal differences by size. You can generate your nodal file according to the
% requirements.

i = 1;
for isubject = slist
    d(i,:) = [pos_template{isubject},1,0.1,isubject];
    i = i + 1;
end

dlmwrite('D:\analysis\WANDER\dipoles.node',d,'\t') % BrainNet

%% extract ERF timecourses from full time period (ERF)

for isubject = slist
    
    fprintf('Loading subject %d \n',isubject);
    
    % cuelocked
    temp = load(['D:\WANDER\data\ERF\s' num2str(isubject) '_cue.mat'],'ERF_cue','ERF_cue_high','ERF_cue_low');
    ERF_cue       = temp.ERF_cue;
    ERF_cue_high  = temp.ERF_cue_high;
    ERF_cue_low   = temp.ERF_cue_low;
    clear temp
    
    % probelocked
    temp = load(['D:\WANDER\data\ERF\s' num2str(isubject) '_cue.mat'],'ERF_probe','ERF_probe_high','ERF_probe_low');
    ERF_probe       = temp.ERF_probe;
    ERF_probe_high  = temp.ERF_probe_high;
    ERF_probe_low   = temp.ERF_probe_low;
    
    % filter
    cfg = [];
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [1 30];
    ERF_cue         = ft_preprocessing(cfg,ERF_cue);
    ERF_cue_high    = ft_preprocessing(cfg,ERF_cue_high);
    ERF_cue_low   	= ft_preprocessing(cfg,ERF_cue_low);
    ERF_probe       = ft_preprocessing(cfg,ERF_probe);
    ERF_probe_high  = ft_preprocessing(cfg,ERF_probe_high);
    ERF_probe_low   = ft_preprocessing(cfg,ERF_probe_low);
    
    % only magnetometers, since dipole is based on those
    cfg = [];
    cfg.channel     = 'MEG*1';
    ERF_cue         = ft_selectdata(cfg,ERF_cue);
    ERF_cue_low     = ft_selectdata(cfg,ERF_cue_low);
    ERF_cue_high    = ft_selectdata(cfg,ERF_cue_high);
    ERF_probe       = ft_selectdata(cfg,ERF_probe);
    ERF_probe_low   = ft_selectdata(cfg,ERF_probe_low);
    ERF_probe_high  = ft_selectdata(cfg,ERF_probe_high);
    
    cfg = [];
    cfg.latency     = 'all';
    cfg.numdipoles  = 1;
    cfg.symmetry    = [];
    cfg.nonlinear   = 'no';  % use a fixed position
    cfg.gridsearch  = 'no';
    cfg.dip.pos     = pos_subject{isubject};
    %     cfg.dip.mom = source_mag{isubject}.dip.mom;
    cfg.vol         = headmodel{isubject};
    cfg.channel     = {'MEG*1'};
    cfg.senstype    = 'meg';
    cfg.grad        = hdr{isubject}.grad;
    
    source_cue_all{isubject}    = ft_dipolefitting(cfg, ERF_cue);
    source_cue_high{isubject}   = ft_dipolefitting(cfg, ERF_cue_high);
    source_cue_low{isubject}    = ft_dipolefitting(cfg, ERF_cue_low);
    source_probe_all{isubject}  = ft_dipolefitting(cfg, ERF_probe);
    source_probe_high{isubject} = ft_dipolefitting(cfg, ERF_probe_high);
    source_probe_low{isubject}  = ft_dipolefitting(cfg, ERF_probe_low);
    
%     clear memory
    source_cue_all{isubject}    = rmfield(source_cue_all{isubject},'cfg');
    source_cue_high{isubject}   = rmfield(source_cue_high{isubject},'cfg');
    source_cue_low{isubject}    = rmfield(source_cue_low{isubject},'cfg');
    source_probe_all{isubject}  = rmfield(source_probe_all{isubject},'cfg');
    source_probe_high{isubject} = rmfield(source_probe_high{isubject},'cfg');
    source_probe_low{isubject}  = rmfield(source_probe_low{isubject},'cfg');
    
    
    
    % clear ERF data for each subject
    clear ERF_cue* ERF_probe*
    
    %     % project to main component
    %
    %     source_cue_all{isubject}.dip.svd    = svdfft(source_cue_all{isubject}.dip.mom,1);
    %     source_cue_high{isubject}.dip.svd   = svdfft(source_cue_high{isubject}.dip.mom,1);
    %     source_cue_low{isubject}.dip.svd    = svdfft(source_cue_low{isubject}.dip.mom,1);
    %     source_probe_all{isubject}.dip.svd    = svdfft(source_probe_all{isubject}.dip.mom(:,end-20000+1:end),1);
    %     source_probe_high{isubject}.dip.svd   = svdfft(source_probe_high{isubject}.dip.mom(:,end-20000+1:end),1);
    %     source_probe_low{isubject}.dip.svd    = svdfft(source_probe_low{isubject}.dip.mom(:,end-20000+1:end),1);
    %
    
    % project to first component
    source_cue_all{isubject}.dip.svd    = source_cue_all{isubject}.dip.mom(1,:);
    source_cue_high{isubject}.dip.svd   = source_cue_high{isubject}.dip.mom(1,:);
    source_cue_low{isubject}.dip.svd    = source_cue_low{isubject}.dip.mom(1,:);
    source_probe_all{isubject}.dip.svd  = source_probe_all{isubject}.dip.mom(1,end-20000+1:end);
    source_probe_high{isubject}.dip.svd = source_probe_high{isubject}.dip.mom(1,end-20000+1:end);
    source_probe_low{isubject}.dip.svd  = source_probe_low{isubject}.dip.mom(1,end-20000+1:end);
    
end

save('D:\WANDER\dipole_sources','source_cue*','source_probe*','-v7.3');





load('/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/WANDER/dipole_sources','source_cue*','source_probe*');



%% put in data matrix
i = 1;
tc_cue          = nan(length(slist),40000);
tc_cue_high     = nan(length(slist),40000);
tc_cue_low      = nan(length(slist),40000);
tc_probe        = nan(length(slist),20000);
tc_probe_high   = nan(length(slist),20000);
tc_probe_low    = nan(length(slist),20000);

for isubject = slist
    tc_cue(i,1:length(source_cue_all{isubject}.dip.svd))            = source_cue_all{isubject}.dip.svd;
    tc_cue_high(i,1:length(source_cue_high{isubject}.dip.svd))      = source_cue_high{isubject}.dip.svd;
    tc_cue_low(i,1:length(source_cue_low{isubject}.dip.svd))        = source_cue_low{isubject}.dip.svd;
    tc_probe(i,1:length(source_probe_all{isubject}.dip.svd))        = source_probe_all{isubject}.dip.svd;
    tc_probe_high(i,1:length(source_probe_high{isubject}.dip.svd))  = source_probe_high{isubject}.dip.svd;
    tc_probe_low(i,1:length(source_probe_low{isubject}.dip.svd))    = source_probe_low{isubject}.dip.svd;
    
    tc_cue(i,1:length(source_cue_all{isubject}.dip.svd))            = source_cue_all{isubject}.dip.mom(1,:);
    tc_cue_high(i,1:length(source_cue_high{isubject}.dip.svd))      = source_cue_high{isubject}.dip.mom(1,:);
    tc_cue_low(i,1:length(source_cue_low{isubject}.dip.svd))        = source_cue_low{isubject}.dip.mom(1,:);
    tc_probe(i,1:length(source_probe_all{isubject}.dip.svd))        = source_probe_all{isubject}.dip.mom(1,end-20000+1:end);
    tc_probe_high(i,1:length(source_probe_high{isubject}.dip.svd))  = source_probe_high{isubject}.dip.mom(1,end-20000+1:end);
    tc_probe_low(i,1:length(source_probe_low{isubject}.dip.svd))    = source_probe_low{isubject}.dip.mom(1,end-20000+1:end);
    
    i = i + 1;
end

%% orient according to correlation with mean

time_cue        = linspace(-1.5,40-1.5,40000) - 0.033;
time_probe      = linspace(-20+1,1,20000) - 0.033 - 1/16; % t=0 at first nonstimulus

latency_cue     = 1600:1750; % define in time
latency_cue     = 1534:1834; % define in time
latency_cue     = 1534:1634; % define in time
latency_probe   = 19000:20000; % define in time

temp_cue_avg    = mean(tc_cue,1);
temp_probe_avg  = mean(tc_probe,1);
flip_cue        = ones(1,size(tc_cue,1));
flip_probe      = ones(1,size(tc_cue,1));

for iteration = 1 : 10
    disp(['orienting iteration',num2str(iteration)]);
    for isubject = 1 : size(tc_cue,1)
        if corr(tc_cue(isubject,latency_cue)',temp_cue_avg(latency_cue)','rows','complete') < 0
            flip_cue(isubject) = -flip_cue(isubject);
            tc_cue(isubject,:) = -tc_cue(isubject,:);
            tc_cue_high(isubject,:) = -tc_cue_high(isubject,:);
            tc_cue_low(isubject,:) = -tc_cue_low(isubject,:);
        end
        if corr(tc_probe(isubject,latency_probe)',temp_probe_avg(latency_probe)','rows','complete') < 0
            flip_probe(isubject) = -flip_probe(isubject);
            tc_probe(isubject,:) = -tc_probe(isubject,:);
            tc_probe_high(isubject,:) = -tc_probe_high(isubject,:);
            tc_probe_low(isubject,:) = -tc_probe_low(isubject,:);
        end
    end
    temp_cue_avg = mean(tc_cue,1);
    temp_probe_avg = mean(tc_probe,1);
end

figure; imagesc(tc_cue(:,latency_cue))

tc_cue_avg      = nanmean(tc_cue);
tc_cue_avg_high = nanmean(tc_cue_high);
tc_cue_avg_low  = nanmean(tc_cue_low);
sd_cue_avg      = std(tc_cue);
sd_cue_avg_high = std(tc_cue_high);
sd_cue_avg_low  = std(tc_cue_low);

%% plot cue-locked ERF
latency_cue         = [find(time_cue   > -0.2,1,'first') : find(time_cue   > 0.8,1,'first')];
latency_probe       = [find(time_probe > -0.4,1,'first') : find(time_probe > 0.6,1,'first')];

tc_probe_avg        = nanmean(tc_probe);
tc_probe_avg_high   = nanmean(tc_probe_high);
tc_probe_avg_low    = nanmean(tc_probe_low);
sd_probe_avg        = std(tc_probe);
sd_probe_avg_high   = std(tc_probe_high);
sd_probe_avg_low    = std(tc_probe_low);

figure;
subplot(1,2,1); hold;
patch([time_cue(latency_cue), time_cue(latency_cue(end:-1:1))],[tc_cue_avg_low(latency_cue) + sd_cue_avg_low(latency_cue)./sqrt(22), tc_cue_avg_low(latency_cue(end:-1:1)) - sd_cue_avg_low(latency_cue(end:-1:1))./sqrt(22)],[0 0 1],'FaceAlpha',.2,'LineStyle','None');
patch([time_cue(latency_cue), time_cue(latency_cue(end:-1:1))],[tc_cue_avg_high(latency_cue) + sd_cue_avg_high(latency_cue)./sqrt(22), tc_cue_avg_high(latency_cue(end:-1:1)) - sd_cue_avg_high(latency_cue(end:-1:1))./sqrt(22)],[1 0 0],'FaceAlpha',.2,'LineStyle','None');
plot(time_cue(latency_cue),tc_cue_avg_low(latency_cue),'b');
plot(time_cue(latency_cue),tc_cue_avg_high(latency_cue),'r');
axis([time_cue(latency_cue(1)),time_cue(latency_cue(end)),-10e-3,10e-3]);
ax = axis;
for t = 0 : 0.0625 : time_cue(latency_cue(end))
    line([t,t],[ax(3),ax(3)-ax(3)/5]);
end

% find peaks
i1 = find(time_cue > 0.030,1,'first');
i2 = find(time_cue > 0.060,1,'first');
[~, peak_i0] = max(tc_cue_avg(i1:i2));
peak_i0 = peak_i0+i1-1;
peak_t0 = time_cue(peak_i0);

i1 = find(time_cue > 0.070,1,'first');
i2 = find(time_cue > 0.100,1,'first');
[~, peak_i1] = min(tc_cue_avg(i1:i2));
peak_i1 = peak_i1+i1-1;
peak_t1 = time_cue(peak_i1);

line([peak_t0,peak_t0],[ax(3),ax(4)],'color','k');
line([peak_t1,peak_t1],[ax(3),ax(4)],'color','k');

[h0,p0,ci0,stats0] = ttest(tc_cue_high(:,peak_i0),tc_cue_low(:,peak_i0),'alpha',0.025)
[h1,p1,ci1,stats1] = ttest(tc_cue_high(:,peak_i1),tc_cue_low(:,peak_i1),'alpha',0.025)

% probe
subplot(1,2,2); hold;
patch([time_probe(latency_probe), time_probe(latency_probe(end:-1:1))],[tc_probe_avg_low(latency_probe) + sd_probe_avg_low(latency_probe)./sqrt(22), tc_probe_avg_low(latency_probe(end:-1:1)) - sd_probe_avg_low(latency_probe(end:-1:1))./sqrt(22)],[0 0 1],'FaceAlpha',.2,'LineStyle','None');
patch([time_probe(latency_probe), time_probe(latency_probe(end:-1:1))],[tc_probe_avg_high(latency_probe) + sd_probe_avg_high(latency_probe)./sqrt(22), tc_probe_avg_high(latency_probe(end:-1:1)) - sd_probe_avg_high(latency_probe(end:-1:1))./sqrt(22)],[1 0 0],'FaceAlpha',.2,'LineStyle','None');
plot(time_probe(latency_probe),tc_probe_avg_low(latency_probe),'b');
plot(time_probe(latency_probe),tc_probe_avg_high(latency_probe),'r');
axis([time_probe(latency_probe(1)),time_probe(latency_probe(end)),-10e-3,10e-3]);
ax = axis;
for t = -0.0625 : -0.0625 : time_probe(latency_probe(1))
    line([t,t],[ax(3),ax(3)-ax(3)/5]);
end

%% find peaks
i1 = find(time_probe > 0.040,1,'first');
i2 = find(time_probe > 0.060,1,'first');
[~, peak_i0] = max(tc_probe_avg(i1:i2));
peak_i0 = peak_i0+i1-1;
peak_t0 = time_probe(peak_i0);

i1 = find(time_probe > 0.1,1,'first');
i2 = find(time_probe > 0.180,1,'first');
[~, peak_i1] = min(tc_probe_avg(i1:i2));
peak_i1 = peak_i1+i1-1;
peak_t1 = time_probe(peak_i1);

i1 = find(time_probe > 0.15,1,'first');
i2 = find(time_probe > 0.25,1,'first');
[~, peak_i2] = max(tc_probe_avg(i1:i2));
peak_i2 = peak_i2+i1-1;
peak_t2 = time_probe(peak_i2);

line([peak_t0,peak_t0],[ax(3),ax(4)],'color','k');
line([peak_t1,peak_t1],[ax(3),ax(4)],'color','k');
line([peak_t2,peak_t2],[ax(3),ax(4)],'color','k');

[h0,p0,ci0,stats0] = ttest(tc_probe_high(:,peak_i0),tc_probe_low(:,peak_i0),'alpha',0.025)
[h1,p1,ci1,stats1] = ttest(tc_probe_high(:,peak_i1),tc_probe_low(:,peak_i1),'alpha',0.025)
[h2,p2,ci2,stats2] = ttest(tc_probe_high(:,peak_i2),tc_probe_low(:,peak_i2),'alpha',0.025)

[h,p,ci,stats] = ttest(tc_probe_high,tc_probe_low,'alpha',0.025)


% [h0,p0,ci0,stats0] = ttest(mean(tc_probe_high(:,firstpeak_i0-10:firstpeak_i0+20),2),mean(tc_probe_low(:,firstpeak_i0-10:firstpeak_i0+10),2),'alpha',0.025)
% [h1,p1,ci1,stats1] = ttest(mean(tc_probe_high(:,firstpeak_i1-10:firstpeak_i1+20),2),mean(tc_probe_low(:,firstpeak_i1-10:firstpeak_i1+10),2),'alpha',0.025)
% [h2,p2,ci2,stats2] = ttest(mean(tc_probe_high(:,firstpeak_i2-10:firstpeak_i2+20),2),mean(tc_probe_low(:,firstpeak_i2-10:firstpeak_i2+10),2),'alpha',0.025)
%
savefig('d:\WANDER\images\article\ERF_dipole');
print('-dsvg','w:\WANDER\images\article\ERF_dipole.svg');

%% FFT

latency_cue     = [find(time_cue > 0, 1, 'first') : find(time_cue < 10, 1, 'last')];
latency_probe   = [find(time_probe > -10, 1, 'first') : find(time_probe > 0, 1, 'first')];

Fs = 1000;
L_cue = length(latency_cue);
f_cue = Fs*(0:(L_cue/2))/L_cue;

L_probe = length(latency_probe);
f_probe = Fs*(0:(L_probe/2))/L_probe;

FFT_cue_avg                 = fft(tc_cue(:,latency_cue)');
FFT_cue_avg_high            = fft(tc_cue_high(:,latency_cue)');
FFT_cue_avg_low             = fft(tc_cue_low(:,latency_cue)');
FFT_cue_avg                 = abs(FFT_cue_avg/L_cue);
FFT_cue_avg                 = FFT_cue_avg(1:L_cue/2+1,:);
FFT_cue_avg(2:end-1,:)      = 2*FFT_cue_avg(2:end-1,:);
FFT_cue_avg_high            = abs(FFT_cue_avg_high/L_cue);
FFT_cue_avg_high            = FFT_cue_avg_high(1:L_cue/2+1,:);
FFT_cue_avg_high(2:end-1,:) = 2*FFT_cue_avg_high(2:end-1,:);
FFT_cue_avg_low             = abs(FFT_cue_avg_low/L_cue);
FFT_cue_avg_low             = FFT_cue_avg_low(1:L_cue/2+1,:);
FFT_cue_avg_low(2:end-1,:)  = 2*FFT_cue_avg_low(2:end-1,:);

FFT_probe_avg                 = fft(tc_probe(:,latency_probe)');
FFT_probe_avg_high            = fft(tc_probe_high(:,latency_probe)');
FFT_probe_avg_low             = fft(tc_probe_low(:,latency_probe)');
FFT_probe_avg                 = abs(FFT_probe_avg/L_probe);
FFT_probe_avg                 = FFT_probe_avg(1:L_probe/2+1,:);
FFT_probe_avg(2:end-1,:)      = 2*FFT_probe_avg(2:end-1,:);
FFT_probe_avg_high            = abs(FFT_probe_avg_high/L_probe);
FFT_probe_avg_high            = FFT_probe_avg_high(1:L_probe/2+1,:);
FFT_probe_avg_high(2:end-1,:) = 2*FFT_probe_avg_high(2:end-1,:);
FFT_probe_avg_low             = abs(FFT_probe_avg_low/L_probe);
FFT_probe_avg_low             = FFT_probe_avg_low(1:L_probe/2+1,:);
FFT_probe_avg_low(2:end-1,:)  = 2*FFT_probe_avg_low(2:end-1,:);

[~, ilow_cue] = find(f_cue>12,1,'first');
[~, ihigh_cue] = find(f_cue>20,1,'first');

[~, ilow_probe] = find(f_probe>12,1,'first');
[~, ihigh_probe] = find(f_probe>20,1,'first');


%% FOR REVIEWER QUESTION ABOUT SUBHARMONICS
[~, ilow_probe] = find(f_probe>0,1,'first');
[~, ihigh_probe] = find(f_probe>20,1,'first');


% plot CUE FFT
figure;
subplot(2,1,1);
plot(f_cue(ilow_cue:ihigh_cue),FFT_cue_avg(ilow_cue:ihigh_cue,:));
axis tight
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

subplot(2,1,2);
plot(f_cue(ilow_cue:ihigh_cue),mean(FFT_cue_avg(ilow_cue:ihigh_cue,:),2))
title('Single-Sided Amplitude Spectrum of X(t)');
axis tight
xlabel('f (Hz)')
ylabel('|P1(f)|')

% plot CUE absolute differences
figure; subplot(3,1,1);
hold;
plot(f_cue(ilow_cue:ihigh_cue),mean(FFT_cue_avg_high(ilow_cue:ihigh_cue,:),2),'r');
plot(f_cue(ilow_cue:ihigh_cue),mean(FFT_cue_avg_low(ilow_cue:ihigh_cue,:),2),'b');
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
axis tight

[h,p,ci,stats] = ttest(FFT_cue_avg_high(ilow_cue:ihigh_cue,:)',FFT_cue_avg_low(ilow_cue:ihigh_cue,:)','alpha',0.025);
subplot(3,1,2);
plot(f_cue(ilow_cue:ihigh_cue),p); axis tight

subplot(3,1,3);
plot(f_cue(ilow_cue:ihigh_cue),h); axis tight

% relative CUE difference
FFT_cue_avg_diff = (FFT_cue_avg_low - FFT_cue_avg_high) ./ (FFT_cue_avg_low + FFT_cue_avg_high);

figure; hold;
plot(f_cue(ilow_cue:ihigh_cue),FFT_cue_avg_diff(ilow_cue:ihigh_cue,:));
title('CUE-LOCKED Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure; hold;
plot(f_cue(ilow_cue:ihigh_cue),mean(FFT_cue_avg_diff(ilow_cue:ihigh_cue,:),2));
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% plot PROBE FFT
figure;
plot(f_probe(ilow_probe:ihigh_probe),FFT_probe_avg(ilow_probe:ihigh_probe,:))
title('PROBE-LOCKED: Single-Sided Amplitude Spectrum of X(t)');
axis tight
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure;
plot(f_probe(ilow_probe:ihigh_probe),mean(FFT_probe_avg(ilow_probe:ihigh_probe,:),2))
title('Single-Sided Amplitude Spectrum of X(t)');
axis tight
xlabel('f (Hz)')
ylabel('|P1(f)|')

% plot PROBE absolute differences
figure;
subplot(2,1,1); hold;
plot(f_probe(ilow_probe:ihigh_probe),mean(FFT_probe_avg_high(ilow_probe:ihigh_probe,:),2),'r');
plot(f_probe(ilow_probe:ihigh_probe),mean(FFT_probe_avg_low(ilow_probe:ihigh_probe,:),2),'b');
% plot(f_probe(ilow_probe:ihigh_probe),mean(FFT_probe_avg_high(ilow_probe:ihigh_probe,:)' ./ (FFT_probe_avg_high(ilow_probe:ihigh_probe,:)' + FFT_probe_avg_low(ilow_probe:ihigh_probe,:)'),1),'r');
% plot(f_probe(ilow_probe:ihigh_probe),mean(FFT_probe_avg_low(ilow_probe:ihigh_probe,:)' ./ (FFT_probe_avg_high(ilow_probe:ihigh_probe,:)' + FFT_probe_avg_low(ilow_probe:ihigh_probe,:)'),1),'b');
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
axis tight
% print('-dsvg','w:\WANDER\images\article\FFT_16Hz_dipole.svg');



%  plot RESPONSE FOR REVIEWER
fig = figure;
% subplot(1,2,1); hold;
plot(f_probe(ilow_probe:ihigh_probe),mean(FFT_probe_avg(ilow_probe:ihigh_probe,:),2),'k');
xlabel('f (Hz)')
ylabel('|P1(f)|')
axis tight
t = xticks;
xticks([0, 8, 10, 16, 20]);
xlim([0, 20]);
box off
exportgraphics(fig, '/network/lustre/iss01/charpier/analyses/stephen.whitmarsh/WANDER/revision_EJN/FFT_16Hz_dipole_low_freq.jpg');

[h,p,ci,stats] = ttest(FFT_probe_avg_high(ilow_probe:ihigh_probe,:)',FFT_probe_avg_low(ilow_probe:ihigh_probe,:)','alpha',0.025);
[h,p,ci,stats] = ttest(FFT_probe_avg_high(ilow_probe:ihigh_probe,:)' ./ (FFT_probe_avg_high(ilow_probe:ihigh_probe,:)' + FFT_probe_avg_low(ilow_probe:ihigh_probe,:)'),FFT_probe_avg_low(ilow_probe:ihigh_probe,:)' ./ (FFT_probe_avg_high(ilow_probe:ihigh_probe,:)' + FFT_probe_avg_low(ilow_probe:ihigh_probe,:)'),'alpha',0.025);


subplot(2,1,2);
plot(f_probe(ilow_probe:ihigh_probe),p); axis tight
q = find(f_probe(ilow_probe:ihigh_probe) == 16);
stats.tstat(q)
stats.df(q)
p(q)

x = f_probe(ilow_probe:ihigh_probe);
y = p;
plot(x, y)

[x; y; stats.tstat]






















subplot(3,1,3);
plot(f_probe(ilow_probe:ihigh_probe),h); axis tight





% relative PROBE difference
FFT_cue_avg_diff = (FFT_cue_avg_low - FFT_cue_avg_high) ./ (FFT_cue_avg_low + FFT_cue_avg_high);

figure; hold;
plot(f_probe(ilow_probe:ihigh_probe),FFT_cue_avg_diff(ilow_probe:ihigh_probe,:));
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure; hold;
plot(f_probe(ilow_probe:ihigh_probe),mean(FFT_cue_avg_diff(ilow_probe:ihigh_probe,:),2));
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


%% FFT welch method

% put in FT format
for isubject = slist
    data_cue{isubject}                  = [];
    data_cue{isubject}.avg              = source_cue_all{isubject}.dip.svd;
    data_cue{isubject}.time             = source_cue_all{isubject}.time;
    data_cue{isubject}.label            = {'dipole'};
    
    data_cue_high{isubject}             = [];
    data_cue_high{isubject}.avg         = source_cue_high{isubject}.dip.svd;
    data_cue_high{isubject}.time        = source_cue_high{isubject}.time;
    data_cue_high{isubject}.label       = {'dipole'};
    
    data_cue_low{isubject}              = [];
    data_cue_low{isubject}.avg          = source_cue_low{isubject}.dip.svd;
    data_cue_low{isubject}.time         = source_cue_low{isubject}.time;
    data_cue_low{isubject}.label        = {'dipole'};
    
    data_probe{isubject}                = [];
    data_probe{isubject}.avg            = source_probe_all{isubject}.dip.svd;
    data_probe{isubject}.time           = source_probe_all{isubject}.time(end-20000+1:end);
    data_probe{isubject}.label          = {'dipole'};
    
    data_probe_high{isubject}           = [];
    data_probe_high{isubject}.avg       = source_probe_high{isubject}.dip.svd;
    data_probe_high{isubject}.time      = source_probe_high{isubject}.time(end-20000+1:end);
    data_probe_high{isubject}.label     = {'dipole'};
    
    data_probe_low{isubject}            = [];
    data_probe_low{isubject}.avg        = source_probe_low{isubject}.dip.svd;
    data_probe_low{isubject}.time       = source_probe_low{isubject}.time(end-20000+1:end);
    data_probe_low{isubject}.label      = {'dipole'};
    
    cfg                                 = [];
    cfg.output                          = 'pow';
    cfg.method                          = 'mtmconvol';
    cfg.keeptrials                      = 'no';
    cfg.taper                           = 'hanning';
    cfg.foi                             = [1:30];
    cfg.t_ftimwin                       = ones(size(cfg.foi))*2;
    cfg.toi                             = [-1 : 0.01 : 10];
    TFR_cue{isubject}                   = ft_freqanalysis(cfg, data_cue{isubject});
    TFR_cue_high{isubject}              = ft_freqanalysis(cfg, data_cue_high{isubject});
    TFR_cue_low{isubject}               = ft_freqanalysis(cfg, data_cue_low{isubject});
    TFR_cue_diff{isubject}              = TFR_cue_high{isubject};
    TFR_cue_diff{isubject}.powspctrm    = (TFR_cue_low{isubject}.powspctrm - TFR_cue_high{isubject}.powspctrm) ./ (TFR_cue_low{isubject}.powspctrm + TFR_cue_high{isubject}.powspctrm);
    
    cfg.toi                             = [-10 : 0.01 : 1];
    
    TFR_probe{isubject}                 = ft_freqanalysis(cfg, data_probe{isubject});
    TFR_probe_high{isubject}            = ft_freqanalysis(cfg, data_probe_high{isubject});
    TFR_probe_low{isubject}             = ft_freqanalysis(cfg, data_probe_low{isubject});
    TFR_probe_diff{isubject}            = TFR_probe_high{isubject};
    TFR_probe_diff{isubject}.powspctrm  = (TFR_probe_low{isubject}.powspctrm - TFR_probe_high{isubject}.powspctrm) ./ (TFR_probe_low{isubject}.powspctrm + TFR_probe_high{isubject}.powspctrm);
    
end

TFR_cue_avg         = ft_freqgrandaverage([],TFR_cue{slist});
TFR_cue_diff_avg    = ft_freqgrandaverage([],TFR_cue_diff{slist});
TFR_probe_avg       = ft_freqgrandaverage([],TFR_probe{slist});
TFR_probe_diff_avg  = ft_freqgrandaverage([],TFR_probe_diff{slist});

figure;
cfg = [];
cfg.colormap = parula(1000);
ft_singleplotTFR(cfg,TFR_cue_avg);

figure;
cfg = [];
ft_singleplotTFR(cfg,TFR_cue_diff_avg);

figure;
plot(TFR_cue_diff_avg.time, squeeze(TFR_cue_diff_avg.powspctrm(1,16,:)))

figure;
cfg = [];
cfg.colormap = parula(1000);
ft_singleplotTFR(cfg,TFR_probe_avg);

figure;
cfg = [];
ft_singleplotTFR(cfg,TFR_probe_diff_avg);

figure;
plot(TFR_probe_diff_avg.time, squeeze(TFR_probe_diff_avg.powspctrm(1,16,:)))

%% get SD for plotting
i = 1;
for isubject = slist
    tc_TFR_cue(i,:)         = TFR_cue{isubject}.powspctrm(1,16,:);
    tc_TFR_cue_high(i,:)    = TFR_cue_high{isubject}.powspctrm(1,16,:);
    tc_TFR_cue_low(i,:)     = TFR_cue_low{isubject}.powspctrm(1,16,:);
    tc_TFR_cue_diff(i,:)    = TFR_cue_diff{isubject}.powspctrm(1,16,:);
    
    tc_TFR_probe(i,:)       = TFR_probe{isubject}.powspctrm(1,16,:);
    tc_TFR_probe_high(i,:)  = TFR_probe_high{isubject}.powspctrm(1,16,:);
    tc_TFR_probe_low(i,:)   = TFR_probe_low{isubject}.powspctrm(1,16,:);
    tc_TFR_probe_diff(i,:)  = TFR_probe_diff{isubject}.powspctrm(1,16,:);
    
    i = i + 1;
end

m_cue       = mean(tc_TFR_cue,1);
sd_cue      = std(tc_TFR_cue,1)/sqrt(22);
m_cue_diff  = mean(tc_TFR_cue_diff,1);
sd_cue_diff = std(tc_TFR_cue_diff,1)/sqrt(22);
m_probe_diff  = mean(tc_TFR_probe_diff,1);
sd_probe_diff = std(tc_TFR_probe_diff,1)/sqrt(22);

figure; hold;
plot(TFR_cue_diff_avg.time,m_cue_diff);
plot(TFR_cue_diff_avg.time,m_cue_diff+sd_cue_diff);
plot(TFR_cue_diff_avg.time,m_cue_diff-sd_cue_diff);

figure; hold;
plot(TFR_probe_diff_avg.time,m_probe_diff);
plot(TFR_probe_diff_avg.time,m_probe_diff+sd_probe_diff);
plot(TFR_probe_diff_avg.time,m_probe_diff-sd_probe_diff);


%% STATISTICS

load('neuromag306mag_neighb_last');

cfg = [];
cfg.channel             = 'all';
cfg.frequency           = 16;
cfg.latency             = [-8 -1];
cfg.avgovertime         = 'yes';
cfg.neighbours          = neighbours;
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.numrandomization    = 1000;
cfg.alpha               = 0.025;
cfg.clusteralpha        = 0.05;
cfg.tail                = 0;
cfg.design(1,:)         = [1:length(slist) 1:length(slist)];
cfg.design(2,:)         = [ones(size(slist)) ones(size(slist))*2];
cfg.uvar                = 1;
cfg.ivar                = 2;
cfg.parameter           = 'powspctrm';
stat                    = ft_freqstatistics(cfg, TFR_probe_high{slist}, TFR_probe_low{slist});












latency = 1633:3000;

% tc_avg      = nanmean(tc(:,latency));
% tc_avg_high = nanmean(tc_high(:,latency));
% tc_avg_low  = nanmean(tc_low(:,latency));

Fs = 1000;
L = length(latency);
f = Fs*(0:(L/2))/L;

% FFT_avg               = fft(tc_avg);
% FFT_avg_high          = fft(tc_avg_high);
% FFT_avg_low           = fft(tc_avg_low);

FFT_avg               = fft(tc(:,latency)');
FFT_avg_high          = fft(tc_high(:,latency)');
FFT_avg_low           = fft(tc_low(:,latency)');

FFT_avg               = abs(FFT_avg/L);
FFT_avg               = FFT_avg(1:L/2+1,:);
FFT_avg(2:end-1,:)    = 2*FFT_avg(2:end-1,:);

FFT_avg_high          = abs(FFT_avg_high/L);
FFT_avg_high          = FFT_avg_high(1:L/2+1,:);
FFT_avg_high(2:end-1,:) = 2*FFT_avg_high(2:end-1,:);

FFT_avg_low           = abs(FFT_avg_low/L);
FFT_avg_low           = FFT_avg_low(1:L/2+1,:);
FFT_avg_low(2:end-1,:)  = 2*FFT_avg_low(2:end-1,:);

[~, ilow] = find(f>1,1,'first');
[~, ihigh] = find(f>30,1,'first');

figure;
plot(f(ilow:ihigh),FFT_avg(ilow:ihigh,:))
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure;
plot(f(ilow:ihigh),mean(FFT_avg(ilow:ihigh,:),2))
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')


figure; hold;
plot(f(ilow:ihigh),mean(FFT_avg_high(ilow:ihigh,:),2),'b');
plot(f(ilow:ihigh),mean(FFT_avg_low(ilow:ihigh,:),2),'r');
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

% relative difference
FFT_avg_diff = (FFT_avg_low - FFT_avg_high) ./ (FFT_avg_low + FFT_avg_high);

figure; hold;
plot(f(ilow:ihigh),FFT_avg_diff(ilow:ihigh,:));
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

figure; hold;
plot(f(ilow:ihigh),mean(FFT_avg_diff(ilow:ihigh,:),2));
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%% extract timecourses from full time period of each trial
for isubject = slist
    
    fprintf('Loading subject %d \n',isubject);
    data_MEG = WANDER_epoch_MEG(isubject,0,rootpath,0);
    artdef   = WANDER_artefact_detection_MEG(isubject,0,rootpath,0);
    
    for iblock = 1 : 4
        
        % select all channels, only correct rejections
        cfg                 = [];
        cfg.channel         = 'MEG';
        cfg.trials          = find(data_MEG{iblock}.trialinfo(:,3) == 4);
        data_MEG{iblock}    = ft_selectdata(cfg,data_MEG{iblock});
        data_MEG{iblock}    = rmfield(data_MEG{iblock},'cfg');
        
        cfg = [];
        cfg.reject                      = 'nan'; % make bugreport, does not work with cfg.artfctdef.reject
        cfg.artfctdef.minaccepttim      = 0;
        cfg.artfctdef                   = artdef{iblock};
        data_MEG{iblock}                 = ft_rejectartifact(cfg,data_MEG{iblock});
        
        for itrial = 1 : size(data_MEG{iblock}.trial,2)
            fprintf('Cutting off beginning of trial %d block %d \n',itrial,iblock);
            data_MEG{iblock}.trial{itrial} = data_MEG{iblock}.trial{itrial}(:,1501:end);
            data_MEG{iblock}.time{itrial}  = data_MEG{iblock}.time{itrial}(1501:end);
            fprintf('Shifting timecourse of trial %d block %d \n',itrial,iblock);
            data_MEG{iblock}.time{itrial}  = data_MEG{iblock}.time{itrial} - data_MEG{iblock}.time{itrial}(end) + 1.000;
        end
        
        % change timecourse according to latency of stimulation
        for itrial = size(data_MEG{iblock}.time,2)
            data_MEG{iblock}.time{itrial} = data_MEG{iblock}.time{itrial} + 0.033;
        end
    end
    
    % concatinate trials over blocks
    data_MEG_append = ft_appenddata([],data_MEG{:});
    
    clear data_MEG
    
    for itrial = 1 : size(data_MEG_append.trial,2)
        
        cfg = [];
        cfg.trials = itrial;
        sel = ft_selectdata(cfg,data_MEG_append);
        sel.avg = sel.trial{1};
        
        cfg = [];
        cfg.latency     = 'all';
        cfg.numdipoles  = 1;
        cfg.symmetry    = [];
        cfg.nonlinear   = 'no';  % use a fixed position
        cfg.gridsearch  = 'no';
        cfg.dip.pos     = pos_subject{isubject};
        %     cfg.dip.mom = source_mag{isubject}.dip.mom;
        cfg.vol         = headmodel{isubject};
        cfg.channel     = {'MEG*1'};
        cfg.senstype    = 'meg';
        cfg.grad        = hdr{isubject}.grad;
        
        temp            = ft_dipolefitting(cfg, sel);
        
        source{isubject}.label = {'dip'};
        source{isubject}.time{itrial} = temp.time;
        source{isubject}.dimord = temp.dimord;
        source{isubject}.dip{itrial} = temp.dip; % for archiving
        source{isubject}.trial{itrial} = temp.dip.mom(1,:);
    end
    trialinfo{isubject} = data_MEG_append.trialinfo;
    
end
clear data_MEG data_MEG_append

save('dipole_timecourse','source','trialinfo','-v7.3');

for isubject = slist
    cfg = [];
    cfg.method = 'mtmconvol';
    cfg.foi = 1:30;
    cfg.toi = -20:0.05:1;
    cfg.taper = 'hanning';
    cfg.t_ftimwin = ones(size(cfg.foi));
    cfg.keeptrials = 'yes';
    cfg.pad='nextpow2';
    cfg.output = 'fourier';
    source_TFR{isubject} = ft_freqanalysis(cfg,source{isubject});
    
end

for isubject = slist
    cfg           = [];
    cfg.method    = 'coh';
    coh{isubject}           = ft_connectivityanalysis(cfg, source_TFR{isubject});
end


%% EXTRACT RMS for SINGLE STIM
load('dipole_timecourse');

for isubject = slist
    
    fprintf('Working on subject %d \n',isubject);
    dat = nan(size(source{isubject}.trial,2),16*10,floor(1000/16)+1); % trial x stim x time
    
    for itrial = 1 : size(source{isubject}.trial,2)
        
        for istim = 1 : 16*10
            start_sample = find(source{isubject}.time{itrial} > 0-(1/16)*istim, 1, 'first');
            end_sample   = start_sample + floor(1000/16);
            
            dat(itrial,istim,:) = source{isubject}.trial{itrial}(start_sample:end_sample);
            disp(['subject: ', num2str(isubject), ', trial: ',num2str(itrial), ', sample: ', num2str(start_sample)]);
        end
    end
    
    r{isubject} = rms(nanmean(dat,2),3);
    %     figure; imagesc(squeeze(mean(dat,2)))
    
end

save('W:\WANDER\data\SSEF\rms','r');















%%
for isubject = slist
    
    for itrial = 1 : size(source_TFR{isubject}.powspctrm,1)
        corr_trial{isubject}(itrial,:) = corr(squeeze(source_TFR{isubject}.powspctrm(itrial,1,16,:)),squeeze(source_TFR{isubject}.powspctrm(itrial,1,:,:))','rows','complete','type','spearman');
    end
    
    F = ceil(2 * tiedrank(trialinfo{isubject}(:,2)) / length(trialinfo{isubject}(:,2)));
    trialinfo{isubject}(:,7) = ones(size(F));
    trialinfo{isubject}((F==2),7)  = 2;
    
    corr_trial_avg(isubject,:)      = nanmean(corr_trial{isubject});
    corr_trial_avg_high(isubject,:) = nanmean(corr_trial{isubject}(trialinfo{isubject}(:,7) == 2,:));
    corr_trial_avg_low(isubject,:)  = nanmean(corr_trial{isubject}(trialinfo{isubject}(:,7) == 1,:));
    corr_trial_avg_diff(isubject,:) = corr_trial_avg_high(isubject,:) - corr_trial_avg_low(isubject,:);
end

figure; imagesc(corr_trial_avg(slist,:)')
figure; imagesc(corr_trial_avg_high(slist,:)')
figure; imagesc(corr_trial_avg_low(slist,:)')
figure; imagesc(corr_trial_avg_diff(slist,:)')

corr_trial_avg_high

plot(mean(corr_trial_avg(slist,:)));
plot(corr_trial{isubject}(1:10,:)')


figure;
i = 1;
for isubject = slist
    subplot(5,5,i);
    % plot
    cfg         = [];
    cfg.zlim    = 'maxabs';
    cfg.layout  = 'neuromag306mag';
    cfg.channel = 'dip';
    cfg.baseline = [0 0];
    %     cfg.baselinetype = 'relative';
    ft_singleplotTFR(cfg,source_TFR{isubject});
    i = i + 1;
end
















time = latency/1000-1.5-0.033;

figure; hold;
plot(time,tc_avg);
axis tight;

figure; hold;
plot(time,tc_avg_high,'b');
plot(time,tc_avg_low,'r');
axis tight







s = svdfft(source_all{isubject}.dip.mom,1);
figure; hold;
plot(s(1:1000));
plot(source_all{isubject}.dip.mom(:,1:1000)')

data_MEG_stim_GA        = ft_timelockgrandaverage([],data_MEG_stim_avg{slist});
data_MEG_stim_high_GA   = ft_timelockgrandaverage([],data_MEG_stim_avg_high{slist});
data_MEG_stim_low_GA    = ft_timelockgrandaverage([],data_MEG_stim_avg_low{slist});
data_MEG_stim_diff_GA   = ft_timelockgrandaverage([],data_MEG_stim_diff{slist});
data_MEG_MI_GA          = ft_timelockgrandaverage([],data_MEG_stim_MI{slist});
data_MEG_MI_high_GA     = ft_timelockgrandaverage([],data_MEG_stim_MI_high{slist});
data_MEG_MI_low_GA      = ft_timelockgrandaverage([],data_MEG_stim_MI_low{slist});
data_MEG_MI_diff_GA     = ft_timelockgrandaverage([],data_MEG_MI_diff{slist});
data_MEG_rms_GA         = ft_timelockgrandaverage([],data_MEG_stim_rms{slist});
data_MEG_rms_high_GA    = ft_timelockgrandaverage([],data_MEG_stim_rms_high{slist});
data_MEG_rms_low_GA     = ft_timelockgrandaverage([],data_MEG_stim_rms_low{slist});
data_MEG_rms_diff_GA    = ft_timelockgrandaverage([],data_MEG_stim_rms_diff{slist});

% plot timecourses
cfg         = [];
cfg.zlim    = 'maxabs';
cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEG*1';
figure; ft_multiplotER(cfg,data_MEG_stim_avg{25});

% shows no clear difference, perhaps because these are averages of absolute
figure; ft_multiplotER(cfg,data_MEG_stim_GA);
figure; ft_multiplotER(cfg,data_MEG_stim_high_GA,data_MEG_stim_low_GA);
figure; ft_multiplotER(cfg,data_MEG_stim_diff_GA);

cfg.layout  = 'neuromag306cmb';
cfg.channel = 'MEG*3';

% shows no clear difference, perhaps because these are averages of absolute
figure; ft_multiplotER(cfg,data_MEG_stim_high_GA,data_MEG_stim_low_GA);

% plot topography
cfg = [];
cfg.zlim    = 'maxabs';
cfg.layout  = 'neuromag306mag';
cfg.channel = {'MEG*1'};

% topography
figure; ft_topoplotER(cfg,data_MEG_stim_GA);

cfg.xlim = [0 : 0.01 : 0.06];
figure; ft_topoplotER(cfg,data_MEG_stim_diff_GA);

figure; ft_topoplotER(cfg,data_MEG_MI_GA);
figure; ft_topoplotER(cfg,data_MEG_MI_diff_GA);
figure; ft_topoplotER(cfg,data_MEG_MI_high_GA);
figure; ft_topoplotER(cfg,data_MEG_MI_low_GA);

figure; ft_topoplotER(cfg,data_MEG_rms_GA);
figure; ft_topoplotER(cfg,data_MEG_rms_diff_GA);
figure; ft_topoplotER(cfg,data_MEG_rms_high_GA);
figure; ft_topoplotER(cfg,data_MEG_rms_low_GA);




cfg         = [];
cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEGMAG';
cfg.zlim    = 'absmax';
figure; ft_multiplotER(cfg,data_MEG_stim_bin{isubject}{:});

cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEGMAG';
cfg.zlim    = 'absmax';
cfg.colorbar = 'yes';
cfg.parameter = 'rms';
figure;
ft_multiplotER(cfg,data_MEG_stim_binhist{isubject});
ft_multiplotER(cfg,data_MEG_stim_binhist_low{isubject},data_MEG_stim_binhist_high{isubject});


cfg         = [];
cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEG*1';
cfg.zlim    = [-0.02 0.02];
figure; ft_topoplotER(cfg,data_MEG_stim_MI_high{isubject},data_MEG_stim_MI_low{isubject});
cfg.zlim    = 'absmax';
figure; ft_topoplotER(cfg,data_MEG_stim_MI_diff{isubject});


% STATISTICS

load('neuromag306mag_neighb_last');

cfg = [];
cfg.channel             = 'MEG*1';
cfg.neighbours          = neighbours;
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.correctm            = 'cluster';
cfg.numrandomization    = 1000;
cfg.alpha               = 0.05;
cfg.clusteralpha        = 0.05;
cfg.tail                = 0;
cfg.design(1,:)         = [1:length(slist) 1:length(slist)];
cfg.design(2,:)         = [ones(size(slist)) ones(size(slist))*2];
cfg.uvar                = 1;
cfg.ivar                = 2;
cfg.parameter           = 'avg';
stat_MI                 = ft_timelockstatistics(cfg, data_MEG_stim_MI_high{slist}, data_MEG_stim_MI_low{slist});
stat_MI                 = ft_timelockstatistics(cfg, data_MEG_stim_MI_high_norm{slist}, data_MEG_stim_MI_low_norm{slist});


stat_rms                = ft_timelockstatistics(cfg, data_MEG_stim_rms_low{slist}, data_MEG_stim_rms_high{slist});


cfg         = [];
cfg.layout  = 'neuromag306mag';
cfg.channel = 'MEG*1';
cfg.highlight = 'on';
cfg.parameter = 'stat';
cfg.highlightchannel = {stat_MI.label{stat_MI.mask}};
% cfg.zlim    = [-0.02 0.02];
figure; ft_topoplotER(cfg,stat_MI);

cfg.highlightchannel = {stat_rms.label{stat_rms.mask}};

cfg.zlim    = 'maxabs';

figure; ft_topoplotER(cfg,stat_rms);

stat_MI = ft_combineplanar([],stat_MI);
cfg         = [];
cfg.layout  = 'neuromag306grad';
cfg.channel = 'MEG*1';
cfg.highlight = 'on';
cfg.parameter = 'stat';
cfg.highlightchannel = {stat_MI.label{stat_MI.mask}};
% cfg.zlim    = [-0.02 0.02];
figure; ft_topoplotER(cfg,stat_MI);


