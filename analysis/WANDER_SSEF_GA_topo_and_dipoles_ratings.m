
load('dipoles.mat'); 
slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD


for isubject = slist
 
    
    % load headmodel
    [mri_segmented{isubject}, headmodel{isubject}, subject_grid{isubject}, template_grid, ~] = WANDER_grid(isubject,1,0);
    
    % load sensor information
    [dataset_task, ~] = WANDER_subjectinfo;
    hdr{isubject} = ft_read_header(dataset_task{isubject,1});
    
    % load timcourse data
    ERF_probe_ratings = WANDER_ERF_ratings(isubject,0,1);

    % filter
    for irating = 1 : 4
        cfg = [];
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [1 30];
        ERF_probe_ratings{irating} = ft_preprocessing(cfg,ERF_probe_ratings{irating} );
    end
    
    cfg = [];
    cfg.latency     = 'all';
    cfg.numdipoles  = 1;
    cfg.symmetry    = [];
    cfg.nonlinear   = 'no';  % use a fixed position
    cfg.gridsearch  = 'no';
    cfg.dip.pos     = source_mag{isubject}.dip.pos;
    cfg.vol         = headmodel{isubject};
    cfg.channel     = {'MEG*1'};
    cfg.senstype    = 'meg';
    cfg.grad = hdr{isubject}.grad;
    dip_1{isubject}  = ft_dipolefitting(cfg, ERF_probe_ratings{1});
    dip_2{isubject}  = ft_dipolefitting(cfg, ERF_probe_ratings{2});
    dip_3{isubject}  = ft_dipolefitting(cfg, ERF_probe_ratings{3});
    dip_4{isubject}  = ft_dipolefitting(cfg, ERF_probe_ratings{4});
    
    clear ERF_probe_ratings
    
    dip_1{isubject}.dip.svd  = dip_1{isubject}.dip.mom(:,end-11000+1:end);
    dip_2{isubject}.dip.svd  = dip_2{isubject}.dip.mom(:,end-11000+1:end);
    dip_3{isubject}.dip.svd  = dip_3{isubject}.dip.mom(:,end-11000+1:end);
    dip_4{isubject}.dip.svd  = dip_4{isubject}.dip.mom(:,end-11000+1:end);
end

for isubject = slist
    dip_1_tc{isubject} = dip_1{isubject};    
    dip_1_tc{isubject}.avg = dip_1_tc{isubject}.dip.mom;    
    dip_1_tc{isubject}.label = {'dip1','dip2','dip3'};
    dip_1_tc{isubject} = rmfield(dip_1_tc{isubject},{'dip','Vdata','Vmodel'});
    dip_2_tc{isubject} = dip_2{isubject};    
    dip_2_tc{isubject}.avg = dip_2_tc{isubject}.dip.mom;    
    dip_2_tc{isubject}.label = {'dip1','dip2','dip3'};
    dip_2_tc{isubject} = rmfield(dip_2_tc{isubject},{'dip','Vdata','Vmodel'});
    dip_3_tc{isubject} = dip_3{isubject};    
    dip_3_tc{isubject}.avg = dip_3_tc{isubject}.dip.mom;    
    dip_3_tc{isubject}.label = {'dip1','dip2','dip3'};
    dip_3_tc{isubject} = rmfield(dip_3_tc{isubject},{'dip','Vdata','Vmodel'});
    dip_4_tc{isubject} = dip_4{isubject};    
    dip_4_tc{isubject}.avg = dip_4_tc{isubject}.dip.mom;    
    dip_4_tc{isubject}.label = {'dip1','dip2','dip3'};
    dip_4_tc{isubject} = rmfield(dip_4_tc{isubject},{'dip','Vdata','Vmodel'});         
end

for isubject = slist
    cfg = [];
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [1 30];
    dip_1_tc{isubject} = ft_preprocessing(cfg,dip_1_tc{isubject} );
    dip_2_tc{isubject} = ft_preprocessing(cfg,dip_2_tc{isubject} );
    dip_3_tc{isubject} = ft_preprocessing(cfg,dip_3_tc{isubject} );
    dip_4_tc{isubject} = ft_preprocessing(cfg,dip_4_tc{isubject} );
end

% normalize data over trials per subject
for isubject = slist
    dip_cat = [dip_1_tc{isubject}.avg(1,:), dip_2_tc{isubject}.avg(1,:), dip_3_tc{isubject}.avg(1,:), dip_4_tc{isubject}.avg(1,:)];
    dip_cat = zscore(dip_cat);
    
    dip_1_norm{isubject} = dip_1_tc{isubject}; 
    dip_2_norm{isubject} = dip_2_tc{isubject}; 
    dip_3_norm{isubject} = dip_3_tc{isubject}; 
    dip_4_norm{isubject} = dip_4_tc{isubject}; 
    
    dip_1_norm{isubject}.avg = dip_cat(1:11000);
    dip_2_norm{isubject}.avg = dip_cat(11001:22000);
    dip_3_norm{isubject}.avg = dip_cat(22001:33000);
    dip_4_norm{isubject}.avg = dip_cat(33001:44000);
end

% 
% for isubject = slist
%     dip_1_tc{isubject}.avg  = dip_1_tc{isubject}.avg(:,end-10000+1:end);
%     dip_2_tc{isubject}.avg  = dip_2_tc{isubject}.avg(:,end-10000+1:end);
%     dip_3_tc{isubject}.avg  = dip_3_tc{isubject}.avg(:,end-10000+1:end);
%     dip_4_tc{isubject}.avg  = dip_4_tc{isubject}.avg(:,end-10000+1:end);
%     dip_1_tc{isubject}.time  = dip_1_tc{isubject}.time(:,end-10000+1:end);
%     dip_2_tc{isubject}.time  = dip_2_tc{isubject}.time(:,end-10000+1:end);
%     dip_3_tc{isubject}.time  = dip_3_tc{isubject}.time(:,end-10000+1:end);
%     dip_4_tc{isubject}.time  = dip_4_tc{isubject}.time(:,end-10000+1:end);
% end

clear dip_1 dip_2 dip_3 dip_4

%% put in data matrix
i = 1;
dat_1 = nan(length(slist),11000);
dat_2 = nan(length(slist),11000);
dat_3 = nan(length(slist),11000);
dat_4 = nan(length(slist),11000);
    
for isubject = slist
    dat_1(i,:) = dip_1_norm{isubject}.avg(1,:);
    dat_2(i,:) = dip_2_norm{isubject}.avg(1,:);
    dat_3(i,:) = dip_3_norm{isubject}.avg(1,:);
    dat_4(i,:) = dip_4_norm{isubject}.avg(1,:);
    i = i + 1;  
end

temp(:,:,1) = dat_1;
temp(:,:,2) = dat_2;
temp(:,:,3) = dat_3;
temp(:,:,4) = dat_4;
dat_avg = nanmean(temp,3);
clear temp

%% orient according to correlation with mean

time      = linspace(-10,1,10000) - 0.033 - 1/16; % t=0 at first nonstimulus

%%%%%%%%%%%%%%%%% ADAPT TO NEW TRIALLENGTH
latency   = 9150:9450; % define in time
%%%%%%%%%%%%%%%%%

temp_avg  = mean(dat_avg,1);
flip      = ones(1,size(dat_avg,1));

for iteration = 1 : 10
    disp(['orienting iteration',num2str(iteration)]);
    for isubject = 1 : size(dat_avg,1)
        if corr(dat_avg(isubject,latency)',temp_avg(latency)','rows','complete') < 0
            flip(isubject) = -flip(isubject);
            dat_avg(isubject,:)     = -dat_avg(isubject,:);
            dat_1(isubject,:)       = -dat_1(isubject,:);
            dat_2(isubject,:)       = -dat_2(isubject,:);
            dat_3(isubject,:)       = -dat_3(isubject,:);
            dat_4(isubject,:)       = -dat_4(isubject,:);
        end        
    end
    temp_avg = mean(dat_avg,1);
end

figure; imagesc(dat_avg(:,latency))

dat_avg_GA = nanmean(dat_avg,1);
dat_std_GA = std(dat_avg,1);
figure; plot(time, dat_avg_GA); % first peak at 9229 = , second peak at 9307 = 

dat_1_GA = nanmean(dat_1,1);
std_1_GA = std(dat_1,1) /sqrt(22);
dat_2_GA = nanmean(dat_2,1);
std_2_GA = std(dat_2,1) /sqrt(22);
dat_3_GA = nanmean(dat_3,1);
std_3_GA = std(dat_3,1) /sqrt(22);
dat_4_GA = nanmean(dat_4,1);
std_4_GA = std(dat_4,1) /sqrt(22);

figure; plot(time, [dat_1_GA; dat_2_GA; dat_3_GA; dat_4_GA]);

[~, it1] = find(time>-0.6,1,'first');
[~, it2] = find(time>0.4,1,'first');
%%
figure; hold;
c = jet(4);
c = c(4:-1:1,:);
fill([time(it1:it2) time(it2:-1:it1)],[dat_1_GA(it1:it2)+std_1_GA(it1:it2), dat_1_GA(it2:-1:it1)-std_1_GA(it2:-1:it1)],c(1,:),'linestyle','none','facealpha',0.5);
fill([time(it1:it2) time(it2:-1:it1)],[dat_2_GA(it1:it2)+std_1_GA(it1:it2), dat_2_GA(it2:-1:it1)-std_2_GA(it2:-1:it1)],c(2,:),'linestyle','none','facealpha',0.5);
fill([time(it1:it2) time(it2:-1:it1)],[dat_3_GA(it1:it2)+std_1_GA(it1:it2), dat_3_GA(it2:-1:it1)-std_3_GA(it2:-1:it1)],c(3,:),'linestyle','none','facealpha',0.5);
fill([time(it1:it2) time(it2:-1:it1)],[dat_4_GA(it1:it2)+std_1_GA(it1:it2), dat_4_GA(it2:-1:it1)-std_4_GA(it2:-1:it1)],c(4,:),'linestyle','none','facealpha',0.5);

plot(time(it1:it2),dat_1_GA(it1:it2),'Color',c(1,:),'Linewidth',1);
plot(time(it1:it2),dat_2_GA(it1:it2),'Color',c(2,:),'Linewidth',1);
plot(time(it1:it2),dat_3_GA(it1:it2),'Color',c(3,:),'Linewidth',1);
plot(time(it1:it2),dat_4_GA(it1:it2),'Color',c(4,:),'Linewidth',1);

axis tight
title('Offset response')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend({'high',' ',' ','low'},'location','northwest');
print('-dsvg','w:\WANDER\images\article\SSEF_offset.svg');
print('-dpng','w:\WANDER\images\article\SSEF_offset.png');


%% find and test peaks
i1 = find(time > 0.040,1,'first');
i2 = find(time > 0.060,1,'first');
[~, peak_i0] = max(dat_avg_GA(i1:i2));
peak_i0 = peak_i0+i1-1;
peak_t0 = time(peak_i0);

i1 = find(time > 0.1,1,'first');
i2 = find(time > 0.180,1,'first');
[~, peak_i1] = min(dat_avg_GA(i1:i2));
peak_i1 = peak_i1+i1-1;
peak_t1 = time(peak_i1);

i1 = find(time > 0.15,1,'first');
i2 = find(time > 0.25,1,'first');
[~, peak_i2] = max(dat_avg_GA(i1:i2));
peak_i2 = peak_i2+i1-1;
peak_t2 = time(peak_i2);


[~,~,stats] = anova1([dat_1(:,peak_i1), dat_2(:,peak_i1), dat_3(:,peak_i1), dat_4(:,peak_i1)]);
[c,~,~,gnames] = multcompare(stats);

[~,~,stats] = anova1([dat_1(:,peak_i0), dat_2(:,peak_i0), dat_3(:,peak_i0), dat_4(:,peak_i0)]);
[c,~,~,gnames] = multcompare(stats);

[~,~,stats] = anova1([dat_1(:,peak_i2), dat_2(:,peak_i2), dat_3(:,peak_i2), dat_4(:,peak_i2)]);
[c,~,~,gnames] = multcompare(stats);

[H,P] = ttest([dat_1(:,peak_i1), dat_4(:,peak_i1)]);
[H,P] = ttest([dat_1(:,peak_i2), dat_4(:,peak_i2)]);

%% FFT
latency = [find(time > -10, 1, 'first') : find(time > 0, 1, 'first')];

Fs = 1000;
L = length(latency);
f = Fs*(0:(L/2))/L;
i = 1;

FFT_1            = fft(dat_1(:,latency)');
FFT_1            = abs(FFT_1/L)';
FFT_1            = FFT_1(:,1:L/2+1)';
FFT_1(2:end-1,:) = 2*FFT_1(2:end-1,:);
FFT_1_avg = mean(FFT_1,2)';
FFT_1_std = std(FFT_1')/sqrt(22);

FFT_2            = fft(dat_2(:,latency)');
FFT_2            = abs(FFT_2/L)';
FFT_2            = FFT_2(:,1:L/2+1)';
FFT_2(2:end-1,:) = 2*FFT_2(2:end-1,:);
FFT_2_avg = mean(FFT_2,2)';
FFT_2_std = std(FFT_2')/sqrt(22);

FFT_3            = fft(dat_3(:,latency)');
FFT_3            = abs(FFT_3/L)';
FFT_3            = FFT_3(:,1:L/2+1)';
FFT_3(2:end-1,:) = 2*FFT_3(2:end-1,:);
FFT_3_avg = mean(FFT_3,2)';
FFT_3_std = std(FFT_3')/sqrt(22);

FFT_4            = fft(dat_4(:,latency)');
FFT_4            = abs(FFT_4/L)';
FFT_4            = FFT_4(:,1:L/2+1)';
FFT_4(2:end-1,:) = 2*FFT_4(2:end-1,:);
FFT_4_avg = mean(FFT_4,2)';
FFT_4_std = std(FFT_4')/sqrt(22);

FFT_avg            = fft(dat_4(:,latency)');
FFT_avg            = abs(FFT_avg/L)';
FFT_avg            = FFT_avg(:,1:L/2+1)';
FFT_avg(2:end-1,:) = 2*FFT_avg(2:end-1,:);
FFT_avg_avg = mean(FFT_avg,2)';
FFT_avg_std = std(FFT_avg')/sqrt(22);


[~, ilow] = find(f>13,1,'first');
[~, ihigh] = find(f>19,1,'first');

figure; hold;
c = jet(4);
c = c(4:-1:1,:);
% fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_1_avg(ilow:ihigh)+FFT_1_std(ilow:ihigh), FFT_1_avg(ihigh:-1:ilow)-FFT_1_std(ihigh:-1:ilow)],c(1,:),'linestyle','none','facealpha',0.5);
% fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_2_avg(ilow:ihigh)+FFT_2_std(ilow:ihigh), FFT_2_avg(ihigh:-1:ilow)-FFT_2_std(ihigh:-1:ilow)],c(2,:),'linestyle','none','facealpha',0.5);
% fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_3_avg(ilow:ihigh)+FFT_3_std(ilow:ihigh), FFT_3_avg(ihigh:-1:ilow)-FFT_3_std(ihigh:-1:ilow)],c(3,:),'linestyle','none','facealpha',0.5);
% fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_4_avg(ilow:ihigh)+FFT_4_std(ilow:ihigh), FFT_4_avg(ihigh:-1:ilow)-FFT_4_std(ihigh:-1:ilow)],c(4,:),'linestyle','none','facealpha',0.5);
plot(f(ilow:ihigh),FFT_1_avg(ilow:ihigh),'Color',c(1,:),'Linewidth',1);
plot(f(ilow:ihigh),FFT_2_avg(ilow:ihigh),'Color',c(2,:),'Linewidth',1);
plot(f(ilow:ihigh),FFT_3_avg(ilow:ihigh),'Color',c(3,:),'Linewidth',1);
plot(f(ilow:ihigh),FFT_4_avg(ilow:ihigh),'Color',c(4,:),'Linewidth',1);
axis tight;
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend({'high',' ',' ','low'});
print('-dsvg','w:\WANDER\images\article\SSEF_power_line.svg');
print('-dpng','w:\WANDER\images\article\SSEF_power_line.png');

figure; hold;
c = jet(4);
c = c(4:-1:1,:);
fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_1_avg(ilow:ihigh)+FFT_1_std(ilow:ihigh), FFT_1_avg(ihigh:-1:ilow)-FFT_1_std(ihigh:-1:ilow)],c(1,:),'linestyle','none','facealpha',0.5);
fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_2_avg(ilow:ihigh)+FFT_2_std(ilow:ihigh), FFT_2_avg(ihigh:-1:ilow)-FFT_2_std(ihigh:-1:ilow)],c(2,:),'linestyle','none','facealpha',0.5);
fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_3_avg(ilow:ihigh)+FFT_3_std(ilow:ihigh), FFT_3_avg(ihigh:-1:ilow)-FFT_3_std(ihigh:-1:ilow)],c(3,:),'linestyle','none','facealpha',0.5);
fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_4_avg(ilow:ihigh)+FFT_4_std(ilow:ihigh), FFT_4_avg(ihigh:-1:ilow)-FFT_4_std(ihigh:-1:ilow)],c(4,:),'linestyle','none','facealpha',0.5);
% plot(f(ilow:ihigh),FFT_1_avg(ilow:ihigh),'Color',c(1,:),'Linewidth',1);
% plot(f(ilow:ihigh),FFT_2_avg(ilow:ihigh),'Color',c(2,:),'Linewidth',1);
% plot(f(ilow:ihigh),FFT_3_avg(ilow:ihigh),'Color',c(3,:),'Linewidth',1);
% plot(f(ilow:ihigh),FFT_4_avg(ilow:ihigh),'Color',c(4,:),'Linewidth',1);
axis tight
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend({'high',' ',' ','low'});
print('-dsvg','w:\WANDER\images\article\SSEF_power_fill.svg');
print('-dpng','w:\WANDER\images\article\SSEF_power_fill.png');


figure; hold;
c = jet(4);
c = c(4:-1:1,:);
fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_1_avg(ilow:ihigh)+FFT_1_std(ilow:ihigh), FFT_1_avg(ihigh:-1:ilow)-FFT_1_std(ihigh:-1:ilow)],c(1,:),'linestyle','none','facealpha',0.5);
fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_2_avg(ilow:ihigh)+FFT_2_std(ilow:ihigh), FFT_2_avg(ihigh:-1:ilow)-FFT_2_std(ihigh:-1:ilow)],c(2,:),'linestyle','none','facealpha',0.5);
fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_3_avg(ilow:ihigh)+FFT_3_std(ilow:ihigh), FFT_3_avg(ihigh:-1:ilow)-FFT_3_std(ihigh:-1:ilow)],c(3,:),'linestyle','none','facealpha',0.5);
fill([f(ilow:ihigh) f(ihigh:-1:ilow)],[FFT_4_avg(ilow:ihigh)+FFT_4_std(ilow:ihigh), FFT_4_avg(ihigh:-1:ilow)-FFT_4_std(ihigh:-1:ilow)],c(4,:),'linestyle','none','facealpha',0.5);
plot(f(ilow:ihigh),FFT_1_avg(ilow:ihigh),'Color',c(1,:),'Linewidth',1);
plot(f(ilow:ihigh),FFT_2_avg(ilow:ihigh),'Color',c(2,:),'Linewidth',1);
plot(f(ilow:ihigh),FFT_3_avg(ilow:ihigh),'Color',c(3,:),'Linewidth',1);
plot(f(ilow:ihigh),FFT_4_avg(ilow:ihigh),'Color',c(4,:),'Linewidth',1);
axis tight
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
legend({'high',' ',' ','low'});
print('-dsvg','w:\WANDER\images\article\SSEF_power_fill_line.svg');
print('-dpng','w:\WANDER\images\article\SSEF_power_fill_line.png');


%% find and test peaks
[y fi] = max(FFT_avg_avg)

[~,~,stats] = anova1([FFT_1(fi,:)', FFT_2(fi,:)', FFT_3(fi,:)', FFT_4(fi,:)']);
[c,~,~,gnames] = multcompare(stats);

[H,P] = ttest([FFT_1(fi,:)', FFT_4(fi,:)'])

%% correlate FFT with peak amplitude
[y fi] = max(FFT_avg_avg)

clear c
for i = 1 : size(FFT_1,2)
    temp = corr([FFT_1(fi,i), dat_1(i,peak_i1); FFT_2(fi,i),  dat_2(i,peak_i1); FFT_3(fi,i), dat_3(i,peak_i1); FFT_4(fi,i), dat_4(i,peak_i1)],'type','spearman');
    c(i) = temp(1,2);
end
[h p] = ttest(c)

clear c
for i = 1 : size(FFT_1,2)
    temp = corr([FFT_1(fi,i), dat_1(i,peak_i2); FFT_2(fi,i),  dat_2(i,peak_i2); FFT_3(fi,i), dat_3(i,peak_i2); FFT_4(fi,i), dat_4(i,peak_i2)],'type','spearman');
    c(i) = temp(1,2);
end
[h p] = ttest(c)

%% correlate mean power between 15-17Hz with peak amplitude 

fi = [find(f>15,1,'first'):find(f<17,1,'last')];
% fi = [find(f>15,1,'first'):find(f>17,1,'first')];

clear c
for i = 1 : size(FFT_1,2)
    temp = corr([mean(FFT_1(fi,i)), dat_1(i,peak_i1); mean(FFT_2(fi,i)),  dat_2(i,peak_i1); mean(FFT_3(fi,i)), dat_3(i,peak_i1); mean(FFT_4(fi,i)), dat_4(i,peak_i1)],'type','spearman');
    c(i) = temp(1,2);
end
[h p] = ttest(c)

clear c
for i = 1 : size(FFT_1,2)
    temp = corr([mean(FFT_1(fi,i)), dat_1(i,peak_i2); mean(FFT_2(fi,i)),  dat_2(i,peak_i2); mean(FFT_3(fi,i)), dat_3(i,peak_i2); mean(FFT_4(fi,i)), dat_4(i,peak_i2)],'type','spearman');
    c(i) = temp(1,2);
end
[h p] = ttest(c)


%% for each freq

clear c
for fi = 1 : size(FFT_1,1)
    for i = 1 : size(FFT_1,2)
        temp = corr([FFT_1(fi,i), dat_1(i,peak_i2); FFT_2(fi,i),  dat_2(i,peak_i2); FFT_3(fi,i), dat_3(i,peak_i2); FFT_4(fi,i), dat_4(i,peak_i2)],'type','spearman');
        c(fi,i) = temp(1,2);
    end
end
figure; plot(f,mean(c,2))
[h p] = ttest(c')

%% for each timepoint, at 16hz
[y fi] = max(FFT_avg_avg)

[y ti] = find(time > 0,1,'first');
ti = 1;

clear c
tii = 1;
for t = ti : length(time)
    for isubject = 1 : size(FFT_1,2)
        temp = corr([FFT_1(fi,isubject), dat_1(isubject,t); FFT_2(fi,isubject),  dat_2(isubject,t); FFT_3(fi,isubject), dat_3(isubject,t); FFT_4(fi,isubject), dat_4(isubject,t)],'type','spearman');
        c(tii,isubject) = temp(1,2);
    end
    tii = tii + 1;
end

figure; plot(time(ti:end),mean(c,2))
[h p] = ttest(c')
figure; plot(time(ti:end),p)

figure; hold;
plot(time(ti:end),h)
plot(time(ti:end),dat_avg_GA(ti:end))

%% save data for correlations with other data
save('FFT_ERF_data_for_correlations','dat*','FFT*','peak*','time','f');

%% now to extract single-trial data from dipole
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

load('dipoles.mat'); 
slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

for isubject = slist
 
    
    % load headmodel
    [mri_segmented{isubject}, headmodel{isubject}, subject_grid{isubject}, template_grid, ~] = WANDER_grid(isubject,1,0);
    
    % load sensor information
    [dataset_task, ~] = WANDER_subjectinfo;
    hdr{isubject} = ft_read_header(dataset_task{isubject,1});
    
    % load timcourse data
    data_epoch_MEG          = WANDER_ICA(isubject,0,1,0);  
    artdef                  = WANDER_artefact_detection_MEG(isubject,0,1,0);
    
    % now lock to end of trial
    for iblock = 1 : 4
        
        % add trialnum
        for itrial = 1 : size(data_epoch_MEG{iblock}.trialinfo,1)
            data_epoch_MEG{iblock}.trialinfo(itrial,7) = itrial + (iblock-1)*50;
        end
        
        % select magnetometers
        cfg = [];
        cfg.channel = 'MEG*1';
        data_epoch_MEG{iblock} = ft_selectdata(cfg,data_epoch_MEG{iblock});
        
        % filter
        cfg = [];
        cfg.bpfilter    = 'yes';
        cfg.bpfreq      = [1 30];
        data_epoch_MEG{iblock} = ft_preprocessing(cfg,data_epoch_MEG{iblock});    
        
        % start trial at first stimulation
        for itrial = 1 : size(data_epoch_MEG{iblock}.trial,2)
            fprintf('Cutting off beginning of trial %d block %d \n',itrial,iblock);
            data_epoch_MEG{iblock}.trial{itrial} = data_epoch_MEG{iblock}.trial{itrial}(:,1500:end);
            data_epoch_MEG{iblock}.time{itrial}  = data_epoch_MEG{iblock}.time{itrial}(1500:end);
        end
        
        % remove sampleinfo which else creates unexpected behaviour in
        % rejectartefact
        data_epoch_MEG{iblock} = rmfield(data_epoch_MEG{iblock},'sampleinfo');
        
        % remove trials with artefacts
        cfg = [];
        cfg.artfctdef                   = artdef{iblock};        
        cfg.artfctdef.reject            = 'complete';
        data_epoch_MEG{iblock}          = ft_rejectartifact(cfg,data_epoch_MEG{iblock});
         
        % lock data to end of trial
        for itrial = 1 : size(data_epoch_MEG{iblock}.trial,2)
            fprintf('Shifting timecourse of trial %d block %d \n',itrial,iblock);
            keeptime = find(data_epoch_MEG{iblock}.time{itrial} <= 0.5);
            data_epoch_MEG{iblock}.time{itrial}  = data_epoch_MEG{iblock}.time{itrial} - data_epoch_MEG{iblock}.time{itrial}(end) + 1.000;           
        end
        
        % select last 10 seconds of trial + second after offset
        cfg = [];
        cfg.latency = [-9.8, 0.6];
        data_epoch_MEG{iblock} = ft_selectdata(cfg,data_epoch_MEG{iblock});    
    end

    % concatinate blocks
    data_comb = ft_appenddata([],data_epoch_MEG{:});
    clear data_epoch_MEG
    
    % extract timecourse from dipole
    dip = data_comb;
    dip.label = {'ori1','ori2','ori3'};

    for itrial = 1 : size(data_comb.trialinfo,1)
        
        cfg                 = [];
        cfg.trials          = itrial;
        data_temp           = ft_selectdata(cfg,data_comb);
        
        data_temp.avg       = data_temp.trial{1};
        data_temp.time      = data_temp.time{1};        
        data_temp           = rmfield(data_temp,'trial');
        
        cfg                 = [];
        cfg.latency         = 'all';
        cfg.numdipoles      = 1;
        cfg.symmetry        = [];
        cfg.nonlinear       = 'no';  % use a fixed position
        cfg.gridsearch      = 'no';
        cfg.dip.pos         = source_mag{isubject}.dip.pos;
        cfg.vol             = headmodel{isubject};
        cfg.channel         = {'MEG*1'};
        cfg.senstype        = 'meg';
        cfg.grad            = hdr{isubject}.grad;
        temp                = ft_dipolefitting(cfg, data_temp);
        dip.trial{itrial}   = temp.dip.mom;
        
        % FFT
        latency = 1 : size(dip.trial{itrial},2)-1000;        
        Fs = 1000;
        L = length(latency);
        f = Fs*(0:(L/2))/L;
        
        trial_FFT            = fft(dip.trial{itrial}(:,latency)');
        trial_FFT            = abs(trial_FFT/L)';
        trial_FFT            = trial_FFT(:,1:L/2+1);
        trial_FFT(2:end-1,:) = 2*trial_FFT(2:end-1,:);
        
        dip.pow{itrial} = trial_FFT;
        dip.freq = f;
    end
        
    fname_ERF_ratings_singletrial = ['W:\WANDER\data\ERF\s' num2str(isubject) '_ratings_singletrial.mat'];
    save(fname_ERF_ratings_singletrial,'dip');
end
% 
% for isubject = slist
%     dip_1_tc{isubject} = dip_1{isubject};    
%     dip_1_tc{isubject}.avg = dip_1_tc{isubject}.dip.mom;    
%     dip_1_tc{isubject}.label = {'dip1','dip2','dip3'};
%     dip_1_tc{isubject} = rmfield(dip_1_tc{isubject},{'dip','Vdata','Vmodel'});
%     dip_2_tc{isubject} = dip_2{isubject};    
%     dip_2_tc{isubject}.avg = dip_2_tc{isubject}.dip.mom;    
%     dip_2_tc{isubject}.label = {'dip1','dip2','dip3'};
%     dip_2_tc{isubject} = rmfield(dip_2_tc{isubject},{'dip','Vdata','Vmodel'});
%     dip_3_tc{isubject} = dip_3{isubject};    
%     dip_3_tc{isubject}.avg = dip_3_tc{isubject}.dip.mom;    
%     dip_3_tc{isubject}.label = {'dip1','dip2','dip3'};
%     dip_3_tc{isubject} = rmfield(dip_3_tc{isubject},{'dip','Vdata','Vmodel'});
%     dip_4_tc{isubject} = dip_4{isubject};    
%     dip_4_tc{isubject}.avg = dip_4_tc{isubject}.dip.mom;    
%     dip_4_tc{isubject}.label = {'dip1','dip2','dip3'};
%     dip_4_tc{isubject} = rmfield(dip_4_tc{isubject},{'dip','Vdata','Vmodel'});         
% end
% 
% for isubject = slist
%     cfg = [];
%     cfg.bpfilter    = 'yes';
%     cfg.bpfreq      = [1 30];
%     dip_1_tc{isubject} = ft_preprocessing(cfg,dip_1_tc{isubject} );
%     dip_2_tc{isubject} = ft_preprocessing(cfg,dip_2_tc{isubject} );
%     dip_3_tc{isubject} = ft_preprocessing(cfg,dip_3_tc{isubject} );
%     dip_4_tc{isubject} = ft_preprocessing(cfg,dip_4_tc{isubject} );
% end
% 
% % normalize data over trials per subject
% for isubject = slist
%     dip_cat = [dip_1_tc{isubject}.avg(1,:), dip_2_tc{isubject}.avg(1,:), dip_3_tc{isubject}.avg(1,:), dip_4_tc{isubject}.avg(1,:)];
%     dip_cat = zscore(dip_cat);
%     
%     dip_1_norm{isubject} = dip_1_tc{isubject}; 
%     dip_2_norm{isubject} = dip_2_tc{isubject}; 
%     dip_3_norm{isubject} = dip_3_tc{isubject}; 
%     dip_4_norm{isubject} = dip_4_tc{isubject}; 
%     
%     dip_1_norm{isubject}.avg = dip_cat(1:11000);
%     dip_2_norm{isubject}.avg = dip_cat(11001:22000);
%     dip_3_norm{isubject}.avg = dip_cat(22001:33000);
%     dip_4_norm{isubject}.avg = dip_cat(33001:44000);
% end
% 
