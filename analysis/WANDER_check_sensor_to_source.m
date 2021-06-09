isubject = 11;

% [mri_segmented, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(isubject,1,0);
cfg.location = [40 -35 55]; % from nr. 11

loc = find(template_grid.pos(:,1) == cfg.location(1) & template_grid.pos(:,2) == cfg.location(2) & template_grid.pos(:,3) == cfg.location(3));
loc = 37117;

fig = figure; 
subplot(2,1,1); hold;
plot(source10_powrpt.pow(loc,:)); 
plot(source11_powrpt.pow(loc,:));
title('max alpha diff voxel');
ylabel('alpha power');
xlabel('observation (trials x time)');
axis tight;
subplot(2,1,2); hold;
plot(FFT.powspctrm(:,41,10));
plot(FFT.powspctrm(:,41,11));
title('max alpha diff sensor');
ylabel('alpha power');
xlabel('observation (trials x time)');
axis tight;
print(fig,'-dpdf','W:/WANDER/images/comparison_source_sens.pdf');

fig = figure;
subplot(2,1,1); hold;
scatter(log(source10_powrpt.pow(loc,:)),log(FFT.powspctrm(:,41,10))','.');
scatter(log(source11_powrpt.pow(loc,:)),log(FFT.powspctrm(:,41,11))','.');
[coef10, pval10] = corr(log(source10_powrpt.pow(loc,:))',log(FFT.powspctrm(:,41,10)));
[coef11, pval11] = corr(log(source10_powrpt.pow(loc,:))',log(FFT.powspctrm(:,41,10)));
lsline
title(['Correlation (source x sensor) 10Hz (blue): ' num2str(coef10) ', 11Hz (orange): ' num2str(coef11) ]);
xlabel('source log(alpha power)');
ylabel('sensor log(alpha power)');

subplot(2,1,2);
hold
[C10,LAGS10] = xcorr(source10_powrpt.pow(loc,:),FFT.powspctrm(:,41,10)');
[C11,LAGS11] = xcorr(source11_powrpt.pow(loc,:),FFT.powspctrm(:,41,11)');
plot(LAGS10,C10);
plot(LAGS11,C11);
title('crosscorrelation alpha power source x sensor, 10Hz (blue), 11Hz (orange)');
print(fig,'-dpdf','W:/WANDER/images/comparison_source_sens_corr.pdf');

figure;
subplot(2,1,1);
title(['Correlation coefficient: ' num2str(coef)]);
scatter(source10_powrpt.pow(loc,:),source11_powrpt.pow(loc,:),'.');

[coef, pval] = corr(source10_powrpt.pow(loc,:)',FFT.powspctrm(:,41,10));
lsline
title(['Correlation coefficient: ' num2str(coef)]);
xlabel('source alpha power');
ylabel('sensor alpha power');
subplot(2,1,2);
[C,LAGS] = xcorr(source10_powrpt.pow(loc,:)',FFT.powspctrm(:,41,10));
plot(LAGS,C);
title('crosscorrelation');


save('W:\WANDER\images\check_sensor_to_source','source10_powrpt','source11_powrpt');
load('W:\WANDER\images\check_sensor_to_source');

cfg.channel = 'MEG1132+1133';
cfg.channel = 41;

 
