slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

rootpath = 1;
force = 0;

for isubject = slist
    [mri_segmented_brain{isubject}, headmodel{isubject}, subject_grid{isubject}, template_grid{isubject}, mri_realigned{isubject}] = WANDER_grid(isubject,rootpath,force);
end


iplot = 1;
figure;
hold on;
for isubject = slist
    subplot(5,5,iplot);
    % make a figure of the single subject headmodel, and grid positions
    ft_plot_vol(headmodel{isubject}, 'edgecolor', 'none', 'facecolor', 'brain','facealpha', 0.4);
    ft_plot_mesh(subject_grid{isubject}.pos(subject_grid{isubject}.inside,:));
    axis tight
    title(num2str(isubject));
    az = 90;
    el = 0;
    view(az, el); 
    iplot = iplot + 1;
end

isubject = 23;
figure;
  ft_plot_vol(headmodel{isubject}, 'edgecolor', 'none', 'facecolor', 'brain','facealpha', 0.4);
    ft_plot_mesh(subject_grid{isubject}.pos(subject_grid{isubject}.inside,:));


% transparent
iplot = 1;
figure;
hold on;
for isubject = slist
    subplot(5,5,iplot);
    mri = imagesc(flipud(squeeze(mri_realigned{isubject}.anatomy(100,:,:))'));
    hold;
    brain = image(flipud(squeeze(mri_segmented_brain{isubject}.brain(100,:,:))')*600);
    alpha = zeros(size(mri.CData));
    alpha(find(brain.CData)) = 0.5;
    set(brain, 'AlphaData', alpha);
    axis square
    iplot = iplot + 1;
    title(isubject);
end
colormap('jet')

% oblique
iplot = 1;
figure;
hold on;
for isubject = slist
    subplot(5,5,iplot);
    mri = (flipud(squeeze(mri_realigned{isubject}.anatomy(100,:,:))'));
    mri(find(flipud(squeeze(mri_segmented_brain{isubject}.brain(100,:,:))'))) = 1000;
    imagesc(mri);
    axis square
    iplot = iplot + 1;
    title(isubject);
end
colormap('jet')

% cut-outs
iplot = 1;
figure;
hold on;
for isubject = slist
    subplot(5,5,iplot);
    mri = (flipud(squeeze(mri_realigned{isubject}.anatomy(100,:,:))'));
    brain = (flipud(squeeze(mri_segmented_brain{isubject}.brain(100,:,:))'));
    mri(~brain) = 0;
    imagesc(mri);
    axis square
    iplot = iplot + 1;
    title(isubject);
end
colormap('jet')

