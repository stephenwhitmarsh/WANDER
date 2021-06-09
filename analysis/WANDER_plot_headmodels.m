

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD
slist = [1 3:5 8:11 13 17:20 22:26]; %without subjects with more than 2 SD

figure;
i = 1;
for isubject = slist
    [mri_segmented_brain{isubject}, headmodel{isubject}, subject_grid{isubject}, template_grid{isubject}, mri_realigned{isubject}] = WANDER_grid(isubject,0);
end

for isubject = slist

    subplot(5,5,i);
    ft_plot_vol(headmodel{isubject}, 'edgecolor', 'none', 'facecolor', 'brain','facealpha', 0.4);
    ft_plot_mesh(subject_grid{isubject}.pos(subject_grid{isubject}.inside,:));
    view(50,0);
    axis tight;
%     saveas(gcf,[outputdir_image filesep num2str(isubject) '_headmodelplot' ],'png');
    
    i = i + 1;
end


for isubject = slist
    WANDER_source_MI(isubject,1,0)
end
