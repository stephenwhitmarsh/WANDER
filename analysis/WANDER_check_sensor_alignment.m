slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD
rootpath = 1;
force = 0;
figure;
i = 1;
for isubject = slist
    
    
    [mri_segmented_brain, headmodel, subject_grid, template_grid, mri_realigned] = WANDER_grid(isubject,rootpath,force);
    
    % load sensor information
    [dataset_task, ~] = WANDER_subjectinfo(rootpath);
    hdr = ft_read_header(dataset_task{isubject,1});
    
    % make a figure of the single subject headmodel, and grid positions
    %     figure; hold on;
    subplot(5,5,i);
    ft_plot_vol(headmodel, 'edgecolor', 'none', 'facecolor', 'brain','facealpha', 0.4);
    %     ft_plot_mesh(subject_grid.pos(subject_grid.inside,:),'vertexsize',2,'vertexmarker','*');
    ft_plot_sens(hdr.grad,'unit','mm');
    view(180,0);
    title(num2str(isubject));
    %     saveas(gcf,['D:\temp\headmodelplot_' num2str(isubject)],'png');
    i = i +1;
    
end
    print('Z:\WANDER\headmodel_sens_plot_overview.png','-dpng','-r600')
