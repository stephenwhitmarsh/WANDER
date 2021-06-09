 TFR_phaselocked = WANDER_TFR_phaselocking_balanced(19,0,'cue',[1 30],1); % max MI
 TFR_phaselocked = WANDER_TFR_phaselocking_balanced(5,0,'cue',[1 30],1); % max diff MI

 cfg = [];
%  cfg.channel =  ft_channelselection(find(stat_MI.mask), stat_MI.label);
 
cfg.avgoverchan = 'yes';
PAC_all = ft_selectdata( cfg,TFR_phaselocked.PPP);
PAC_high = ft_selectdata( cfg,TFR_phaselocked.PPP_high);
PAC_low = ft_selectdata( cfg,TFR_phaselocked.PPP_low);

MI_all = ft_selectdata( cfg,TFR_phaselocked.MI{1});
MI_high = ft_selectdata( cfg,TFR_phaselocked.MI_high{1});
MI_low = ft_selectdata( cfg,TFR_phaselocked.MI_low{1});

PAC_high.powspctrm_norm = PAC_high.powspctrm ./ mean(PAC_high.powspctrm,3);
PAC_low.powspctrm_norm = PAC_low.powspctrm ./ mean(PAC_low.powspctrm,3);



%%

fig = figure;
subplot(1,2,1);

colormap(redblue(256));
imagesc(squeeze(PAC_all.powspctrm),[0.8 1.2]);
set(gca,'XTick',[1:18])
title('all')
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal');
colorbar

subplot(1,2,2);
plot(MI_all.avg);
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal');
xlim([0.5 30.5]);
ylim([0 3e-3]);
view([90 -90]);

savefig(fig,'w:\WANDER\images\poster\MI_PAC_s19');

%%

figure;
colormap(redblue(256));
imagesc(squeeze(PAC_all.powspctrm),[0.7 1.3]);
set(gca,'XTick',[1:18])
title('all')
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal');

figure;
subplot(2,2,1);
colormap(redblue(256));

imagesc(squeeze(PAC_high.powspctrm_norm),[0.7 1.3]);
title('low attention')
colorbar('westoutside');
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal','XTick',[1:18]);

subplot(2,2,2);
plot(MI_high.avg);
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal');
xlim([0.5 30.5]);
ylim([0 6e-3]);
view([90 -90]);

subplot(2,2,3);
imagesc(squeeze(PAC_low.powspctrm_norm),[0.7 1.3]);
title('high attention')
colorbar('westoutside');
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal','XTick',[1:18]);

subplot(2,2,4);
plot(MI_low.avg);
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal');
xlim([0.5 30.5]);
ylim([0 6e-3]);
view([90 -90]);

%%
figure;
 colormap(redblue(256));

subplot(1,4,1);
imagesc(squeeze(PAC_high.powspctrm_norm),[0.7 1.3]);
% title('low attention')
colorbar('westoutside');
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal','XTick',[1:18]);

subplot(1,4,2);
imagesc(squeeze(PAC_low.powspctrm_norm),[0.7 1.3]);
% title('high attention')
colorbar('westoutside');
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal','XTick',[1:18]);

subplot(1,4,3);
plot(MI_low.avg,'r'); hold;
plot(MI_high.avg,'b');
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal');
xlim([0.5 30.5]);
ylim([0 6e-3]);
view([90 -90]);

subplot(1,4,4);
plot((MI_high.avg - MI_low.avg) ./ (MI_high.avg + MI_low.avg),'b');
set(gca,'PlotBoxAspectRatio',[1,1,1],'YDir','normal');
xlim([0.5 30.5]);
ylim([-0.3 0.3]);
view([90 -90]);

