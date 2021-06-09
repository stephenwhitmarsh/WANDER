

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

restingstate    = 0;
rootpath        = 1;
force           = 0;
latency         = 'all';
timing          = 'cue';

for isubject = slist
    
    % load TFR
%     temp = WANDER_TFR(isubject,0,'cue',1,0);
%     TFR{isubject} = temp{1};
%     TFR{isubject}.powspctrm = [temp{1}.powspctrm; temp{2}.powspctrm; temp{3}.powspctrm; temp{4}.powspctrm];
%      
    % load heart beat data
    [HEF_avg{isubject},HEF_avg_low{isubject},HEF_avg_high{isubject},rating_ibi_avg(isubject,:)] = WANDER_Heart(isubject,0,timing,rootpath);
    IBI{isubject}       = HEF_avg{isubject}.trialinfo;
    
    % add trialcount discounting blocknr
    x = abs(diff(IBI{isubject}(:,18)));
    indx = find(x > 10);
    for i = 1 : length(IBI{isubject})
        m = find(i<= indx,1,'first');
        if isempty(m) 
            m = 4;
        end
        IBI{isubject}(i,19) = (m-1)*50 + IBI{isubject}(i,18);
    end

    for itrial = 1 : 200
        trial_indx = find(IBI{isubject}(:,19) == itrial);
        if isempty(trial_indx)
            fprintf('Could not find trial %d in subject %d \n',itrial,isubject);
        else
        IBI_trial_avg{isubject}(itrial)   = nanmean(IBI{isubject}(trial_indx,17));
        trialinfo_new{isubject}(itrial,:) = IBI{isubject}(trial_indx(1),:);
        end
    end
        
    F = ceil(2 * tiedrank(IBI{isubject}(:,2)) / length(IBI{isubject}(:,2)));
    IBI_low_avg(isubject)       = nanmean(IBI{isubject}(F == 1,17));   
    IBI_high_avg(isubject)      = nanmean(IBI{isubject}(F == 2,17));
    IBI_low_std(isubject)       = nanstd(IBI{isubject}(F == 1,17));   
    IBI_high_std(isubject)      = nanstd(IBI{isubject}(F == 2,17));    
    IBI_diff_avg(isubject)      = (IBI_low_avg(isubject) - IBI_high_avg(isubject)) / (IBI_low_avg(isubject) + IBI_high_avg(isubject)) ;
    IBI_diff_std(isubject)      = (IBI_low_std(isubject) - IBI_high_std(isubject)) / (IBI_low_std(isubject) + IBI_high_std(isubject));
    IBI_low_count(isubject)     = size(find(F == 1),1);
    IBI_high_count(isubject)    = size(find(F == 2),1);
end


fig = figure; 
subplot(2,2,1);
errorbar([IBI_high_avg(slist)', IBI_low_avg(slist)'],[IBI_high_std(slist)', IBI_low_std(slist)']); axis tight;
title('IBI per condition');
xlabel('Subject');
legend({'Low','High'});

subplot(2,2,2);
errorbar([mean(IBI_high_avg(slist)), mean(IBI_low_avg(slist))],[mean(IBI_high_std(slist)), mean(IBI_low_std(slist))],'o'); xlim([0,3]); 
title('IBI per condition');
xlabel('Attention');

subplot(2,2,3);
errorbar([IBI_diff_avg(slist)'],[IBI_diff_std(slist)']); axis tight; ax = axis;
line([ax(1), ax(2)],[0 0],'linestyle',':','color','k');
title('Relatuve change IBI');
xlabel('Subject');

subplot(2,2,4);
errorbar([mean(IBI_diff_avg(slist))],[mean(IBI_diff_std(slist))],'o'); xlim([0,2]); ax = axis;
line([ax(1), ax(2)],[0 0],'linestyle',':','color','k');
title('Relative change IBI');
xlabel('Attention');

set(fig,'PaperSize',[11*2 8.5*2 ]);
set(fig,'PaperPosition',[0 0 11*4 8.5*2]);
set(fig,'PaperOrientation','landscape');
set(fig,'Position',[50 50 1200 800]);
print(fig,'-dpdf',['d:\analysis\WANDER\images\behaviour\IBI.pdf']);
print(fig,'-dpng',['d:\analysis\WANDER\images\behaviour\IBI.png']);




    % Only correct rejections
    cfg                                = [];
    cfg.trials                         = find(TFR{isubject}.trialinfo(:,3) == 4);
    TFR{isubject}                      = ft_selectdata(cfg,TFR{isubject});
    
    % statistics
    load('D:\analysis\WANDER\scripts\grad','grad');
    neigh_cmb = load('neuromag306cmb_neighb');
    neigh_mag = load('neuromag306mag_neighb_last');
    
    Nobs                    = size(TFR{isubject});
    cfg                     = [];
    cfg.channel             = 'MEG*1';
    % posterior sensors
    % cfg.channel = {'MEG0732+0733', 'MEG0742+0743', 'MEG1522+1523', 'MEG1532+1533', 'MEG1542+1543', 'MEG1612+1613', 'MEG1622+1623', 'MEG1632+1633', 'MEG1642+1643', 'MEG1712+1713', 'MEG1722+1723', 'MEG1732+1733', 'MEG1742+1743', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843', 'MEG1912+1913', 'MEG1922+1923', 'MEG1932+1933', 'MEG1942+1943', 'MEG2012+2013', 'MEG2022+2023', 'MEG2032+2033', 'MEG2042+2043', 'MEG2112+2113', 'MEG2122+2123', 'MEG2132+2133', 'MEG2142+2143', 'MEG2212+2213', 'MEG2222+2223', 'MEG2232+2233', 'MEG2242+2243', 'MEG2312+2313', 'MEG2322+2323', 'MEG2332+2333', 'MEG2342+2343', 'MEG2412+2413', 'MEG2422+2423', 'MEG2432+2433', 'MEG2442+2443', 'MEG2512+2513', 'MEG2522+2523', 'MEG2532+2533', 'MEG2542+2543', 'MEG2622+2623', 'MEG2632+2633', 'MEG2642+2643'};
    % cfg.channel = {'MEG0132+0133', 'MEG0142+0143', 'MEG0212+0213', 'MEG0222+0223', 'MEG0232+0233', 'MEG0242+0243', 'MEG0412+0413', 'MEG0422+0423', 'MEG0432+0433', 'MEG0442+0443', 'MEG0632+0633', 'MEG0642+0643', 'MEG0712+0713', 'MEG0722+0723', 'MEG0732+0733', 'MEG0742+0743', 'MEG1042+1043', 'MEG1112+1113', 'MEG1122+1123', 'MEG1132+1133', 'MEG1142+1143', 'MEG1312+1313', 'MEG1322+1323', 'MEG1332+1333', 'MEG1342+1343', 'MEG1432+1433', 'MEG1442+1443', 'MEG1512+1513', 'MEG1522+1523', 'MEG1532+1533', 'MEG1542+1543', 'MEG1612+1613', 'MEG1622+1623', 'MEG1632+1633', 'MEG1642+1643', 'MEG1712+1713', 'MEG1722+1723', 'MEG1732+1733', 'MEG1742+1743', 'MEG1812+1813', 'MEG1822+1823', 'MEG1832+1833', 'MEG1842+1843', 'MEG1912+1913', 'MEG1922+1923', 'MEG1932+1933', 'MEG1942+1943', 'MEG2012+2013', 'MEG2022+2023', 'MEG2032+2033', 'MEG2042+2043', 'MEG2112+2113', 'MEG2122+2123', 'MEG2132+2133', 'MEG2142+2143', 'MEG2212+2213', 'MEG2222+2223', 'MEG2232+2233', 'MEG2242+2243', 'MEG2312+2313', 'MEG2322+2323', 'MEG2332+2333', 'MEG2342+2343', 'MEG2412+2413', 'MEG2422+2423', 'MEG2432+2433', 'MEG2442+2443', 'MEG2512+2513', 'MEG2522+2523', 'MEG2532+2533', 'MEG2542+2543', 'MEG2612+2613', 'MEG2622+2623', 'MEG2632+2633', 'MEG2642+2643'};
    % cfg.channel = {'MEG1522+1523', 'MEG1532+1533', 'MEG1542+1543', 'MEG1612+1613', 'MEG1622+1623', 'MEG1632+1633', 'MEG1642+1643', 'MEG1712+1713', 'MEG1722+1723', 'MEG1732+1733', 'MEG1742+1743', 'MEG1832+1833', 'MEG1842+1843', 'MEG1912+1913', 'MEG1922+1923', 'MEG1932+1933', 'MEG1942+1943', 'MEG2012+2013', 'MEG2022+2023', 'MEG2032+2033', 'MEG2042+2043', 'MEG2112+2113', 'MEG2122+2123', 'MEG2132+2133', 'MEG2142+2143', 'MEG2232+2233', 'MEG2242+2243', 'MEG2312+2313', 'MEG2322+2323', 'MEG2332+2333', 'MEG2342+2343', 'MEG2412+2413', 'MEG2422+2423', 'MEG2432+2433', 'MEG2442+2443', 'MEG2512+2513', 'MEG2522+2523', 'MEG2532+2533', 'MEG2542+2543', 'MEG2622+2623', 'MEG2632+2633', 'MEG2642+2643'};
    % cfg.channel = {'MEG0731', 'MEG0741', 'MEG1521', 'MEG1531', 'MEG1541', 'MEG1611', 'MEG1621', 'MEG1631', 'MEG1641', 'MEG1711', 'MEG1721', 'MEG1731', 'MEG1741', 'MEG1811', 'MEG1821', 'MEG1831', 'MEG1841', 'MEG1911', 'MEG1921', 'MEG1931', 'MEG1941', 'MEG2011', 'MEG2021', 'MEG2031', 'MEG2041', 'MEG2111', 'MEG2121', 'MEG2131', 'MEG2141', 'MEG2211', 'MEG2221', 'MEG2231', 'MEG2241', 'MEG2311', 'MEG2321', 'MEG2331', 'MEG2341', 'MEG2411', 'MEG2421', 'MEG2431', 'MEG2441', 'MEG2511', 'MEG2521', 'MEG2531', 'MEG2541', 'MEG2621', 'MEG2631', 'MEG2641'};
    cfg.latency             = [1 9];
    cfg.avgovertime         = 'yes';
    cfg.avgoverchan         = 'no';
    cfg.parameter           = 'powspctrm';
    cfg.method              = 'analytic';
    cfg.statistic           = 'ft_statfun_correlationT';
    % cfg.correctm            = 'cluster';
    % cfg.clusteralpha        = 0.005;
    % cfg.clusterstatistic    = 'maxsum';
    % cfg.minnbchan           = 2;
    cfg.tail                = 0;
    % cfg.clustertail         = 0;
    cfg.alpha               = 0.05;
    % cfg.numrandomization    = 2000;
    cfg.neighbours          = neigh_mag.neighbours;
    % cfg.neighbours          = neigh_mag.neighbours;
    cfg.design(1,:)         = IBI_trial_avg{isubject};
    cfg.ivar                = 1;
    stat{isubject}          = ft_freqstatistics(cfg,TFR{isubject});
    
    
    % plot topos MI values
    cfg             = [];
    cfg.layout      = 'neuromag306mag';
    cfg.parameter   = 'stat';
    cfg.channel     = 'MEG*1';
    cfg.interactive = 'yes';
    cfg.zlim        = 'maxabs';
    % cfg.ylim        = [10 11];
    figure;  ft_topoplotER(cfg,stat{isubject});
end
