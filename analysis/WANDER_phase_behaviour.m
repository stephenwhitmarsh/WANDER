function WANDER_phase_ratings

addpath d:/fieldtrip/
addpath('D:/analysis/WANDER/scripts/');
addpath('D:/analysis/WANDER/scripts/CircStat2012a/');
addpath('D:/analysis/WANDER/scripts/PhaseOppositionCode/');
ft_defaults

slist = [1:5 8:13 15:20 22:26]; %without subjects with more than 2 SD

clear rating_* phase_end phase_begin phase_rand rating

for isubject = [1:26]
    
    % prepare data
    data_EGG = WANDER_filter_EGG(isubject,0,rootpath,0);
    for iblock = 1 : 4
        cfg = [];
        cfg.channel = [FFT.max_chan 'phase'];
        cfg.trials  = data_EGG{iblock}.trialinfo(:,3) == 4;
        data_EGG{iblock} = ft_selectdata(cfg,data_EGG{iblock});
    end
    
    i = 1;
    for iblock = 1 : 4
        for itrial = 1 : 50
            try % since we have less that 50 trials after rejecting incorrect trials
                phase_end{isubject}(i)    = data_EGG{iblock}.trial{itrial}(1,end-1000);
                phase_begin{isubject}(i)  = data_EGG{iblock}.trial{itrial}(1,1500);
                phase_rand{isubject}(i)   = data_EGG{iblock}.trial{itrial}(1,randi([1500,size(data_EGG{iblock}.trial{itrial},2)]));                
                rating{isubject}(i)       = data_EGG{iblock}.trialinfo(itrial,2);
                i = i + 1;
            catch
            end
        end
    end
    
    rating_rand{isubject} = rating{isubject}(randperm(size(rating{isubject},2)));
    
    % prepare phase bins
    phase_nrbins = 21;                                               % The number of phase bins
    phase_bins   = -pi:2*pi/phase_nrbins:pi;                        % The extreme phases of each bin
    phase_axis   = (phase_bins(1:end-1) + phase_bins(2:end)) / 2;   % The midpoint of each phase bin
    
    % bin rating in phase bins
    clear phase_indx edges phase_indx
    [phase_count_end{isubject},           edges,phase_indx_end{isubject}]          = histcounts(phase_end{isubject},   phase_bins);  % bin phase in bins, check usage e.g. with: [count, edges, index] = histcounts([1 2 3 1 2 3 6 7 8 4 2 1],[1 2 3]);
    [phase_count_begin{isubject},         edges,phase_indx_begin{isubject}]        = histcounts(phase_begin{isubject}, phase_bins);  % bin phase in bins, check usage e.g. with: [count, edges, index] = histcounts([1 2 3 1 2 3 6 7 8 4 2 1],[1 2 3]);
    [phase_count_rand{isubject},          edges,phase_indx_rand{isubject}]         = histcounts(phase_rand{isubject},  phase_bins);  % bin phase in bins, check usage e.g. with: [count, edges, index] = histcounts([1 2 3 1 2 3 6 7 8 4 2 1],[1 2 3]);
    [phase_count_rating_rand{isubject},   edges,phase_indx_rating_rand{isubject}]  = histcounts(phase_end{isubject},   phase_bins);  % bin phase in bins, check usage e.g. with: [count, edges, index] = histcounts([1 2 3 1 2 3 6 7 8 4 2 1],[1 2 3]);
   
    % to avoid confound of averaging 0s
    phase_indx_end{isubject}(phase_count_end{isubject} == 0)                 = NaN;
    phase_indx_begin{isubject}(phase_count_begin{isubject} == 0)             = NaN;
    phase_indx_rand{isubject}(phase_count_rand{isubject} == 0)               = NaN;
    phase_indx_rating_rand{isubject}(phase_count_rating_rand{isubject} == 0) = NaN;

    % calculations per bin
    for ibin = unique(phase_indx_end{isubject})
        rating_bincount_end{isubject}(ibin)  = size(find(phase_indx_end{isubject} == ibin),2);
        rating_binmean_end{isubject}(ibin)   = nanmean(rating{isubject}(phase_indx_end{isubject} == ibin));
        rating_binmedian_end{isubject}(ibin) = nanmedian(rating{isubject}(phase_indx_end{isubject} == ibin));
        rating_binstd_end{isubject}(ibin)    = nanstd(rating{isubject}(phase_indx_end{isubject} == ibin));
    end
    
    % calculations per bin
    for ibin = unique(phase_indx_begin{isubject})
        rating_bincount_begin{isubject}(ibin)  = size(find(phase_indx_begin{isubject} == ibin),2);
        rating_binmean_begin{isubject}(ibin)   = nanmean(rating{isubject}(phase_indx_begin{isubject} == ibin));
        rating_binmedian_begin{isubject}(ibin) = nanmedian(rating{isubject}(phase_indx_begin{isubject} == ibin));
        rating_binstd_begin{isubject}(ibin)    = nanstd(rating{isubject}(phase_indx_begin{isubject} == ibin));
    end
        
    % calculations per bin
    for ibin = unique(phase_indx_rand{isubject})
        rating_bincount_rand{isubject}(ibin)  = size(find(phase_indx_rand{isubject} == ibin),2);
        rating_binmean_rand{isubject}(ibin)   = nanmean(rating{isubject}(phase_indx_rand{isubject} == ibin));
        rating_binmedian_rand{isubject}(ibin) = nanmedian(rating{isubject}(phase_indx_rand{isubject} == ibin));
        rating_binstd_rand{isubject}(ibin)    = nanstd(rating{isubject}(phase_indx_rand{isubject} == ibin));
    end
    
    % calculations per bin
    for ibin = unique(phase_indx_rating_rand{isubject})
        rating_bincount_rating_rand{isubject}(ibin)  = size(find(phase_indx_rating_rand{isubject} == ibin),2);
        rating_binmean_rating_rand{isubject}(ibin)   = nanmean(rating_rand{isubject}(phase_indx_rating_rand{isubject} == ibin));
        rating_binmedian_rating_rand{isubject}(ibin) = nanmedian(rating_rand{isubject}(phase_indx_rating_rand{isubject} == ibin));
        rating_binstd_rating_rand{isubject}(ibin)    = nanstd(rating_rand{isubject}(phase_indx_rating_rand{isubject} == ibin));
    end
    
    % shift to max value
    [~,max_indx] = max(rating_binmean_end{isubject});
    shift = ceil(phase_nrbins/2) - max_indx;
    
    % shift = 0;
    rating_bincount_end_shift{isubject}             = circshift(rating_bincount_end{isubject},[0,shift]);
    rating_binmean_end_shift{isubject}              = circshift(rating_binmean_end{isubject},[0,shift]);
    rating_binmedian_end_shift{isubject}            = circshift(rating_binmedian_end{isubject},[0,shift]);
    rating_binstd_end_shift{isubject}               = circshift(rating_binstd_end{isubject},[0,shift]);
    
    [~,max_indx] = max(rating_binmean_begin{isubject});
    shift = ceil(phase_nrbins/2) - max_indx;
    % shift = 0;
    rating_bincount_begin_shift{isubject}           = circshift(rating_bincount_begin{isubject},[0,shift]);
    rating_binmean_begin_shift{isubject}            = circshift(rating_binmean_begin{isubject},[0,shift]);
    rating_binmedian_begin_shift{isubject}          = circshift(rating_binmedian_begin{isubject},[0,shift]);
    rating_binstd_begin_shift{isubject}             = circshift(rating_binstd_begin{isubject},[0,shift]);
    
    [~,max_indx] = max(rating_binmean_rand{isubject});
    shift = ceil(phase_nrbins/2) - max_indx;
    % shift = 0;
    rating_bincount_rand_shift{isubject}            = circshift(rating_bincount_rand{isubject},[0,shift]);
    rating_binmean_rand_shift{isubject}             = circshift(rating_binmean_rand{isubject},[0,shift]);
    rating_binmedian_rand_shift{isubject}           = circshift(rating_binmedian_rand{isubject},[0,shift]);
    rating_binstd_rand_shift{isubject}              = circshift(rating_binstd_rand{isubject},[0,shift]);
    
    [~,max_indx] = max(rating_binmean_rating_rand{isubject});
    shift = ceil(phase_nrbins/2) - max_indx;
    % shift = 0;
    rating_bincount_rating_rand_shift{isubject}     = circshift(rating_bincount_rating_rand{isubject},[0,shift]);
    rating_binmean_rating_rand_shift{isubject}      = circshift(rating_binmean_rating_rand{isubject},[0,shift]);
    rating_binmedian_rating_rand_shift{isubject}    = circshift(rating_binmedian_rating_rand{isubject},[0,shift]);
    rating_binstd_rating_rand_shift{isubject}       = circshift(rating_binstd_rating_rand{isubject},[0,shift]);

    close all
end
clear EGG

% Plot phase/power distribution with error bars
figure; hold
for isubject = [1:26] 
    try
        subplot(5,6,isubject); hold;
        errorbar(phase_axis, rating_binmean_end_shift{isubject}, rating_binstd_end_shift{isubject}./sqrt(rating_bincount_end_shift{isubject}));
        axis tight
    catch
    end
end

% Plot phase/power distribution with error bars - with phase locked to end,
% beginning or random period in trial
figure;
for isubject = [1:26]
    subplot(5,6,isubject); hold;
    errorbar(phase_axis, rating_binmean_end{isubject},      rating_binstd_end{isubject}./sqrt(rating_bincount_end{isubject}));
    errorbar(phase_axis, rating_binmean_begin{isubject},    rating_binstd_begin{isubject}./sqrt(rating_bincount_begin{isubject}));
    errorbar(phase_axis, rating_binmean_rand{isubject},     rating_binstd_rand{isubject}./sqrt(rating_bincount_rand{isubject}));
    axis tight
end

% investigate circular distribution
for isubject = 1:26
    m                                               = quantile(rating{isubject},0.5);
    F                                               = ceil(2 * tiedrank(rating{isubject}) / length(rating{isubject}));
    rating_split{isubject}                          = false(size(F));
    rating_split{isubject}(F==2)                    = true;
    circ{isubject}                                  = circ_stats(phase_end{isubject});
    circ_stat_h{isubject}                           = circ_stats(phase_end{isubject}(rating_split{isubject}));
    circ_stat_l{isubject}                           = circ_stats(phase_end{isubject}(~rating_split{isubject}));
    circ_r_h(isubject)                              = circ_r(phase_end{isubject}(rating_split{isubject})');
    circ_r_l(isubject)                              = circ_r(phase_end{isubject}(~rating_split{isubject})');
end

figure;
hist([circ_r_h' circ_r_l']);
mean([circ_r_h' circ_r_l'])

for isubject = 1:26
    [pval, z] = circ_rtest(rating_binmean_end{isubject},rating_bincount_end{isubject},phase_axis)
    [pval, z] = circ_rtest(rating_binmean_begin{isubject},rating_bincount_begin{isubject},phase_axis)
    [pval, z] = circ_rtest(rating_binmean_rand{isubject},rating_bincount_rand{isubject},phase_axis) 
end

save everything

% Use VanRullen 2016
for isubject = [1:26]
    [p_circWW(isubject), p_POS(isubject), p_zPOS(isubject)] = PhaseOpposition(phase_end{isubject}(rating_split{~isubject}),phase_end{isubject}(rating_split{isubject}),50000);
end
result = combine_pvalues(p_POS,2,1)

% Use circ toolbox
%
%   [rho pval ts] = circ_corrcl(alpha, x)
%   Correlation coefficient between one circular and one linear random
%   variable.
%
%   Input:
%     alpha   sample of angles in radians
%     x       sample of linear random variable
%
%   Output:
%     rho     correlation coefficient
%     pval    p-value

for isubject = [1:26]
    [rho(isubject) pval(isubject)] = circ_corrcl(phase_end{isubject}, rating{isubject});
end

figure;
for isubject = [1:26]
    subplot(5,6,isubject);
    % polarhistogram('BinEdges',phase_bins,'BinCounts',rating_bincount_end{isubject},'FaceAlpha',0.5);
    % polarhistogram(phase_end{isubject}(rating_split{isubject}),'normalization','count','FaceAlpha',0.5,'FaceColor',[0 1 0]); hold;
    % polarhistogram(phase_end{isubject}(~rating_split{isubject}),'normalization','count','FaceAlpha',0.5,'FaceColor',[1 0 0]);
%     polarhistogram(phase_end{isubject}(rating_split{isubject}),'normalization','count','DisplayStyle','stairs'); hold; axis tight
%     polarhistogram(phase_end{isubject}(~rating_split{isubject}),'normalization','count','DisplayStyle','stairs'); axis tight
%     circ_plot(phase_end{isubject}(rating_split{isubject})','pretty','ro',true,'linewidth',2,'color','r');
    circ_plot(phase_rand{isubject}(rating_split{isubject})','hist',[],30,true,true,'linewidth',2,'color','r'); hold;
end

figure;
for isubject = [1:26]
    subplot(5,6,isubject);
    circ_plot(phase_end{isubject}(rating_split{~isubject})','hist',[],30,true,true,'linewidth',2,'color','r'); hold;
end



