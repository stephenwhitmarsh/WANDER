function [trl, event] = WANDER_trialfun_cuelocked(cfg)

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% search for "trigger" events
value  = [event(find(strcmp('STI101', {event.type}))).value]';
sample = [event(find(strcmp('STI101', {event.type}))).sample]';

% remove extra triggers
% for n = 16:-1:8
%     value(value>=2^n) = value(value>=2^n) - 2^n;
% end

stim_pre            = -cfg.trialdef.stim_pre  * hdr.Fs;
stim_post           =  cfg.trialdef.stim_post * hdr.Fs;
stim_indx           = find(ismember(value,cfg.trialdef.stim));
onset_indx          = find(ismember(value,cfg.trialdef.onset));
rating_indx         = find(ismember(value,cfg.trialdef.rating));
offset_indx         = find(ismember(value,cfg.trialdef.offset));
performance_indx    = find(ismember(value,cfg.trialdef.performance));

% to fix bug in subject 17: in part 2, trial 14 & in part 3, trial 41
d = find(diff(rating_indx)==1);
if ~isempty(d)
    rating_indx(d+1) = [];
end

trl  = [];
for itrial = 1:length(onset_indx)
    if itrial > 1 || onset_indx(1) < offset_indx(1) % in case file starts in middle of trial
        try
            last_stim    = find(stim_indx < offset_indx(itrial),1,'last');
            rating       = value(rating_indx(itrial)) - 10;
            rating_RT    = (sample(rating_indx(itrial)) - sample(onset_indx(itrial))) / hdr.Fs;
            performance  = value(performance_indx(itrial)) - 20;
            stimduration = (sample(offset_indx(itrial)) - sample(onset_indx(itrial))) / hdr.Fs;
            trl_new = [sample(onset_indx(itrial))+stim_pre sample(stim_indx(last_stim))+stim_post stim_pre value(onset_indx(itrial)) rating performance stimduration rating_RT itrial];
            trl = [trl; trl_new];
        catch
        end
    end
end;

