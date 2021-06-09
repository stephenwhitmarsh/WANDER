function [trl, event] = WANDER_trialfun_restingstate(cfg)

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = ft_read_event(cfg.dataset);

% search for "trigger" events
value  = [event(find(strcmp('STI101', {event.type}))).value]';
sample = [event(find(strcmp('STI101', {event.type}))).sample]';

if isempty(find(value==1))
    onset = 1;
else
    onset = sample(value == 1);
end

if hdr.nSamples < 12.5*60*hdr.Fs
    offset = hdr.nSamples;
else
    offset = 12.5*60*hdr.Fs;
end

trl = [onset offset 0];


