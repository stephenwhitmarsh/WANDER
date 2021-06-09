
function [data_EGG_epoched] = WANDER_redefine_EGG_to_probe(isubject,force,rootpath)

if isempty(rootpath)
    rootpath = 'd:\analysis\WANDER\data\';
end

fname_probelocked = [rootpath filesep 'trial\s' num2str(isubject) '_epoched_EGG_probelocked.mat'];

if exist(fname_probelocked,'file') && force ~= 1
    fprintf('Returning probelocked EGG data\n');
    load(fname_probelocked);
else
    fprintf('Probelocked EGG data not found, creating it now! \n');
    data_EGG_epoched = WANDER_epoch_EGG(isubject,0);
    
    for ipart = 1 : 4
        for itrial = 1 : size(data_EGG_epoched{ipart}.trial,2)
            
            fprintf('Resizing and ordering EGG trial %d of block %d\n',itrial,ipart);
            
            % cut off 1.5 seconds at the beginning,
            data_EGG_epoched{ipart}.time{itrial}            = data_EGG_epoched{ipart}.time{itrial}(:,1501:end);
            data_EGG_epoched{ipart}.trial{itrial}           = data_EGG_epoched{ipart}.trial{itrial}(:,1501:end);
            data_EGG_epoched{ipart}.sampleinfo(itrial,1)    = data_EGG_epoched{ipart}.sampleinfo(itrial,1) + 1500;
            
            % cut off 1.0 seconds at the end
            data_EGG_epoched{ipart}.time{itrial}            = data_EGG_epoched{ipart}.time{itrial}(:,1:end-1000);
            data_EGG_epoched{ipart}.trial{itrial}           = data_EGG_epoched{ipart}.trial{itrial}(:,1:end-1000);
            data_EGG_epoched{ipart}.sampleinfo(itrial,2)    = data_EGG_epoched{ipart}.sampleinfo(itrial,2) - 1000;
            
            % reorder time axis to end
            data_EGG_epoched{ipart}.time{itrial}  = data_EGG_epoched{ipart}.time{itrial} - data_EGG_epoched{ipart}.time{itrial}(end);
        end
    end
    
    fprintf('Saving data. This can take a while.\n');
    save(fname_probelocked,'data_EGG_epoched','-v7.3');
end