
function [data_probelocked_MEG data_probelocked_EGG] = WANDER_redefine_to_probe(isubject,force,rootpath)

if isempty(rootpath)
    rootpath = 'd:\analysis\WANDER\data\';
end

fname_probelocked = [rootpath filesep 'trial\s' num2str(isubject) '_epoched_MEG_probelocked.mat'];

if exist(fname_probelocked,'file') && force ~= 1
    fprintf('Returning probelocked MEG data\n');
    load(fname_probelocked);
else
    fprintf('Probelocked MEG data not found, creating it now! \n');
    
    temp = WANDER_ICA_round2(isubject,0);
    data_ICA = temp.data_ICA;
    clear temp
    
    data_EGG_epoched = WANDER_epoch_EGG(isubject,0);

    for ipart = 1 : 4
        for itrial = 1 : size(data_ICA{ipart}.trial,2)
            
            fprintf('Resizing and ordering trial %d of block %d\n',itrial,ipart);
            
            % cut off 1.5 seconds at the beginning,
            data_ICA{ipart}.time{itrial}                    = data_ICA{ipart}.time{itrial}(:,1501:end);
            data_ICA{ipart}.trial{itrial}                   = data_ICA{ipart}.trial{itrial}(:,1501:end);
            data_ICA{ipart}.sampleinfo(itrial,1)            = data_ICA{ipart}.sampleinfo(itrial,1) + 1500;
            
            data_EGG_epoched{ipart}.time{itrial}            = data_EGG_epoched{ipart}.time{itrial}(:,1501:end);
            data_EGG_epoched{ipart}.trial{itrial}           = data_EGG_epoched{ipart}.trial{itrial}(:,1501:end);
            data_EGG_epoched{ipart}.sampleinfo(itrial,1)    = data_EGG_epoched{ipart}.sampleinfo(itrial,1) + 1500;
            
            % cut off 1.0 seconds at the end
            data_ICA{ipart}.time{itrial}                    = data_ICA{ipart}.time{itrial}(:,1:end-1000);
            data_ICA{ipart}.trial{itrial}                   = data_ICA{ipart}.trial{itrial}(:,1:end-1000);
            data_ICA{ipart}.sampleinfo(itrial,2)            = data_ICA{ipart}.sampleinfo(itrial,2) - 1000;
            
            data_EGG_epoched{ipart}.time{itrial}            = data_EGG_epoched{ipart}.time{itrial}(:,1:end-1000);
            data_EGG_epoched{ipart}.trial{itrial}           = data_EGG_epoched{ipart}.trial{itrial}(:,1:end-1000);
            data_EGG_epoched{ipart}.sampleinfo(itrial,2)    = data_EGG_epoched{ipart}.sampleinfo(itrial,2) - 1000;
            
            % reorder time axis to end
            end_time = data_ICA{ipart}.time{itrial}(end);
            data_ICA{ipart}.time{itrial}  = data_ICA{ipart}.time{itrial} - end_time;
            data_EGG_epoched{ipart}.time{itrial}  = data_EGG_epoched{ipart}.time{itrial} - data_EGG_epoched{ipart}.time{itrial}(end);
            
            
        end
    end
    
    fprintf('Saving EGG, MEG and artefact data. This can take a while.\n');
    data_probelocked_MEG = data_ICA;
    data_probelocked_EGG = data_EGG_epoched;
    save(fname_probelocked,'data_probelocked_MEG','data_EGG_epoched','-v7.3');
 
end