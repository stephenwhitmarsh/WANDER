function launch_WANDER_beamformer_batch


addpath /shared/projects/project_wander/WANDER/scripts/
addpath /shared/projects/project_wander/WANDER/fieldtrip/
ft_defaults

jobID   = getenv('SLURM_JOB_ID');
arrayID = getenv('SLURM_ARRAY_JOB_ID');
taskID  = getenv('SLURM_ARRAY_TASK_ID');

fprintf('jobID %d arrayID %d taskID %d\n',jobID, arrayID,str2num(taskID));
isubject = str2num(taskID);

force = 1;
rootpath = 0;    

% Calculate common filter for individual alpha freq
% WANDER_common_filter_DICS_individual_onlyalpha(isubject,rootpath,force)

% Source analysis median split
WANDER_source_individual_onlyalpha_GRAD(isubject,force,rootpath)
