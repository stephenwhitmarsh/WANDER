function launch_WANDER_batch


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

WANDER_timelock_SSEF_noEGG(isubject,rootpath,force)