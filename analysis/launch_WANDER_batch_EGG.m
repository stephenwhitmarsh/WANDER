function launch_WANDER_batch

addpath('/home/swhitmarsh/WANDER/scripts/');
addpath('/home/swhitmarsh/WANDER/fieldtrip/');

jobID   = getenv('SLURM_JOB_ID');
arrayID = getenv('SLURM_ARRAY_JOB_ID');
taskID  = getenv('SLURM_ARRAY_TASK_ID');

fprintf('jobID %d arrayID %d taskID %d\n',jobID, arrayID,str2num(taskID));
isubject = str2num(taskID);

force = 1;
timing = 'cue';
latency = 'all';
rootpath = 0;
restingstate = 1;
WANDER_filter_EGG_param(isubject,force,rootpath,restingstate)
